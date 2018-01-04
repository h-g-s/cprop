#include <assert.h>
#include <math.h>
#include <string.h>
#include "cprop_lp.h"
#include "macros.h"
#include "containers.h"

CProp *cprop_create_from_mip( LinearProgram *mip, char verbose )
{
    int n = lp_cols( mip );

    char *integer;
    ALLOCATE_VECTOR( integer, char, n );
    double *lb, *ub, *coef;
    int *idx;
    ALLOCATE_VECTOR( idx, int, n );
    ALLOCATE_VECTOR( lb, double, 3*n );
    ub = lb + n;
    coef = ub + n;

    for ( int i=0 ; (i<n) ; ++i )
        integer[i] = lp_is_integer( mip, i );
    for ( int i=0 ; (i<n) ; ++i )
        lb[i] = lp_col_lb( mip, i );
    for ( int i=0 ; (i<n) ; ++i )
        ub[i] = lp_col_ub( mip, i );

    StrV *names = strv_create( 256 );
    for ( int i=0 ; (i<n) ; ++i )
    {
        char cname[256];
        strv_push_back( names, lp_col_name( mip, i, cname ) );
    }


    CProp *cprop = cprop_create( n, integer, lb, ub, (const char**)strv_ptr( names ) );

    cprop_set_verbose( cprop, verbose );

    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
        if (!cprop_feasible(cprop))
        {
            if (verbose)
            {
                char rname[256];
                lp_row_name( mip, i, rname );
                printf("Adding constraint %s made the problem infeasible.\n", rname );
                printf("CPROP msg: %s\n", cprop_inf_msg(cprop) );                
            }
            goto END;
        }
    }

    cprop_conclude_pre_processing( cprop );    

END:    
    free( integer );
    free( idx );
    free( lb );
    strv_free( &names );

    return cprop;
}

LinearProgram *cprop_preprocess( const CProp *cprop, LinearProgram *mip, char removeVars, int origCols[] )
{
    if (!cprop_feasible(cprop))
        return NULL;

    int nFixed = 0;
    for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        if (cprop_fixed_at_pre_proc(cprop)[i])
            ++nFixed;

    if (nFixed==0)
        return NULL;

    LinearProgram *rmip = lp_create();

    int nCols = lp_cols(mip) - ( removeVars ? nFixed : 0 );

    double *lb, *ub, *obj;
    ALLOCATE_VECTOR( lb, double, nCols*3 );
    ub = lb+nCols;
    obj = ub+nCols;
    char *integer;
    ALLOCATE_VECTOR( integer, char, nCols );
    StrV *cnames = strv_create( 256 );

    int *ppCol;
    ALLOCATE_VECTOR( ppCol, int, lp_cols(mip) );

    double objconst = 0.0;
    if (removeVars)
    {
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
            ppCol[i] = -1;        

        assert( origCols );
        int icol = 0;
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        {
            if (cprop_fixed_at_pre_proc(cprop)[i])
            {
                if (cprop_get_lb( cprop, i )>=0.99)
                    objconst += cprop_get_lb( cprop, i )*lp_obj_coef(mip)[i];
            }
            else
            {
                lb[icol] =  cprop_get_lb( cprop, i );
                ub[icol] =  cprop_get_ub( cprop, i );
                obj[icol] =  lp_obj_coef(mip)[i];
                integer[icol] = lp_is_integer(mip,i);
                origCols[icol] = i;
                ppCol[i] = icol;
                char name[256];
                strv_push_back( cnames, lp_col_name(mip,i,name) );
                ++icol;
            }
        }
    }
    else
    {
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        {
            ppCol[i] = i;
            lb[i] =  cprop_get_lb( cprop, i );
            ub[i] =  cprop_get_ub( cprop, i );
            obj[i] =  lp_obj_coef(mip)[i];
            integer[i] = lp_is_integer(mip,i);
            origCols[i] = i;
            char name[256];
            strv_push_back( cnames, lp_col_name(mip,i,name) );
        }
    }

    lp_add_cols( rmip, nCols, obj, lb, ub, integer, strv_ptr(cnames) );
    if (fabs(objconst)>1e-15)
    {
        char cname[] = "objConst";
        lp_add_col( rmip, objconst, 1.0, 1.0, 1, cname, 0, NULL, NULL );
    }
    
    int *idx; ALLOCATE_VECTOR( idx, int, lp_cols(rmip) );
    double *coef; ALLOCATE_VECTOR( coef, double, lp_cols(rmip) );

    for ( int i=0 ; (i<cprop_n_rows(cprop)) ; ++i )
    {
        const int *cpidx = cprop_idx( cprop, i );
        int nz = cprop_nz( cprop, i );
        for ( int j=0 ; (j<nz) ; ++j )
            idx[j] = ppCol[cpidx[j]];

        memcpy( coef, cprop_coef(cprop,i), sizeof(double)*nz );
        double rhs = cprop_rhs( cprop, i );

        if (cprop_is_equality(cprop,i))
        {
            char rname[256]; 
            strcpy( rname, cprop_row_name( cprop, i )  );
            lp_add_row( rmip, nz, idx, coef, rname, 'E', rhs );
            ++i;
        }
        else
        {
            char sense = 'L';
            int nneg = 0;
            for ( int j=0 ; (j<nz) ; ++j )
                if (coef[j]<-1e-15)
                    ++nneg;
            if ( nneg > ((int)ceil(((double)nz)/2.0)) )
            {
                sense = 'G';
                rhs *= -1.0;
                for ( int j=0 ; (j<nz) ; ++j )
                    coef[j] *= -1.0;
            }

            lp_add_row( rmip, nz, idx, coef, cprop_row_name( cprop,i ), sense, rhs );
        }
    }

    free( idx );
    free( coef );
    free( lb );
    free( integer );
    strv_free( &cnames );
    free( ppCol );

    return rmip;
}
