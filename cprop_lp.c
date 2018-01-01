#include <assert.h>
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

    if (removeVars)
    {
        assert( origCols );
        int icol = 0;
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        {
            if (!cprop_fixed_at_pre_proc(cprop)[i])
            {
                lb[icol] =  cprop_get_lb( cprop, i );
                ub[icol] =  cprop_get_ub( cprop, i );
                obj[icol] =  lp_obj_coef(mip)[i];
                integer[icol] = lp_is_integer(mip,i);
                origCols[icol] = i;
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
    
    free( lb );
    free( integer );
    strv_free( &cnames );
    
    return rmip;
}
