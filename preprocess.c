#include <stdio.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include "macros.h"
#include "memory.h"
#include "lp.h"
#include "cprop.h"
#include "containers.h"
#include "strutils.h"

int main( int argc, char **argv )
{
    if (argc<2)
    {
        fprintf( stderr, "Enter instance name.\n" );
        exit( EXIT_FAILURE );
    }

    LinearProgram *mip = lp_create();
    clock_t cla = clock();
    lp_read( mip, argv[1] );
    clock_t clb = clock();
    double secread = ((double)clb-cla) / ((double)CLOCKS_PER_SEC);

    printf("Read in %.4f seconds\n", secread );
    
    char feasible = 1;
    clock_t stopt = clock();
    int status = lp_optimize_as_continuous( mip );
    clock_t edopt = clock();
    double secopt = ((double)edopt-stopt) /  ((double)CLOCKS_PER_SEC);

    double obj = DBL_MAX;
    if (status!=LP_OPTIMAL)
    {
        feasible = 0;

        if (status!=LP_INFEASIBLE)
        {
            fprintf( stderr, "Error while optimizing. Status %d\n", status );
            exit( EXIT_FAILURE );
        }

    }
    else
    {
        feasible = 1;
        obj = lp_obj_value( mip );
        printf("Initial LP relaxation solved in %.4f seconds, obj is %g\n", secopt, lp_obj_value(mip) );
    } 
    
    int n = lp_cols( mip );
    
    char *integer;
    ALLOCATE_VECTOR( integer, char, n );
    double *lb, *ub, *coef;
    int *idx;
    ALLOCATE_VECTOR( idx, int, n );
    ALLOCATE_VECTOR( lb, double, n );
    ALLOCATE_VECTOR( ub, double, n );
    ALLOCATE_VECTOR( coef, double, n );
        
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
    
    clock_t cpst = clock();

    CProp *cprop = cprop_create( n, integer, lb, ub, (const char**)strv_ptr( names ) );

    cprop_set_verbose( cprop, 1 );
        
    clock_t cped;
    double cpsecs;
    char madeInf;
    double obj2;
    int nImpl;
    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
        if (!cprop_feasible(cprop))
        {
            char rname[256];
            lp_row_name( mip, i, rname );
            printf("Adding constraints %s made the problem infeasible.\n", rname );
            printf("CPROP msg: %s\n", cprop_inf_msg(cprop) );
            goto END;
        }
    }
    
    cprop_conclude_pre_processing( cprop );
    cped = clock();
    cpsecs = ((double)cped-cpst) / ((double)CLOCKS_PER_SEC);
    printf("CPROP in %.4f seconds.\n", cpsecs );

    nImpl = 0;
    for ( int j=0 ; (j<lp_cols(mip)) ; ++j )
    {
        if (cprop_get_lb(cprop,j)>=lp_col_lb(mip,j)+1e-10 || cprop_get_ub(cprop,j)<=lp_col_ub(mip,j)-1e-10)
        {
            lp_set_col_bounds( mip, j, cprop_get_lb(cprop,j), cprop_get_ub(cprop,j));
            nImpl++;
        }
    }
    printf("CPROP fixed %d variables.\n", nImpl );

    madeInf = 0;
    obj2 = DBL_MAX;
    if (nImpl)
    {
        int st2 = lp_optimize_as_continuous( mip );
        if (feasible)
            madeInf = (st2 == LP_INFEASIBLE);
        if (st2==LP_OPTIMAL)
            obj2 = lp_obj_value( mip );

        lp_write_lp( mip, "pp.lp" );

    }
    else
        obj2 = obj;
    

    // writing summary
    {
        char probName[256];
        getFileName( probName, argv[1] );

        char exists = 0;

        // checking if file exists
        {
            FILE *fe = fopen("summary.csv", "r");
            exists = ( fe != NULL );
            if (exists)
                fclose(fe);
        }

        FILE *f = fopen("summary.csv", "a");
        if (!exists)
            fprintf( f, "instance,nImpl,time,madeInf,obj1,obj2\n" );

        fprintf( f, "%s,%d,%.4f,%d,%g,%g\n", probName, nImpl, cpsecs, madeInf, obj, obj2 );

        fclose(f);
    } 
END:
    cprop_free( &cprop );
    lp_free( &mip );

}

