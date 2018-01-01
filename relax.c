#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cprop.h"
#include "lp.h"
#include "memory.h"
#include "containers.h"

int dive_and_propagate( LinearProgram *mip, CProp *cprop, int maxTries, CPCuts *gcuts );

int main( int argc, char **argv )
{
    if (argc<2)
    {
        fprintf( stderr, "Enter instance name.\n" );
        exit( 1 );
    }

    LinearProgram *mip = lp_create();
    lp_read( mip, argv[1] );

    lp_set_print_messages( mip, 0 );

    // getting variables info
    
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
    
    CProp *cprop = cprop_create( n, integer, lb, ub, (const char**)strv_ptr( names ) );
    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
        if (!cprop_feasible(cprop))
            goto END;
    }
    
    
    // *updating bounds
    cprop_update_bound( cprop, 0, 1.0, 1.0 );
    int status = cprop_update_bound( cprop, 1, 1.0, 1.0 );
    printf("st %d\n", status);
    cprop_save_impl_graph( cprop, "imppl.dot");
    
    

    CPCuts *gcuts = cpc_create( MAX(lp_cols(mip), lp_rows(mip) ) );
    
    int nNewCuts = 1, it = 1;
    while (nNewCuts)
    {
        int status = lp_optimize_as_continuous( mip );
        if ( status == LP_INFEASIBLE )
        {
            fprintf( stderr, "INFEASIBLE LP. Check correctness of cuts.\n" );
            abort();
        }
        if ( status != LP_OPTIMAL )
        {
            fprintf( stderr, "Status not optimal obtained (%d).\n", status );
            abort();
        }
        printf("> iteration %d obj: %g\n", it++, lp_obj_value(mip) );

        nNewCuts= dive_and_propagate( mip, cprop, 3, gcuts );
    }

END:
    lp_free( &mip );
    cprop_free( &cprop );

    exit( 0 );
}

int dive_and_propagate( LinearProgram *mip, CProp *cprop, int maxTries, CPCuts *gcuts )
{
    const double *x = lp_x(mip);

    int nCuts = 0;

    Vec_int *scols = vec_int_create_cap( lp_cols(mip) );
    Vec_IntPair *spairs = vec_IntPair_create_cap( lp_cols(mip) );
    for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        vec_int_push_back( scols, i );

    for ( int t=0 ; (t<maxTries) ; ++t )
    {
        cprop_clear( cprop );
        vec_int_shuffle( scols, spairs );
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        {
            int iv = vec_int_get( scols, i );
            if (!lp_is_binary(mip, iv ))
                continue;

            /*
            if ( x[iv]<=1e-7 || x[iv]>=0.99999 )
                continue;*/

            double r = floor( x[iv]+0.5 );
            
            int status = cprop_update_bound( cprop, iv, r, r );
            if (status==-1)
            {
                assert( cprop_feasible( cprop )==0 );
                const CPCuts *ccCuts = cprop_cut_pool(cprop);
                for ( int j=0 ; (j<cpc_n_cuts(ccCuts)); ++j )
                {
                    double viol = 0;
                    int nz = cpc_nz(ccCuts, j);
                    const int *idx = cpc_idx( ccCuts, j );
                    const double *coef = cpc_coef( ccCuts, j );
                    double rhs = cpc_rhs( ccCuts, j );
                    for ( int l=0 ; (l<nz) ; ++l )
                    {
                        int ivc = idx[l];
                        viol += x[ivc] * coef[l];
                    }
                    viol -= rhs;
                    if (viol >= 0.001)
                    {
                        int new = cpc_add_cut( gcuts, cpc_nz(ccCuts,j), cpc_idx(ccCuts,j), cpc_coef(ccCuts,j), cpc_rhs(ccCuts,j) );
                        if (new)
                        {
                            ++nCuts;
                            static int idcut = 0;
                            printf("  > cut %d viol %g : ", idcut, viol );
                            char cname[256];
                            for ( int k=0 ; (k<cpc_nz(ccCuts,j)) ; ++k )
                                printf("%+g %s ", cpc_coef(ccCuts,j)[k], lp_col_name(mip,cpc_idx(ccCuts,j)[k], cname) );
                            printf("<= %g\n", cpc_rhs(ccCuts,j) );

                            char rname[256];
                            sprintf( rname, "cpc%10d", ++idcut );

                            lp_add_row( mip, cpc_nz(ccCuts,j), cpc_idx(ccCuts,j), cpc_coef(ccCuts,j), rname, 'L', cpc_rhs(ccCuts,j) );

                        }
                    }
                }
                break;
            }
        }
    }

    vec_int_free( &scols );
    vec_IntPair_free( &spairs );

    return nCuts;
}

