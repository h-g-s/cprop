#include <stdio.h>
#include <stdlib.h>
#include "cprop.h"
#ifdef __cplusplus
extern "C" {
#include "lp.h"
}
#else
#include "lp.h"
#endif
#include "macros.h"
#include "containers.h"

int main( int argc, char **argv )
{
    LinearProgram *mip = lp_create();
    lp_read( mip, argv[1] );
    
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

        cprop_set_verbose( cprop, 1 );
        
    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
    }
    
    
    // simulating branching
    printf("Fixing x1=1\n");
    cprop_update_bound( cprop, 0, 1.0, 1.0 );
    
    if (cprop_n_implications(cprop))
    {
        printf("Last Operation Produced %d implications (variables): ", cprop_n_implications(cprop) );
        for ( int i=0 ; (i<cprop_n_implications(cprop)) ; ++i )
            printf( "%d ", cprop_implied_var(cprop, i) );
        printf("\n");
        
    }
    
    printf("Fixing x2=1\n");
    cprop_update_bound( cprop, 1, 1.0, 1.0 );
    printf("\n\n");
    printf("\n\n");
    if (cprop_n_implications(cprop))
    {
        printf("Last Operation Produced %d implications (variables): ", cprop_n_implications(cprop) );
        for ( int i=0 ; (i<cprop_n_implications(cprop)) ; ++i )
        {
            char cname[256];
            printf( "%s ", lp_col_name( mip, cprop_implied_var(cprop,i), cname) );
            
        }
        printf("\n");
        
    }

    cprop_save_impl_graph( cprop, "impl.dot");

     
    cprop_free( &cprop );
    
    strv_free( &names );
    free( lb );
    free( ub );
    free( integer );
    free( idx );
    free( coef );
    lp_free( &mip );
}
