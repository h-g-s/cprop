#include <stdio.h>
#include <stdlib.h>
#include "cprop.h"
#include "lp.h"
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
        
    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
    }
    
    // simulating branching
    cprop_update_bound( cprop, 0, 1.0, 1.0 );
    cprop_update_bound( cprop, 1, 1.0, 1.0 );
    
    cprop_free( &cprop );
    
    strv_free( &names );
    free( lb );
    free( ub );
    free( integer );
    free( idx );
    free( coef );
    lp_free( &mip );
}
