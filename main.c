/*! \file main.c
    \brief example of use of the constraint propagation engine in a depth-first-search based
    branch and bound
*/


#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "lp.h"
/*
#ifdef __cplusplus
extern "C" {
#include "lp.h"
}
#else
#endif*/

#include "macros.h"
#include "containers.h"
#include "cprop.h"
#include "cp_cuts.h"

int maxDepth = INT_MAX;

double *best = NULL;
double bestObj = DBL_MAX;

// explores a node of the branch and bound tree
void exploreNode( LinearProgram *mip, int depth, CProp *cprop );

void printIdentDepth( int depth );

// returns the index of the most fractional var in mip, or -1 if all variables are integral
int mostFractionalVar( LinearProgram *mip );

int main( int argc, char **argv )
{
    LinearProgram *mip = lp_create();
    
    if (argc<3)
    {
        fprintf( stderr, "usage: \n\tcprop lpFileName maxDepth\n\n");
        exit( EXIT_FAILURE );
    }
    
    lp_read( mip, argv[1] );
    lp_set_print_messages( mip, 0 );

    maxDepth = atoi( argv[2] );

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

    cprop_set_verbose( cprop, 0 );
        
    // adding constraints
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        char rname[256];
        cprop_add_constraint( cprop, nz, idx, coef, lp_sense(mip,i), lp_rhs(mip,i), lp_row_name(mip, i, rname) );
        if (!cprop_feasible(cprop))
            goto END;
    }

    exploreNode( mip, 0, cprop );


END:

     
    cprop_free( &cprop );
    
    strv_free( &names );
    free( lb );
    free( ub );
    free( integer );
    free( idx );
    free( coef );
    lp_free( &mip );
    if (best)
        free(best);
}

int mostFractionalVar( LinearProgram *mip )
{
    const double *x = lp_x( mip );
    double mostFrac = -1.0;
    int jmf = -1;


    for ( int j=0 ; (j<lp_cols(mip)) ; ++j )
    {
        if (!lp_is_integer(mip,j))
            continue;

        // zero
        if ( fabs(x[j]) <= 1e-10 )
            continue;

        double down = floor( x[j] );
        double up = ceil( x[j] );
        const double distDown = x[j] - down;
        const double distUp = up-x[j];

        const double dist = MIN( distDown, distUp );

        if ( fabs(dist)<=1e-6 ) // not fractional
            continue;

        if ( dist > mostFrac )
        {
            jmf = j;
            mostFrac = dist;
        }
    }

    return jmf;
}


void exploreNode( LinearProgram *mip, int depth, CProp *cprop )
{
    if (depth>maxDepth)
        return;

    double objValue = DBL_MAX;

    int status = lp_optimize_as_continuous( mip );
    switch (status)
    {
    case LP_OPTIMAL:
        objValue = lp_obj_value(mip);
        goto PROCESS_NODE;
        break;
    case LP_INFEASIBLE:
        printIdentDepth( depth );
        printf("INFEASIBLE lp.\n");
        return;
    }
    return;



    int jf;
PROCESS_NODE:
    printIdentDepth( depth );
    printf("node obj val: %g", lp_obj_value(mip) );
    if (objValue+1e-8>=bestObj)
    {
        printf(" pruned by bound\n");
        return;
    }
    printf("\n");

    jf = mostFractionalVar( mip );
    if (jf==-1)
    {
        printIdentDepth( depth );
        printf("INTEGER FEASIBLE solution with cost %g found\n", objValue );
        if (lp_obj_value(mip)<bestObj)
        {
            if (!best)
            {
                ALLOCATE_VECTOR( best, double, lp_cols(mip) );
            }
            memcpy( best, lp_x(mip), sizeof(double)*lp_cols(mip));
            bestObj = lp_obj_value(mip);
        }
        return;
    }

    const double *x = lp_x( mip );

    double newBound[] = { ceil(x[jf]), floor(x[jf]) };
    double fvar = x[jf];

    CProp *back = NULL;
    /* branching */
    for ( int b=0 ; b<2 ; ++b )
    {
        printIdentDepth( depth );
        /* best may have improved since last branch */
        if (objValue+1e-8>=bestObj)
        {
            printf("pruned by bound\n");
            return;
        }
        char cname[256];
        const double newB = newBound[b];
        printf("Branching %s%s%g (frac %g)\n", lp_col_name(mip,jf,cname), (!b) ? ">=" : "<=" , newB, fvar );

        if (lp_is_binary(mip, jf))  // validating in cprop first
        {
            int nCutsBefore = cpc_n_cuts( cprop_cut_pool(cprop) );
#ifdef DEBUG
            back = cprop_clone( cprop );
#endif
            cprop_update_bound( cprop, jf,  newB, newB );
            if (!cprop_feasible( cprop ))
            {
                printIdentDepth( depth );
                printf("INFEASIBILITY DETECTED with cprop while branching\n");
                
                int nNewCuts = cpc_n_cuts( cprop_cut_pool(cprop) ) - nCutsBefore;
                if ( nNewCuts > 0 )
                {
                    printIdentDepth( depth );
                    printf("%d new cuts:\n", nNewCuts );
                    const CPCuts *cp = cprop_cut_pool(cprop);
                    for ( int icut=nCutsBefore ; icut<cpc_n_cuts( cp ) ; ++icut )
                    {
                        int nz = cpc_nz( cp, icut );
                        const int *idx = cpc_idx( cp, icut );
                        const double *coef = cpc_coef( cp, icut );
                        char sense = cpc_sense( cp, icut );
                        double rhs = cpc_rhs( cp, icut );
                        printIdentDepth( depth );
                        for ( int j=0 ; j<nz ; ++j )
                            printf("%+g %s ", coef[j], lp_col_name(mip, idx[j], cname) );
                        printf("%s %g\n", sense=='E' ? "=" : sense == 'L' ? "<=" : ">=", rhs );
                    }
                }
                
                cprop_undo( cprop );
#ifdef DEBUG
                if (!cprop_equals(cprop,back))
                    abort();
                cprop_free( &back );
#endif
                continue;
            }
            else
            {
                // checking if there are other implied bounds
                if (cprop_n_implications(cprop))
                {
                    printIdentDepth(depth);
                    printf("CProp Implications: ");
                    int i;
                    for (i=0 ; (i<cprop_n_implications(cprop) && i<5) ; ++i )
                    {
                        int iv = cprop_implied_var( cprop, i );
                        printf("%s=%g ", lp_col_name(mip, iv, cname), cprop_get_lb(cprop,iv) );
                        assert( fabs(cprop_get_lb(cprop,iv)-cprop_get_ub(cprop,iv))<=1e-10 );
                        lp_fix_col( mip, iv, cprop_get_lb(cprop,iv) );
                    }
                    if (cprop_n_implications(cprop)>5)
                        printf("... (more %d)", cprop_n_implications(cprop)-5 );
                    printf("\n");
                }
            }
        }

        double oldBound = 0.0;
        if (!b)
        {
            oldBound = lp_col_lb( mip, jf );
            lp_set_col_bounds( mip, jf, newB, lp_col_ub(mip,jf) );
        }
        else
        {
            oldBound = lp_col_ub( mip, jf );
            lp_set_col_bounds( mip, jf, lp_col_lb(mip,jf), newB );
        }

        exploreNode( mip, depth+1, cprop );

        for (int i=0 ; (i<cprop_n_implications(cprop)) ; ++i )
        {
            int iv = cprop_implied_var( cprop, i );
            lp_set_col_bounds( mip, iv, 0.0, 1.0 );
         }
        cprop_undo( cprop );
#ifdef DEBUG
        if (lp_is_binary(mip, jf))
        {
            assert(back!=NULL);
            if (!cprop_equals(cprop,back))
                abort();
            cprop_free( &back );
        }
#endif

        if (!b)
            lp_set_col_bounds( mip, jf, oldBound, lp_col_ub(mip,jf) );
        else
            lp_set_col_bounds( mip, jf, lp_col_lb(mip,jf), oldBound );
    } // branching up and down
    
} // node exploration


void printIdentDepth( int depth )
{
    for ( int i=0 ; (i<depth) ; ++i )
        printf("   ");
}

