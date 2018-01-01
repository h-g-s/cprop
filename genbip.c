#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "lp.h"
#include "macros.h"

#define MAXROWS 999999

int main( )
{
    double bestEval = 0.0;

    int rr[MAXROWS];
    for ( int i=0 ; (i<MAXROWS) ; ++i )
        rr[i] = i;

    int n=8;
    int m=9;
    int A[m][n];
    int b[m];

    char cname[][64] = { "x0", "x1", "x2", "x3",
                         "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19" };

    char **names = malloc( sizeof(char*)*n );
    for ( int i=0 ; i<n ; ++i )
        names[i] = &cname[i][0];

    double lb[] = { 0.0, 0.0, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0};
    double ub[] = { 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 , 1.0, 1.0, 1.0, 1.0 };
    double obj[] = { 0.0, 0.0, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 , 1.0, 1.0, 1.0, 1.0 };

    char integer[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    LinearProgram *bip = lp_create();

    lp_add_cols( bip, n, obj, lb, ub, integer, names );


    for ( int p=0 ; (p<99900000) ; ++p )
    {
        printf("generating bip %d\n", p );

        for ( int i=0 ; (i<m) ; ++i )
        {
            for ( int j=0 ; (j<n) ; ++j )
            {
                char appears = INT_RANDOM(5) ? 0 : 1;
                int sign = INT_RANDOM( 2 ) ? -1 : 1;
                A[i][j] = appears ? sign * (1+INT_RANDOM(10)) : 0;
                int signRHS  = INT_RANDOM( 2 ) ? -1 : 1;
                b[i] = signRHS * INT_RANDOM(10);
            }
        }

        int nr = 0;
        for ( int i=0 ; (i<m) ; ++i )
        {
            int idx[n]; double coef[n];
            int nz = 0;
            for ( int j=0 ; (j<n) ; ++j )
            {
                if (A[i][j])
                {
                    idx[nz] = j;
                    coef[nz] = A[i][j];
                    ++nz;
                }
            }
            if (!nz) continue;
            char rname[256];
            sprintf( rname, "r%d", nr++ );
            lp_add_row( bip, nz, idx, coef, rname, 'L', b[i] );
        }

        if (nr==0)
            goto NEXT_BIP;

        for ( int i=0 ; i<n ; ++i )
            obj[i] = INT_RANDOM(20)-10;

        lp_set_obj( bip, obj );

            //lp_write_lp( bip, "bug" );
            //lp_write_lp( bip, "bipa" ); exit(1);
        int status = lp_optimize_as_continuous( bip );
        if ( status != LP_OPTIMAL )
            goto NEXT_BIP;

        int obj1 = lp_obj_value( bip );

        int maxRCuts[LP_CUT_TYPES] = { 0, 0, 10, 10, 0, 10, 10, 10 };
        int newCuts = lp_strengthen_with_cuts( bip,  maxRCuts );

        int st = lp_optimize( bip );
        if (st != LP_FEASIBLE && st != LP_OPTIMAL)
            goto NEXT_BIP;

        double obj2 = lp_obj_value( bip );
        //assert( obj2+1e-10 >= obj1 );

        double eval = (obj2-obj1)*100 + newCuts;

        if (eval>bestEval)
        {
            bestEval = eval;
            lp_write_lp( bip, "bip" );
        }


NEXT_BIP:
        lp_remove_rows( bip, lp_rows(bip), rr );

    }
    
    lp_free( &bip );
}

