#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "lp.h"
#include "macros.h"
#include "cprop_lp.h"

#define MAXROWS 999999

int main( )
{

    srand(11544123);
    double bestEval = 0.0;

    int rr[MAXROWS];
    for ( int i=0 ; (i<MAXROWS) ; ++i )
        rr[i] = i;

    int n=20;
    int m=13;
    int A[m][n];
    int b[m];

    char cname[][64] = { "x0", "x1", "x2", "x3",
        "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19" };

    char **names =(char **) malloc( sizeof(char*)*n );
    for ( int i=0 ; i<n ; ++i )
        names[i] = &cname[i][0];

    double lb[] = { 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0};
    double ub[] = { 1.0, 1.0, 1.0, 1.0, 
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 , 1.0, 1.0, 1.0, 1.0 };
    double obj[] = { 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 , 1.0, 1.0, 1.0, 1.0 };

    char integer[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };



    for ( int p=0 ; (p<9999999) ; ++p )
    {
        LinearProgram *bip = lp_create();

        lp_add_cols( bip, n, obj, lb, ub, integer, names );

        double eval = 0;
        double obj1 = 0, obj2 = 0;
        int status = LP_INFEASIBLE, st = LP_INFEASIBLE;
        int newCuts = 0;
        int maxRCuts[LP_CUT_TYPES] = { 0, 0, 0, 0, 0, 10, 10, 0 };
        printf("generating bip %d\n", p );

        for ( int i=0 ; (i<m) ; ++i )
        {
GENLHS:
            {
                double sumNeg = 0.0, sumPos = 0.0;
                int nz = 0;
                for ( int j=0 ; (j<n) ; ++j )
                {
                    char appears = INT_RANDOM(3) ? 0 : 1;

                    int sign = INT_RANDOM( 2 ) ? -1 : 1;
                    A[i][j] = appears ? sign * (1+INT_RANDOM(10)) : 0;
                    if (abs(A[i][j])>0)
                    {
                        if ( A[i][j] <= -1e-10 )
                            sumNeg += A[i][j];
                        else
                            if ( A[i][j] >= 1e-10 )
                                sumPos += A[i][j];
                        ++nz;
                    }
                }
                if ((nz<=1) || (nz>6))
                    goto GENLHS;

                double dif = sumPos - sumNeg;

                if ( dif < 3 )
                    goto GENLHS;
                assert( fabs( floor( dif+0.5 ) - dif )<1e-10 ); // should be integral

                b[i] = sumNeg+INT_RANDOM(dif);
                printf(".");
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

        {
            CProp *cprop = cprop_create_from_mip( bip, 0 );
            char feas = cprop_feasible(cprop);
            // counting implications
            int nfix = 0;
            for ( int i=0 ; (i<n) ; ++i )
                if (fabs(cprop_get_lb(cprop,i)-cprop_get_ub(cprop,i))<1e-10)
                    ++nfix;

            cprop_free( &cprop );

            if (  !feas )
                goto NEXT_BIP;

            if (nfix==0 || nfix>10)
                goto NEXT_BIP;
        }

        status = lp_optimize_as_continuous( bip );
        if ( status != LP_OPTIMAL )
            goto NEXT_BIP;

        obj1 = lp_obj_value( bip );

        newCuts = lp_strengthen_with_cuts( bip,  maxRCuts );

        st = lp_optimize( bip );
        if (st != LP_FEASIBLE && st != LP_OPTIMAL)
            goto NEXT_BIP;

        obj2 = lp_obj_value( bip );
        //assert( obj2+1e-10 >= obj1 );

        eval = (obj2-obj1) + newCuts;

        if (eval>bestEval)
        {
            bestEval = eval;
            lp_write_lp( bip, "bip" );
        }


NEXT_BIP:
        lp_free( &bip );

    }

}

