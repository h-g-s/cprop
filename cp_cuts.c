#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "cp_cuts.h"
#include "containers.h"
#include "memory.h"

struct _CPCuts
{
    int hashSize;
    
    Vec_int *cutNz;
    Vec_int *cutStart;
    Vec_int *cutIdx;
    Vec_double *cutCoef;
    Vec_char *cutSense;
    Vec_double *cutRHS;

    // cuts per hash bucket
    Vec_int **cutsBucket;
};


/* checks if two doubles are significantly different */
static char dbl_diff( const double v1, const double v2 );

static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157, 821, 257, 263, 389, 397, \
                                        457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 9, 53, 59  };

static unsigned int cutHash( int nz, const int idx[], char sense );

CPCuts *cpc_create( int hashSize )
{
    CPCuts *res;
    ALLOCATE( res, CPCuts );

    res->hashSize = hashSize;
    ALLOCATE_VECTOR_INI( res->cutsBucket, Vec_int*, hashSize  );

    res->cutNz = vec_int_create();
    res->cutStart = vec_int_create();
    res->cutIdx = vec_int_create();
    res->cutCoef = vec_double_create();
    res->cutSense = vec_char_create();
    res->cutRHS = vec_double_create();

    return res;
}

int cpc_n_cuts( const CPCuts *cp )
{
    return vec_int_size( cp->cutStart );
}

char cpc_add_cut( CPCuts *cp, int nz, const int idx[], const double coef[], char sense, double rhs )
{
    int cutBucket = cutHash( nz, idx, sense ) % cp->hashSize;

    if (cp->cutsBucket[cutBucket] == NULL)
        cp->cutsBucket[cutBucket] = vec_int_create();
    else
    {
        /* checking if this is not a duplicate cut */
        const Vec_int *cb = cp->cutsBucket[cutBucket];
        int ic;
        for ( ic=0 ; (ic<vec_int_size(cb)) ;  )
        {
            int idxCut = vec_int_get( cb, ic );

            int onz = vec_int_get( cp->cutNz, idxCut );

            {
                if (onz != nz)
                    goto CHECK_NEXT_CUT_BUCKET;

                char osense = vec_char_get( cp->cutSense, idxCut );
                if (osense != sense)
                    goto CHECK_NEXT_CUT_BUCKET;

                double orhs = vec_double_get( cp->cutRHS, idxCut );
                if (dbl_diff( rhs, orhs ))
                    goto CHECK_NEXT_CUT_BUCKET;

                int start = vec_int_get( cp->cutStart, idxCut );

                const int *oidx = vec_int_getp( cp->cutIdx, start );
                for ( int i=0 ; (i<nz) ; ++i )
                    if (idx[i]!=oidx[i])
                        goto CHECK_NEXT_CUT_BUCKET;

                const double *ocoef = vec_double_getp( cp->cutCoef, start );
                for ( int i=0 ; (i<onz) ; ++i )
                    if (dbl_diff( coef[i], ocoef[i]))
                        goto CHECK_NEXT_CUT_BUCKET;
                }

            /* arrived here, exactly equal to another cut */
            return 0;

CHECK_NEXT_CUT_BUCKET:
            ++ic;
        } // all cuts in this bucket
    } // another cuts in this bucket

    vec_int_push_back( cp->cutsBucket[cutBucket], vec_int_size( cp->cutNz) );

    vec_int_push_back( cp->cutStart, vec_int_size(cp->cutNz)==0 ? 0 : vec_int_last(cp->cutNz)+vec_int_last(cp->cutStart) );
    vec_int_push_back( cp->cutNz, nz );
    vec_char_push_back( cp->cutSense, sense );
    vec_double_push_back( cp->cutRHS, rhs );

    for ( int i=0 ; (i<nz) ; ++i )
        vec_int_push_back( cp->cutIdx, idx[i] );

    for ( int i=0 ; (i<nz) ; ++i )
        vec_double_push_back( cp->cutCoef, coef[i] );

    return 1;
}

static unsigned int cutHash( int nz, const int idx[], char sense )
{
    unsigned int res = 0;
    int nVals = sizeof(hashval) / sizeof(unsigned int);
    assert( nVals >= 2 );
    int i=0;
    for ( ; (i<nz) ; ++i )
        res += idx[i]*hashval[i%nVals];

    res += hashval[(++i)%nVals]*nz;

    res += hashval[(++i)%nVals]*sense;

    return res;
}

static char dbl_diff( const double v1, const double v2 )
{
    return (fabs(v1-v2)>=1e-11);
}

int cpc_nz( const CPCuts *cp, int idxCut )
{
    return vec_int_get( cp->cutNz, idxCut );
}

int *cpc_idx( const CPCuts *cp, int idxCut )
{
    return vec_int_getp( cp->cutIdx, vec_int_get( cp->cutStart, idxCut ) );
}

double *cpc_coef( const CPCuts *cp, int idxCut )
{
    return vec_double_getp( cp->cutCoef, vec_int_get( cp->cutStart, idxCut ) );
}

char cpc_sense( const CPCuts *cp, int idxCut )
{
    return vec_char_get( cp->cutSense, idxCut );
}

double cpc_rhs( const CPCuts *cp, int idxCut )
{
    return vec_double_get( cp->cutRHS, idxCut );
}

void cpc_clear( CPCuts *cpcuts )
{
    vec_int_clear(cpcuts->cutNz);
    vec_int_clear(cpcuts->cutStart);
    vec_int_clear(cpcuts->cutIdx);
    vec_double_clear(cpcuts->cutCoef);
    vec_char_clear(cpcuts->cutSense);
    vec_double_clear(cpcuts->cutRHS);

    for ( int i=0 ; (i<cpcuts->hashSize) ; ++i )
        if (cpcuts->cutsBucket[i])
            vec_int_free( &cpcuts->cutsBucket[i] );

}

void cpc_free( CPCuts **pcp )
{
    CPCuts *cp = *pcp;
    
    for ( int i=0 ; (i<cp->hashSize) ; ++i )
        if (cp->cutsBucket[i])
            vec_int_free( &cp->cutsBucket[i] );

    free( cp->cutsBucket);

    vec_int_free( &cp->cutNz );
    vec_int_free( &cp->cutStart );
    vec_int_free( &cp->cutIdx );
    vec_double_free( &cp->cutCoef );
    vec_double_free( &cp->cutRHS );
    vec_char_free( &cp->cutSense );

    free( cp );

    *pcp = NULL;
}

