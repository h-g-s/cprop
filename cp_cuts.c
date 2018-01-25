#define _GNU_SOURCE
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
    Vec_double *cutRHS;

    // cuts per hash bucket
    Vec_int **cutsBucket;
};


/* checks if two doubles are significantly different */
static char dbl_diff( const double v1, const double v2 );

static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157, 821, 257, 263, 389, 397, \
                                        457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 9, 53, 59  };

static unsigned int cutHash( int nz, const int idx[] );

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
    res->cutRHS = vec_double_create();

    return res;
}

int cpc_n_cuts( const CPCuts *cp )
{
    return vec_int_size( cp->cutStart );
}

struct CutNZs
{
    int idx;
    double coef;
};

static int cmp_cutnz( const void *v1, const void *v2, void *extra )
{
    const struct CutNZs *cnz1 = (const struct CutNZs *) v1;
    const struct CutNZs *cnz2 = (const struct CutNZs *) v2;

    return cnz2->idx - cnz1->idx;
}

char cpc_add_cut( CPCuts *cp, int nz, const int _idx[], const double _coef[], double rhs )
{
    int *sidx = NULL;
    double *scoef = NULL;
    struct CutNZs *scutnz = NULL;

    char addedCut = 0;

    /* checking if cut indexes are sorted */
    char sorted = True;
    for ( int i=0 ; (i<nz-1) ; ++i )
    {
        if ( _idx[i+1] < _idx[i] )
        {
            sorted = False;
            break;
        }
    }

    if (!sorted)
    {
        ALLOCATE_VECTOR( scutnz, struct CutNZs, nz );

        for ( int i=0 ; (i<nz) ; ++i )
        {
            scutnz[i].idx = _idx[i];
            scutnz[i].coef = _coef[i];
        }

        qsort_r( scutnz, nz, sizeof(struct CutNZs), cmp_cutnz, NULL );

        ALLOCATE_VECTOR( sidx, int, nz );
        ALLOCATE_VECTOR( scoef, double, nz );

        for ( int i=0 ; (i<nz) ; ++i )
            sidx[i] = scutnz[i].idx;

        for ( int i=0 ; (i<nz) ; ++i )
            scoef[i] = scutnz[i].coef;
    }

    const int *idx = sorted ? _idx : sidx;
    const double *coef =  sorted ? _coef : scoef;

    int cutBucket = cutHash( nz, idx ) % cp->hashSize;

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
            addedCut = 0;
            goto RETURN_POINT;

CHECK_NEXT_CUT_BUCKET:
            ++ic;
        } // all cuts in this bucket
    } // another cuts in this bucket

    vec_int_push_back( cp->cutsBucket[cutBucket], vec_int_size( cp->cutNz) );

    vec_int_push_back( cp->cutStart, vec_int_size(cp->cutNz)==0 ? 0 : vec_int_last(cp->cutNz)+vec_int_last(cp->cutStart) );
    vec_int_push_back( cp->cutNz, nz );
    vec_double_push_back( cp->cutRHS, rhs );

    for ( int i=0 ; (i<nz) ; ++i )
        vec_int_push_back( cp->cutIdx, idx[i] );

    for ( int i=0 ; (i<nz) ; ++i )
        vec_double_push_back( cp->cutCoef, coef[i] );

    addedCut = 1;

RETURN_POINT:

    if (sidx)
    {
        free( sidx );
        free( scoef );
        free( scutnz );
    }

    return addedCut;
}

static unsigned int cutHash( int nz, const int idx[] )
{
    unsigned int res = 0;
    int nVals = sizeof(hashval) / sizeof(unsigned int);
    assert( nVals >= 2 );
    int i=0;
    for ( ; (i<nz) ; ++i )
        res += idx[i]*hashval[i%nVals];

    res += hashval[(++i)%nVals]*nz;

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

    free( cp );

    *pcp = NULL;
}

