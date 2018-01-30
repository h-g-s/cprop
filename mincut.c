#include <string.h>
#include <assert.h>
#include <limits.h>
#include "mincut.h"
#include "macros.h"
#include "containers.h"

struct MinCArc
{
    int v;
    int cap;
    int rpos; // position of reverse arc
    char original; // if arc exists in original graph
};

struct _MinCut
{
    int n;

    // start[u] indicates the starting position outbond arcs from u
    int *start;

    struct MinCArc *arcs;

    // maps to original node indexes,
    // just to output the solution
    int *orig;

    // incidence vector and list of visited nodes
    char *ivVisited;
    int nVisited;
    int *visited;
    
    // queue of unvisited nodes
    int *queue;

    // path from t to s
    int *parent;

    // to quickly search for an (u,v) arc
    Dict_int *arcDict;
    
    // source node
    int s;

    // target node
    int t;

    // minimum cut answer
    int nCut;
    int *cutU;
    int *cutV;
};

static char *arcName( char *str, int u, int v )
{
    sprintf( str, "(%d,%d)",u,v );
    return str;
}

static int arcPos( const MinCut *minc, int u, int v )
{
    char str[256];
    const Dict_int *arcDict = minc->arcDict;
    const char *aname = arcName( str, u, v );
    return dict_int_get( arcDict, aname );
}

MinCut *minc_create( int nArcs, const int _tail[], const int _head[], const int _cap[], int s, int t )
{
    assert( s!=t );

    int maxN = -1;
    for ( int i=0 ; (i<nArcs) ; ++i )
        maxN = MAX( maxN, _tail[i] );
    for ( int i=0 ; (i<nArcs) ; ++i )
        maxN = MAX( maxN, _head[i] );

    // maps original nodes to new pre-processed nodes
    int *ppnode;
    ALLOCATE_VECTOR( ppnode, int, (maxN+1)  );
    for ( int i=0 ; (i<maxN+1) ; ++i )
        ppnode[i] = -1;

    int n=0;

    for ( int i=0 ; (i<nArcs) ; ++i )
        if (ppnode[_tail[i]]==-1)
            ppnode[_tail[i]] = n++;
    for ( int i=0 ; (i<nArcs) ; ++i )
        if (ppnode[_head[i]]==-1)
            ppnode[_head[i]] = n++;

    // may be larger due to insertion of reverse arcs
    int *tail;
    int *head;
    int *cap;
    ALLOCATE_VECTOR( tail, int, (6*nArcs) );
    head = tail + 2*nArcs;
    cap = head + 2*nArcs;

    // inverse arcs must be added if needed
    Dict_int *arcDict = dict_int_create( nArcs*2, -1 );

    // storing positions of existing arcs
    char str[256];
    for ( int i=0 ; i<nArcs ; ++i )
    {
        tail[i] = ppnode[_tail[i]];
        head[i] = ppnode[_head[i]];
        if (tail[i] == head[i] )
        {
            fprintf( stderr, "minc ERROR: arc (%d,%d) specified (self arc).\n", _tail[i], _head[i] );
            abort();
        }
        cap[i] = _cap[i];
        const char *aname = arcName(str, tail[i], head[i] );
        if ( dict_int_get( arcDict, aname ) != -1 )
        {
            fprintf( stderr, "minc ERROR: arc (%d,%d) specified twice.\n", _tail[i], _head[i] );
            abort();
        }
        dict_int_set( arcDict, aname, i );
    }

    int *orig;
    ALLOCATE_VECTOR( orig, int, n );

    for ( int i=0 ; (i<maxN+1) ; ++i )
        if (ppnode[i]!=-1)
            orig[ppnode[i]] = i;

    MinCut *minc;
    ALLOCATE( minc, MinCut );

    minc->s = ppnode[s];
    minc->t = ppnode[t];

    free( ppnode );

    // adding missing arcs
    int nOrigArcs = nArcs;
    for ( int i=0 ; i<nOrigArcs ; ++i )
    {
        int u = tail[i];
        int v = head[i];
        const char *rname = arcName(str, v, u );

        if ( dict_int_get( arcDict, rname ) == -1 )
        {
            dict_int_set( arcDict, rname, nArcs );
            tail[nArcs] = v;
            head[nArcs] = u;
            cap[nArcs] = 0;
            ++nArcs;
        }
    }

    int *start, *nNeigh;
    struct MinCArc *arcs;
    ALLOCATE_VECTOR( arcs, struct MinCArc, nArcs );

    ALLOCATE_VECTOR( start, int, (n+1) );
    ALLOCATE_VECTOR_INI( nNeigh, int, n );

    // counting neighbors per node
    for ( int i=0 ; (i<nArcs) ; ++i )
        nNeigh[tail[i]]++;

    // setting up start
    start[0] = 0;
    for ( int i=1 ; (i<n+1) ; ++i )
        start[i] = start[i-1] + nNeigh[i-1];

    memset( nNeigh,  0, sizeof(int)*n );

    // storing arcs in positions
    // organized by tail
    for ( int i=0 ; (i<nArcs) ; ++i )
    {
        const int ctail = tail[i];
        const int chead = head[i];
        const int pos = start[ctail]+nNeigh[ctail];
        arcs[pos].v = chead;
        arcs[pos].cap = cap[i];
        const char *aname = arcName( str, ctail, chead );
        dict_int_set( arcDict, aname, pos );

        ++(nNeigh[ctail]);
        
    }

    free( nNeigh );

    // filling reverse arcs positions
    for ( int u=0 ; (u<n) ; ++u )
    {
        for ( int j=start[u] ; j<start[u+1] ; ++j )
        {
            int v = arcs[j].v;
            const char *rname = arcName( str, v, u );
            int rpos = dict_int_get( arcDict, rname );
            assert( rpos >= 0 && rpos < nArcs && rpos != j && rpos >= start[v] );
            assert( arcs[rpos].v == u );
            arcs[j].rpos = rpos;
        } // all v
    } // all us 


    minc->n = n;
    minc->orig = orig;
    minc->arcs = arcs;
    minc->start = start;


    ALLOCATE_VECTOR_INI( minc->ivVisited, char, minc->n );
    ALLOCATE_VECTOR( minc->visited, int, minc->n );
    minc->nVisited = 0;

    ALLOCATE_VECTOR( minc->queue, int, n );
    ALLOCATE_VECTOR( minc->parent, int, n );

    minc->arcDict = arcDict;

    free( tail );

    ALLOCATE_VECTOR( minc->cutU, int, 2*n );
    minc->cutV = minc->cutU + n;

    minc->nCut = 0;
    
    /* info about original arcs */
    for ( int i=0 ; i<nArcs ; ++i )
        arcs[i].original = arcs[i].cap>0;

    return minc;
}

static void addVisited( MinCut *minc, int node )
{
#ifdef DEBUG
    assert( minc->ivVisited[node]==False );
#endif
    minc->ivVisited[node] = True;
    minc->visited[minc->nVisited++] = node;
}

static void clearVisited( MinCut *minc )
{
    for ( int i=0 ; (i<minc->nVisited) ; ++i )
        minc->ivVisited[minc->visited[i]] = False;
    minc->nVisited = 0;
#ifdef DEBUG
    for ( int i=0 ; (i<minc->n) ; ++i )
    {
        assert( minc->ivVisited[i]==False );
    }
#endif
}

static char bfs( MinCut *minc )
{
    minc->nVisited = 0;

    const int s = minc->s;
    const int t = minc->t;

    int *queue = minc->queue;
    int *parent = minc->parent;
    const int *start = minc->start;
    char *ivVisited = minc->ivVisited;
    const struct MinCArc *arcs = minc->arcs;


    queue[0] = s;
    int nQueue = 1;
    addVisited( minc, s );
    parent[s] = -1;

    while ( nQueue>0 )
    {
        int u = queue[--nQueue];

        // exploring neighbors of u
        for ( int p=start[u] ; p<start[u+1] ; ++p )
        {
            int v = arcs[p].v;

            if (ivVisited[v]==False && arcs[p].cap>0 )
            {
                queue[nQueue++] = v;
                parent[v] = u;
                addVisited( minc, v );
            }
        }
    }

    char reachedT = ivVisited[t];

    // clear visited
    clearVisited( minc );

    return reachedT;
}

static void dfs( MinCut *minc, int s )
{
    const int *start = minc->start;
    char *ivVisited = minc->ivVisited;

    addVisited( minc, s );
    
    // checking neighbors
    for ( int j=start[s] ; (j<start[s+1]) ; ++j )
    {
        const struct MinCArc *arc = minc->arcs+j;
        if ( ivVisited[arc->v]==False && arc->cap )
            dfs( minc, arc->v );
    }
}

int minc_optimize( MinCut *minc )
{
    const int s = minc->s;
    const int t = minc->t;
    const int *parent = minc->parent;
    struct MinCArc *arcs = minc->arcs;
    const int *start = minc->start;

    int totalFlow = 0;
    while ( bfs( minc ) )
    {
        int flow = INT_MAX;
       
        for ( int v=t; (v!=s) ; v=parent[v] )
        {
            int u = parent[v];

            int apos = arcPos( minc, u, v );
            flow = MIN( flow, arcs[apos].cap );
        } // checking path capacity
        assert( flow > 0 );
        
        totalFlow += flow;

        // updating residual capacities
        for ( int v=t; (v!=s) ; v=parent[v] )
        {
            int u = parent[v];
            int apos = arcPos( minc, u, v );
            arcs[apos].cap -= flow;

            int rpos = arcPos( minc, v, u );
            arcs[rpos].cap += flow;
        }

    } // while found a path

    if (totalFlow)
    {
        dfs( minc, minc->s );

        const char *ivVisited = minc->ivVisited;

        // checking arc cuts
        for ( int u=0 ; (u<minc->n) ; ++u )
        {
            if (!ivVisited[u])
                continue;

            for ( int pos=start[u] ; (pos<start[u+1]) ; ++pos )
            {
                const struct MinCArc *arc = arcs+pos;
                if (ivVisited[arc->v] || !arc->original)
                    continue;

                minc->cutU[minc->nCut] = minc->orig[u];
                minc->cutV[minc->nCut] =  minc->orig[arc->v];

                ++minc->nCut;
            } // destination side
        } // source side 
    }

    return totalFlow;
}

int minc_n_cut( MinCut *minc )
{
    return minc->nCut;
}


int minc_cut_arc_source( MinCut *minc, int i )
{
    return minc->cutU[i];
}


int minc_cut_arc_destination( MinCut *minc, int i )
{
    return minc->cutV[i];
}


void minc_free( MinCut **_minc )
{
    MinCut *minc = *_minc;

    free( minc->start );
    free( minc->arcs );
    free( minc->orig );
    free( minc->ivVisited );
    free( minc->visited );
    free( minc->queue );
    free( minc->parent );
    dict_int_free( &minc->arcDict );
    free( minc->cutU );
    free( minc );

    *_minc = NULL;
}

