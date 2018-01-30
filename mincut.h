#ifndef MINCUT
#define MINCUT

typedef struct _MinCut MinCut;

/** @brief creates a min cut solver
 * @param nArcs number of arcs
 * @param _tail vector with arc sources
 * @param _head vector with arc destinations
 * @param _cap vector with arc capacities
 * @param s source
 * @param t destination
 * @return minimum cut solver
 **/
MinCut *minc_create( int nArcs, const int tail[], const int head[], const int _cap[], int s, int t );


/** @brief solver the min cut problem
 * @param minc minc solver object
 * @return capacity of of min cut (max flow)
 */
int minc_optimize( MinCut *minc );


/** @brief number of arcs in minimum cut
 * @param minc minc solver object
 * @return number of arcs in minimum cut
 **/
int minc_n_cut( MinCut *minc );


/** @brief returns the i-th arc source in the cut
 * @param minc minc solver object
 * @return source of the i-th arc in the cut
 **/
int minc_cut_arc_source( MinCut *minc, int i );


/** @brief returns the i-th arc destination in the cut
 * @param minc minc solver object
 * @return source of the i-th arc in the cut
 **/
int minc_cut_arc_destination( MinCut *minc, int i );


/** @brief frees memory of mincut solver
 **/
void minc_free( MinCut **_minc );

#endif

