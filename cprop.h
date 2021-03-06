#ifndef CPROP
#define CPROP

/*! \file cprop.h
    \brief constraint propagation engine API
    
  A simple, solver independent, Constraint Propagation engine for binary programs based on the technical report 
  "Nogood Learning for Mixed Integer Programming", from Sandhol, T. and Shields, R. (2006) 
  [https://www.cs.cmu.edu/~sandholm/nogoodsForMip.techReport06.pdf].
 
*/



typedef struct _CProp CProp;


#include "cp_cuts.h"


/** @brief creates a new CProp object
 *
 *  Creates a new Constraint Propagation object
 *  entering information about all variables.   
 *
 *  @param cols number of variables
 *  @param integer indicates if each variable is integer or not
 *  @param lb lower bound of each variable
 *  @param ub upper bound of each variable
 *  @param name names of variables, useful for debuging
 *  @return constraint propagation object
 **/
CProp *cprop_create( int cols, const char integer[], const double lb[], const double ub[], const char **name );


/** @brief creates a copy of a constraint propagation object 
 *
 *  @param CProp constraint propagation object
 *  @return constraint propagation object
 *
 **/
CProp *cprop_clone( const CProp *_cprop );


/** @brief adds a new constraint
 *
 *  Adds a new constraint and performs constraint propagation (CP). Returns 0 
 *  if not implication was found, -1 if CP indicates that the problem 
 *  became unfeasible or >0 if new implications (changed variable bounds)
 *  were generated. Complexity is O(nz) (unless there are new implications in bounds).
 *
 *  @param CProp constraint propagation object
 *  @param nz non-zero variables in this constraint
 *  @param idx indice of each non-zero variable
 *  @param coef coefficient of each non-zero variable
 *  @param sense E equal G greater-or-equal L less-or-equal
 *  @param rhs right-hand-side
 *  @param rname row name
 **/
int cprop_add_constraint( CProp *cprop, int nz, const int idx[], const double coef[], char sense, double rhs, const char rname[] );


/** @brief number of constraints stored in cprop
 *
 *  @param CProp constraint propagation object
 **/
int cprop_n_rows( const CProp *cprop );


/** @brief number of columns (variables) in cprop
 *
 *  @param CProp constraint propagation object
 **/
int cprop_n_cols( const CProp *cprop );


/** @brief number of non-zeros in constraint irow
 *
 *  @param CProp constraint propagation object
 *  @param irow row index
 **/
int cprop_nz( const CProp *cprop, int irow );


/** @brief column indexes in constraint irow
 *
 *  @param CProp constraint propagation object
 *  @param irow row index
 **/
const int *cprop_idx( const CProp *cprop, int row );


/** @brief column coefficients in constraint irow
 *
 *  @param CProp constraint propagation object
 *  @param irow row index
 **/
const double *cprop_coef( const CProp *cprop, int row );


/** @brief number of constraints stored in cprop
 *
 *  @param CProp constraint propagation object
 *  @param irow row index
 **/
double cprop_rhs( const CProp *cprop, int row );


/** @brief should be called after all initial constraints are added
 *
 *  @param CProp constraint propagation object
 **/
void cprop_conclude_pre_processing( CProp *cprop );


/** @brief returns a vector indicating if column was fixed at pre-processing
 * 
 */
const char *cprop_fixed_at_pre_proc( const CProp *cprop );


/** @brief returns how many variables were fixed at the pre-processing
 * 
 */
int cprop_n_fixed_at_pre_proc( const CProp *cprop );



/** @brief updates the bounds of a variable
 *
 * Updates the bounds of variable j to [l, u] and performs Constraint Propagation.
 * returns -1 if the problem became unfeasible, 0 if no implications were found
 * or >0 if constraint propagation changed additional bounds
 * all the effects of this bound change can be undone calling cprop_undo. 
 *
 * @param CProp constraint propagation object
 * @param j variable index
 * @param l lower bound
 * @param u upper bound
 **/
int cprop_update_bound( CProp *cprop, int j, double l, double u );


/** @brief undo the last bound changed and all its implications
 *
 * Undo the last bound changed and all its implications. Can be
 * called multiple times (unlimited undo) to go back to the initial state.
 *
 * @param cprop constraint propagation object
 *
 * */
void cprop_undo( CProp *cprop );


/** @brief returns 1 if problem is still feasible, zero otherwise 
 * 
 * @param cprop constraint propagation object
 **/
char cprop_feasible( const CProp *cprop );


/** @brief message explaining how infeasibility was detected
 * 
 * @param cprop constraint propagation object
 * 
 **/
const char *cprop_inf_msg( const CProp *cprop );


/** @brief returns how many implications the last change produced
 * 
 * 
 * @param cprop constraint propagation object
 * @return returns how many implications the last action (add constraint or change bound) produced
 * 
 * */
int cprop_n_implications( const CProp *cprop );

/** @brief returns the i-th variable which bound was implied in the last operation
 * 
 * returns the i-th variable which bound was implied in the last operation
 * 
 * @param cprop constraint propagation object
 * @param i the i-th variable which a bound was changed was implied in the last operation
 * 
 **/
int cprop_implied_var( const CProp *cprop, int i );


/** @brief returns the current lower bound for column j
 * 
 * @param cprop constraint propagation object
 * @param j index of variable
 * */
double cprop_get_lb( const CProp *cprop, int j );


/** @brief returns the current upper bound for column j
 * 
 * @param cprop constraint propagation object
 * @param j index of variable
 * 
 * */
double cprop_get_ub( const CProp *cprop, int j );


/** @brief checks if variable j is fixed or not
 * 
 * @param cprop constraint propagation object
 * @param j index of variable
 */
char cprop_fixed( const CProp *cprop,  int j );

/// Implication Graph Node Types
enum IGNType { 
    EOne,      /*!< variable fixed at one */
    EZero,     /*!< variable fixed at zero */
    Infeasible /*!< infeasible node */
};


/** @brief returns the node type 
 * 
 * Returns the node type of a node in the implication graph
 *
 * @param cprop constraint propagation object
 * @param nodeId id of the node in the graph
 **/
enum IGNType cprop_impl_graph_node_type( const CProp *cprop, int nodeId );


/** @brief returns the variable which a node refers
 * 
 * Returns the variable thar a node in the implication graph refers to
 *
 * @param cprop constraint propagation object
 * @param nodeId id of the node in the graph
 **/
int cprop_impl_graph_node_var( const CProp *cprop, int nodeId );


/** @brief returns the node id 
 * 
 * Returns the node id for a node of type ntype related to column col
 *
 * @param cprop constraint propagation object
 * @param ntype node type
 * @param col variable (column) index
 **/
int cprop_impl_graph_node_id( const CProp *cprop, enum IGNType ntype , int col );


/** @brief number of incident arcs in implication graph node
 *
 * Returns the number of incident arcs in the implication graph for a node with id nodeId
 * col is ignored if a node of type Infeasible is informed 
 *
 * @param cprop Constraint Propagation Object
 * @param nodeId nodeId
 **/
int cprop_impl_graph_in_d( const CProp *cprop, int nodeId );


/** @brief returns the origin of the i-th incident arc to nodeId
 *
 * Returns the nodeId of the origin of the i-th incident arc in node nodeId
 *
 * @param cprop Constraint Propagation Object
 * @param nodeid nodeId
 **/
int cprop_impl_graph_in_neigh( const CProp *cprop, int nodeId, int i );


/** @brief number of outbount arcs in implication graph node
 *
 * Returns the number of outbond arcs in the implication graph for a node with id nodeId
 * col is ignored if a node of type Infeasible is informed 
 *
 * @param cprop Constraint Propagation Object
 * @param nodeId nodeId
 **/
int cprop_impl_graph_out_d( const CProp *cprop, int nodeId );


/** @brief returns the destination of the i-th outbound arc to nodeId
 *
 * Returns the nodeId of the destination of the i-th outbound arc in node nodeId
 *
 * @param cprop Constraint Propagation Object
 * @param nodeid nodeId
 **/
int cprop_impl_graph_out_neigh( const CProp *cprop, int nodeId, int i );


/** @brief saves the implication graph in the DOT file format (graphviz)
 *
 * @param cprop constraint propagation object
 * @param fName file name (.dot extension recommended)
 **/
void cprop_save_impl_graph( const CProp *cprop, const char *fName );



/** @brief activates detailed printing of implications as they are discovered
 *
 * Activates detailed printing of implications and bound changes as they are discovered.
 * Useful to debug.
 * 
 * @param cprop Constraint Propagation Object
 * @param verbose if messages will be printed on the screen
 *  
 **/
void cprop_set_verbose( CProp *cprop, char verbose );


/// Diving strategy
enum DVCStrat { 
    DVCMostFrac,      /*!< starts from the most fractional variables */
    DVCRandom         /*!< randomized diving */
};

/// Diving direction
enum DVDir {
    DVDClosestInteger,    /*!< closest integer */
    DVDBoth               /*!< both */
};


/* @brief enters fractional solution information. whith this information only violated cuts are stored
 *  
 * @param cprop Constraint Propagation Object
 **/
void cprop_enter_relax_sol( CProp *cprop, const double x[] );


/* @brief checks the effect of fixing a variable in a value, generating violated cuts if possible
 * @param cprop Constraint Propagation Object
 * @param var variable index
 * @param value value which variable will be temporaly fixed
 * @return number of new violated cuts found
 **/
int cprop_probe( CProp *cprop, int var, int value );


/* @brief probes the influence of fixing variables in zero or one, trying to generate cuts
 *  
 * @param cprop Constraint Propagation Object
 **/
int cprop_probe_and_cut( CProp *cprop, const double x[], int maxTries, enum DVCStrat dvstrat, char onlyFrac, enum DVDir direction );


/** @brief returns the cut pool with all cuts discovered so far
 *
 * @param cprop Constraint Propagation Object
 */
const CPCuts *cprop_cut_pool( CProp *cprop );


/** @brief compares two constraint propagation objects
 *
 **/
char cprop_equals( const CProp *cprop1, const CProp *cprop2 );


/** @brief goes back to the original state
 *
 **/
void cprop_clear( CProp *cprop );



/** @brief returns the name of row 0
 * @param cprop Constraint Propagation Object
 * @param irow row index
 * @return row name
 */
const char *cprop_row_name( const CProp *cprop, int i );



/** @brief checks if constraint irow is part of an equality (with row irow+1)
 * @param cprop Constraint Propagation Object
 * 
 **/
char cprop_is_equality( const CProp *cprop, int irow );


/** @brief returns the number of non-zeros in constraint irow
 * @param cprop Constraint Propagation Object
 * @param irow row index
 * @return number of variables in row irow
 */
int cprop_nz( const CProp *cprop, int irow );


/** @brief returns the column indexes in row irow
 * @param cprop Constraint Propagation Object
 * @param irow row index
 * @return column indexes vector
 */
const int *cprop_idx( const CProp *cprop, int row );


/** @brief returns the column cofficients in row irow
 * @param cprop Constraint Propagation Object
 * @param irow row index
 * @return column coefficients vector
 */
const double *cprop_coef( const CProp *cprop, int row );


/** @brief returns the right-and-side of row irow
 * @param cprop Constraint Propagation Object
 * @param irow row index
 * @return right hand side
 */

double cprop_rhs( const CProp *cprop, int row );



int cprop_n_rows_col( const CProp *cprop, int col );
int *cprop_rows_col( const CProp *cprop, int col );


/** @brief frees memory of cprop object 
 *
 * Frees memory of cprop object and sets it to NULL
 *
 **/
void cprop_free( CProp **cprop );



void cprop_print_impl( const CProp *cprop );


#endif

