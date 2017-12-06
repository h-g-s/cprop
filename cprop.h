#ifndef CPROP
#define CPROP

typedef struct _CProp CProp;


/** @brief creates a new CProp object
 *
 *  creates a new Constraint Propagation object
 *  entering information about all variables
 *
 *  @param cols number of variables
 *  @param integer indicates if each variable is integer or not
 *  @param lb lower bound of each variable
 *  @param ub upper bound of each variable
 *  @param name (optional) names of variables, useful for debuging
 **/
CProp *cprop_create( int cols, const char integer[], const double lb[], const double ub[], const char **name );


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


/** @brief updates the bounds of a variable
 *
 * Updates the bounds of variable j to [l, u] and performs Constraint Propagation.
 * returns -1 if the problem became unfeasible, 0 if no implications were found
 * or >0 if constraint propagation changed additional bounds
 * all the effects of this bound change can be undone calling cprop_undo. 
 *
 * @param CProp constraint propagation object
 * @param j variable index
 * @l lower bound
 * @u upper bound
 **/
int cprop_update_bound( CProp *cprop, int j, double l, double u );


/** @brief undo the last bound changed and all its implications
 *
 * Undo the last bound changed and all its implications. Can be
 * called multiple times (unlimited undo) 
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


/** @brief returns how many implications the last change produced
 * 
 * returns how many implications the last action (add constraint or change bound) produced
 * 
 * @param cprop constraint propagation object
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

/** @brief frees memory of cprop object 
 *
 * Frees memory of cprop object and sets it to NULL
 *
 **/
void cprop_free( CProp **cprop );

#endif
