#include "lp.h"
#include "cprop.h"

/** @brief creates a constraint propagation object from mip
 **/
CProp *cprop_create_from_mip( LinearProgram *mip, char verbose );

/** @brief preprocess mip using constraint propagation
 * @param mip mixed integer program with some binary variables
 * @param removeVars after pre-processing variables with implied bounds can be removed or just fixed
 * @param origCols if columns were removed, keeps relationship between new and old variables
 * @return the pre-processed problem or NULL if no implication was discovered or if cprop indicates that the problem is infeasible
 */
LinearProgram *cprop_preprocess( const CProp *cprop, LinearProgram *mip, char removeVars, int origCols[] );
