#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "cprop.h"
#include "memory.h"
#include "macros.h"
#include "containers.h" 


// large values, summation does not gives an overflow
const double oo  = 1e23;
#define VEPS  1e-10
#define FIXED( lb, ub ) ( ub-lb <= VEPS )

enum ConstraintType
{ 
    ConstraintBV,    // constraint involving only binary variables
    ConstraintINNB,  // constraint involving only integer, >= 0 bounded variables (< oo)
    ConstraintOther
} ;



struct _CProp
{
    int cols;

    double *lb;
    double *ub;

    /* original bounds */
    double *olb;
    double *oub;

    char *integer;
    char *binary;

    char **cname;

    /* temporary area */
    int *idx;
    double *coef;

    double *pv;
    double *nv;
    /* incidence vector indicating the processing of rows */
    Vec_char *ivr;
    /* stack of rows to be processed after an update */
    Vec_int *stkRows;

    /* all constraints stored in the format
     * ax <= b */
    Vec_int *vrstart;
    Vec_int *vrnz;
    Vec_int *vridx;
    Vec_double *vrcoef;
    Vec_double *vrrhs;
    Vec_char *vrtype;
    // in which rows a column appears
    V2D_int *rowsCol;
    Vec_int *vrnamest; // start of each row name
    Vec_char *vrnames; // row names

    /* undo information */
    Vec_int *vnimpl;
    Vec_int *vcolimpl;
    Vec_double *voldlb;
    Vec_double *voldub;
    Vec_IntPair *lastArcsIG; // last nodes added to the implication graph
    Vec_int *nNewArcsIG; // number of nodes added in the last operation to the implication graph



    char feasible;

    // Input arcs of implication graph
    // positions j=0 .. cols-1  indicate nodes for
    // x_j=1, positions j=cols .. 2*cols-1 indicate nodes 
    // for x_j = 0, finally, node 2*cols indicate infeasibility
    ISet **implGIn;

    int nimpl;

    Vec_char *msgInf;

    char verbose;
};

int cprop_process_constraint_binary_variables( CProp *cprop, int irow );

int cprop_process_constraint( CProp *cprop, int irow );

/* reports an infeasible constraint considering the minimum value minLhs on the left-hand-side 
 * also adds nodes to the implicating graph if there are fixed variables */
void cprop_report_infeasible_constraint( CProp *cprop, int irow, double minLhs );

#define HASH_SIZE_IMPLG 128

/* maximum number of variables printed 
 * when some information of the row is displayed */
#define MAX_VAR_ROW_PRINT 6

void cprop_add_row_name( CProp *cprop, const char rname[] );
const char *cprop_row_name( const CProp *cprop, int i );
void cprop_add_arc_impl_g( CProp *cprop, enum IGNType ntSource, int colSource, enum IGNType ntDest, int colDest );
void cprop_add_msg_inf( CProp *cprop, const char *msg );
/* trimmed description of constraint to show in messages */
char *cprop_constraint_descr( const CProp *cprop, int irow, char *descr );
/* trimmed description of bounds in variables */
char *cprop_constraint_descr_variable_bounds( const CProp *cprop, int irow, char *descr );


/* stored constraints in the format
 * ax <= b */
void cprop_add_row( CProp *cprop, int nz, const int idx[], const double coef[], double rhs, double mult, const char rname[] );
int cprop_nz( const CProp *cprop, int irow );
const int *cprop_idx( const CProp *cprop, int row );
const double *cprop_coef( const CProp *cprop, int row );
double cprop_rhs( const CProp *cprop, int row );
int cprop_n_rows( const CProp *cprop );

int cprop_n_rows_col( const CProp *cprop, int col );
int *cprop_rows_col( const CProp *cprop, int col );

CProp *cprop_create( int cols, const char integer[], const double lb[], const double ub[], const char **name )
{
    CProp *cprop;
    ALLOCATE( cprop, CProp );

    cprop->cname = NULL;

    cprop->feasible = 1;
    cprop->nimpl = 0;

    cprop->cols = cols;
    ALLOCATE_VECTOR( cprop->lb, double, cols );
    ALLOCATE_VECTOR( cprop->ub, double, cols );
    ALLOCATE_VECTOR( cprop->olb, double, cols );
    ALLOCATE_VECTOR( cprop->oub, double, cols );

    ALLOCATE_VECTOR( cprop->integer, char, cols );
    ALLOCATE_VECTOR( cprop->binary, char, cols );
    memcpy( cprop->lb, lb, sizeof(double)*cols );
    memcpy( cprop->ub, ub, sizeof(double)*cols );
    memcpy( cprop->integer, integer, sizeof(char)*cols );

    int j;
    // preventing under/overflow
    for ( j=0 ; (j<cols) ; ++j )
        cprop->lb[j] = MAX( cprop->lb[j], -oo );
    for ( j=0 ; (j<cols) ; ++j )
        cprop->ub[j] = MIN( cprop->ub[j], oo );
    for ( j=0 ; (j<cols) ; ++j )
        cprop->olb[j] = MAX( cprop->lb[j], -oo );
    for ( j=0 ; (j<cols) ; ++j )
        cprop->oub[j] = MIN( cprop->ub[j], oo );

    for ( j=0 ; (j<cols) ; ++j )
        cprop->binary[j] = integer[j] && lb[j] >= -VEPS && ub[j] <= 1.0+VEPS;

    if (name)
    {
        int strspace = 0;
        int *charsname;
        ALLOCATE_VECTOR_INI( charsname, int, cols );

        for ( j=0 ; (j<cols) ; ++j )
        {
            charsname[j] = strlen( name[j] );
            strspace += charsname[j];
        }

        ALLOCATE_VECTOR( cprop->cname, char *, cols );
        ALLOCATE_VECTOR( cprop->cname[0], char, (strspace+cols) );

        for ( j=1 ; (j<cols) ; ++j )
            cprop->cname[j] = cprop->cname[j-1] + charsname[j-1] + 1;

        for ( j=0 ; (j<cols) ; ++j )
            strcpy( cprop->cname[j], name[j] );

        free( charsname );
    }

    /* temporary area in cprop */
    ALLOCATE_VECTOR( cprop->idx, int, cols );
    ALLOCATE_VECTOR( cprop->coef, double, cols );
    ALLOCATE_VECTOR( cprop->pv, double, cols );
    ALLOCATE_VECTOR( cprop->nv, double, cols );
    FILL( cprop->pv, 0, cols, 0.0 );
    FILL( cprop->nv, 0, cols, 0.0 );
    cprop->ivr = vec_char_create();
    cprop->stkRows = vec_int_create();
    cprop->vrnames = vec_char_create();
    cprop->vrnamest = vec_int_create();

    /* all constraints stored in the format
     * ax <= b */
    /* to store constraints */
    cprop->vrstart = vec_int_create();
    cprop->vrnz = vec_int_create();
    cprop->vridx = vec_int_create();
    cprop->vrcoef = vec_double_create();
    cprop->vrrhs = vec_double_create();
    cprop->vrtype = vec_char_create();
    // in which rows a column appears
    cprop->rowsCol = v2d_create( cols );

    /* undo information */
    cprop->vnimpl = vec_int_create();
    cprop->vcolimpl = vec_int_create();
    cprop->voldlb = vec_double_create();
    cprop->voldub = vec_double_create();

    ALLOCATE_VECTOR( cprop->implGIn, ISet *, 2*cols+1 );
    for ( int j=0 ; (j<2*cols+1) ; ++j )
        cprop->implGIn[j] = NULL;



    cprop->lastArcsIG = vec_IntPair_create();
    cprop->nNewArcsIG = vec_int_create();

    cprop->verbose = 0;

    cprop->msgInf = NULL;

    return cprop;
}


int cprop_process_constraint_binary_variables( CProp *cprop, int irow )
{
    int nz = cprop_nz( cprop, irow );
    const int *idx = cprop_idx( cprop, irow );
    const double *coef = cprop_coef( cprop, irow );
    double rhs = cprop_rhs( cprop, irow );

    int nimpl = 0;

    // minimum positive value and maximum negative value in the LHS
    double minP = 0.0, maxN = 0.0;

    double *pv = cprop->pv; 
    double *nv = cprop->nv;
    double *ub = cprop->ub;
    double *lb = cprop->lb;
    const char **cname = (const char**) cprop->cname;
    Vec_int *vcolimpl = cprop->vcolimpl;
    Vec_double *voldlb = cprop->voldlb;
    Vec_double *voldub = cprop->voldub;

    int nFixed = 0;

    for ( int j=0 ; j<nz ; ++j )
    {        
        pv[j] = nv[j] = 0.0;

        if (coef[j]>=VEPS)
        {
            pv[j] = coef[j]*lb[idx[j]];
            minP += pv[j];
        }
        else
        {
            if (coef[j]<=-VEPS)
            {
                nv[j] = (-1.0*coef[j])*ub[idx[j]];
                maxN += nv[j];
            }
            else
            {
                fprintf( stderr, "constraint %s with coefficient too close to zero: %g\n", cprop_row_name(cprop, irow), coef[j] );
                abort();
            }
        }

        minP = MIN( oo, minP );
        maxN = MIN( oo, maxN );
    } // all non-zeros

    for ( int j=0 ; j<nz ; ++j )
        pv[j] = MIN( pv[j], oo );

    for ( int j=0 ; j<nz ; ++j )
        nv[j] = MIN( nv[j], oo );

    for ( int j=0 ; j<nz ; ++j )
        if (FIXED(lb[idx[j]], ub[idx[j]]))
            ++nFixed;

    /* checking if constraint if unfeasible */
    if ( minP-maxN >= rhs+VEPS )
    {
        cprop_report_infeasible_constraint( cprop, irow, minP-maxN );

        return -1;
    }

    // computing upper bound of each var in this constraint
    for ( int j=0 ; j<nz ; ++j )
    {
        if (FIXED( lb[idx[j]], ub[idx[j]] ))
            continue;

        const double uij = rhs - minP + pv[j] + maxN - nv[j];

        // N+ set
        if ( coef[j] > VEPS )
        {
            if ( uij <= -VEPS )
            {
                char msg[1024], strConstr[512], strFixed[512] = "";
                cprop_constraint_descr( cprop, irow, strConstr );
                sprintf( msg, "No value for variable %s can satisfy constraint %s.\n", cname[idx[j]], strConstr );
                cprop_add_msg_inf( cprop, msg );
                if (nFixed)
                {
                    char strComment[64] = "";
                    cprop_constraint_descr_variable_bounds( cprop, irow, strFixed );
                    if (strstr(strFixed,"*"))
                        strcpy( strComment, " (* implied bound)");
                    sprintf( msg, "Fixed variables: %s %s\n", strFixed, strComment );
                    cprop_add_msg_inf( cprop, msg );
                    if (cprop->verbose)
                        printf("%s", msg );
                }

                // all other fixed variables in this constraint
                for ( int k=0 ; (k<nz) ; ++k )
                {
                    if (k!=j)
                    {
                        cprop_add_arc_impl_g( cprop, 
                                lb[idx[k]]>=0.98 ? EOne : EZero, idx[k], 
                                Infeasible, 0 );
                    }
                }

                return -1;
            }
            else
            {
                if ( coef[j] >= uij+VEPS )
                {
                    vec_int_push_back( vcolimpl, idx[j] );
                    vec_double_push_back( voldlb, lb[idx[j]] );
                    vec_double_push_back( voldub, ub[idx[j]] );
                    ub[idx[j]] = 0.0;
                    if (cprop->verbose)
                        printf(" impl -> %s=0\n", cname[idx[j]] );

                    // all other fixed variables in this constraint
                    for ( int k=0 ; (k<nz) ; ++k )
                    {
                        if (k!=j)
                        {
                            cprop_add_arc_impl_g( cprop, 
                                    lb[idx[k]]>=0.98 ? EOne : EZero, idx[k], 
                                    EZero, idx[k] );
                        }
                    }

                    ++nimpl;
                }
            }
        }
        else
        {
            // set N-
            if ( coef[j] <= -VEPS )
            {
                if ( coef[j] >= uij+VEPS )
                {
                    cprop_report_infeasible_constraint( cprop, irow, 0.0 );

                    return -1;
                }
                else
                {
                    if ( uij <= -VEPS )
                    {
                        vec_int_push_back( vcolimpl, idx[j] );
                        vec_double_push_back( voldlb, lb[idx[j]] );
                        vec_double_push_back( voldub, ub[idx[j]] );
                        ++nimpl;
                        lb[idx[j]] = 1.0;
                        if (cprop->verbose)
                            printf("  impl %s=1\n", cname ? cname[idx[j]] : "" );

                        // all other fixed variables in this constraint
                        for ( int k=0 ; (k<nz) ; ++k )
                        {
                            if (k!=j)
                            {
                                cprop_add_arc_impl_g( cprop, 
                                        lb[idx[k]]>=0.98 ? EOne : EZero, idx[k], 
                                        EOne, idx[j] );
                            }
                        }
                    }
                }
            } // negative
        } // not positive
    }

    return nimpl;
}

void cprop_report_infeasible_constraint( CProp *cprop, int irow, double minLhs )
{
    char strRDes[256];
    cprop_constraint_descr( cprop, irow, strRDes );
    char msg[1024];
    double rhs = cprop_rhs( cprop, irow );
    sprintf( msg, "Constraint %s impossible to satisfy (%g <= %g)\n", strRDes, minLhs, rhs );

    if (cprop->verbose)
        printf("%s\n", msg );

    cprop_add_msg_inf( cprop, msg );

    /* checking fixed variables */
    int nz = cprop_nz( cprop, irow );
    const int *idx = cprop_idx( cprop, irow );
    int nFixed = 0;
    const double *ub = cprop->ub;
    const double *lb = cprop->lb;
    for ( int j=0 ; j<nz ; ++j )
        if (FIXED(lb[idx[j]], ub[idx[j]]))
            ++nFixed;

    if ( nFixed == 0 )
        return;

    char strFix[512];
    cprop_constraint_descr_variable_bounds( cprop, irow, strFix );
    char strComment[256] = "";
    if (strstr(strFix,"*"))
        strcpy( strComment, " (* implied bound)" );
    sprintf( msg, "Fixed variables: %s\n", strFix );
    cprop_add_msg_inf( cprop, msg );

    if (cprop->verbose)
        printf("%s", msg );
    for ( int j=0 ; (j<nz) ; ++j )
    {
        if (FIXED( lb[idx[j]], ub[idx[j]]))
        {
            cprop_add_arc_impl_g( cprop, 
                    lb[idx[j]]>=0.98 ? EOne : EZero, idx[j], 
                    Infeasible, 0 );
        }
    }
}

/* considers constraint in for format ax <= b */
int cprop_process_constraint( CProp *cprop, int irow )
{
    double rhs = cprop_rhs( cprop, irow );

    /* rhs too large */
    if ( rhs >= oo )
        return 0;

    /* all fixed ? */
    int nz = cprop_nz( cprop, irow );
    const int *idx = cprop_idx( cprop, irow );
    const double *coef = cprop_coef( cprop, irow );
    double *ub = cprop->ub;
    double *lb = cprop->lb;
    int nFixed = 0;
    for ( int j=0 ; j<nz ; ++j )
        if (FIXED(lb[idx[j]], ub[idx[j]]))
            ++nFixed;

    if ( nFixed == nz  )
    {
        double lhs = 0.0;
        for ( int j=0 ; j<nz ; ++j )
            lhs += ub[idx[j]]*coef[j];

        if ( lhs >= rhs+VEPS )
            cprop_report_infeasible_constraint( cprop, irow, lhs );
    }

    switch ( vec_char_get(cprop->vrtype, irow) )
    {
        case ConstraintBV:
            return cprop_process_constraint_binary_variables( cprop, irow );
            break;
        case ConstraintINNB:
            return 0;
    }

    return 0;
}

int cprop_add_constraint( CProp *cprop, int nz, const int idx[], const double coef[], char sense, double rhs, const char rname[] )
{
    int res1 = 0, res2 = 0;

    if (!cprop->feasible)
    {
        fprintf( stderr, "Cannot add more constraints to infeasible program.\n" );
        abort();
    }

    Vec_int *vcolimpl = cprop->vcolimpl;

    int implBoundsStart = vec_int_size( vcolimpl );

    {
        double mult = toupper(sense) == 'G' || toupper(sense) == 'E' ? -1.0 : 1.0;
        int irow = cprop_n_rows( cprop );
        char rName[256];
        sprintf( rName, "%s%s", rname, toupper(sense) == 'E' ? "p1" : "" );
        cprop_add_row( cprop, nz, idx, coef, rhs, mult, rName );
        res1 = cprop_process_constraint( cprop, irow );
        if (cprop->verbose)
            printf("\n");

    }
    if (toupper(sense)=='E')
    {
        int irow = cprop_n_rows( cprop );
        char rName[256];
        sprintf( rName, "%s%s", rname, toupper(sense) == 'E' ? "p2" : "" );
        cprop_add_row( cprop, nz, idx, coef, rhs, 1.0, rName );
        res2 = cprop_process_constraint( cprop, irow );
        if (cprop->verbose)
            printf("\n");
    }

    cprop->nimpl = vec_int_size( vcolimpl ) - implBoundsStart;

    if ( res1==-1 || res2==-1 )
    {
        cprop->feasible = 0;
        return -1;
    }

    return res1 + res2;
}

void cprop_try_add_stk_rows( int nr, const int rows[], Vec_int *stkr, Vec_char *ivr )
{
    for ( int i=0 ; (i<nr) ; ++i )
    {
        if ( vec_char_get( ivr, rows[i] ) == 0 )
        {
            vec_int_push_back( stkr, rows[i] );
            vec_char_set( ivr, rows[i], 1 );
        }
    }
}


int cprop_update_bound( CProp *cprop, int j, double l, double u )
{
    if (!cprop->feasible)
    {
        fprintf( stderr, "Cannot change bounds on infeasible program. Use cprop_undo to undo last changes.\n" );
        abort();
    }

    cprop->nimpl = 0;

    double *lb = cprop->lb;
    double *ub = cprop->ub;

    /* stores the number of implications since a bound changes */
    Vec_int *vnimpl = cprop->vnimpl;    

    /* columns with bounds changed and their old bounds */
    Vec_int *vcolimpl = cprop->vcolimpl;
    Vec_double *voldlb = cprop->voldlb;
    Vec_double *voldub = cprop->voldub;

    int implBoundsStart = vec_int_size( vcolimpl );
    int nArcsIGStart = vec_IntPair_size( cprop->lastArcsIG );

    /* storing original bound of this column */
    vec_int_push_back( vcolimpl, j );
    vec_double_push_back( voldlb, lb[j] );
    vec_double_push_back( voldub, ub[j] );

    /* rows to be processed after implications */
    Vec_int *stkRows = cprop->stkRows;
    Vec_char *ivr = cprop->ivr;

    lb[j] = l;
    ub[j] = u;

    if ( vec_char_size( ivr ) < cprop_n_rows( cprop ) )
        vec_char_resize( ivr, cprop_n_rows( cprop), 0 );

#ifdef DEBUG
    for ( int i=0 ; (i<vec_char_size(ivr)) ; ++i )
        assert( vec_char_get(ivr,i) == 0);
#endif

    int implBefore = vec_int_size( vcolimpl );


    /* processing initial rows of the first column */
    {
        int nr = cprop_n_rows_col( cprop, j );
        const int *rcs = cprop_rows_col( cprop, j );

        for (int i=0 ; (i<nr) ; ++i )
        {
            int newImpl = cprop_process_constraint( cprop, rcs[i] );
            if ( newImpl == -1 )
            {
                cprop->feasible = 0;
                cprop->nimpl = vec_int_size( vcolimpl ) - implBoundsStart - 1;
                vec_int_push_back( vnimpl, vec_int_size( vcolimpl ) - implBoundsStart );
                vec_int_push_back( cprop->nNewArcsIG, vec_IntPair_size(cprop->lastArcsIG)-nArcsIGStart );
                return -1;
            }
        }
    }

    // processing constraints of columns with new implied bounds
    int nNewImpl;
    while ( (nNewImpl=(vec_int_size( vcolimpl ) - implBefore)) )
    {
        vec_int_clear( stkRows );

        for ( int i=implBefore ; i<implBefore+nNewImpl ; ++i )
        {
            // checking rows of this new column
            j = vec_int_get( vcolimpl, i );
            int nr = cprop_n_rows_col( cprop, j );
            const int *rcs = cprop_rows_col( cprop, j );
            cprop_try_add_stk_rows( nr, rcs, stkRows, ivr );
        }

        implBefore = vec_int_size( vcolimpl );

        int newImpl = 0;
        /* computed new set of rows, processing them */
        for ( int i=0 ; (i<vec_int_size(stkRows)) ; i++ )
        {
            newImpl = cprop_process_constraint( cprop, vec_int_get( stkRows, i ) );
            if ( newImpl == -1 )
                goto DONE_PROCESSING_ROWS;
        }
DONE_PROCESSING_ROWS:
        /* clearing incidence vector */
        for ( int l=0 ; (l<vec_int_size(stkRows)) ; ++l )
            vec_char_set( ivr, vec_int_get( stkRows, l ), 0 );

        if ( newImpl == -1 )
        {
            cprop->feasible = 0;
            cprop->nimpl = vec_int_size( vcolimpl ) - implBoundsStart - 1;
            vec_int_push_back( vnimpl, vec_int_size( vcolimpl ) - implBoundsStart );
            vec_int_push_back( cprop->nNewArcsIG, vec_IntPair_size(cprop->lastArcsIG)-nArcsIGStart );
            return -1;
        }
    }

    vec_int_push_back( vnimpl, vec_int_size( vcolimpl ) - implBoundsStart );
    vec_int_push_back( cprop->nNewArcsIG, vec_IntPair_size(cprop->lastArcsIG)-nArcsIGStart );
    cprop->nimpl = vec_int_size( vcolimpl ) - implBoundsStart - 1;
    return vec_int_size( vcolimpl ) - implBoundsStart - 1;
}

void cprop_add_row( CProp *cprop, int nz, const int idx[], const double coef[], double rhs, double mult, const char rname[] )
{
    rhs = MIN( rhs, oo );
    rhs = MAX( rhs, -oo );

    Vec_int *vrstart = cprop->vrstart;
    Vec_int *vrnz = cprop->vrnz;
    Vec_int *vridx = cprop->vridx;
    Vec_double *vrcoef = cprop->vrcoef;
    Vec_double *vrrhs = cprop->vrrhs;
    V2D_int *rowsCol = cprop->rowsCol;
    const char *integer = cprop->integer;
    const double *lb = cprop->lb; 
    const double *ub = cprop->ub;

    int start = vec_int_size( vrstart ) ? vec_int_last( vrstart ) +  vec_int_last( vrnz ) : 0;
    vec_int_push_back( vrstart, start );
    vec_int_push_back( vrnz, nz );

    /* number of binary, integer and unbounded variables in this constraint */
    int nBin = 0, nInt = 0, nNeg = 0, nPUnb = 0;
    for ( int j=0 ; (j<nz) ; ++j )
    {
        assert( idx[j] >= 0 && idx[0] < cprop->cols );
        vec_int_push_back( vridx, idx[j] );
        double v = MIN( coef[j], oo );
        v = MAX( -oo, v )*mult;
        vec_double_push_back( vrcoef, v );

        if (integer[idx[j]])
        {
            if (lb[idx[j]]>= -1e-10 && ub[idx[j]]<=1+1e+10)
                nBin++;
            else
            {
                nInt++;
                if (lb[idx[j]]<=-1e-10)
                    nNeg++;

                if (ub[idx[j]]>oo)
                    nPUnb++;
            }
        }
    }

    enum ConstraintType consType = ConstraintOther;
    if ( nBin == nz )
        consType = ConstraintBV;
    else
        if ( nInt==nz && nNeg==0 && nPUnb==0 )
            consType = ConstraintINNB;

    vec_char_push_back( cprop->vrtype, consType );

    vec_double_push_back( vrrhs, rhs*mult );
    for ( int j=0 ; j<nz ; ++j )
        v2d_int_row_push_back( rowsCol, idx[j], vec_int_size(vrnz)-1 );

    cprop_add_row_name( cprop, rname );
}

int cprop_nz( const CProp *cprop, int irow )
{
    return vec_int_get( cprop->vrnz, irow );
}

const int *cprop_idx( const CProp *cprop, int row )
{
    return vec_int_getp( cprop->vridx,
            vec_int_get( cprop->vrstart, row ) );
}

const double *cprop_coef( const CProp *cprop, int row )
{
    return vec_double_getp( cprop->vrcoef,
            vec_int_get( cprop->vrstart, row ) );
}

int cprop_n_rows( const CProp *cprop )
{
    return vec_int_size( cprop->vrnz );
}

double cprop_rhs( const CProp *cprop, int row )
{
    return vec_double_get( cprop->vrrhs, row );
}

int cprop_n_rows_col( const CProp *cprop, int col )
{
    return v2d_int_row_size( cprop->rowsCol, col );
}

int *cprop_rows_col( const CProp *cprop, int col )
{
    return v2d_int_row_ptr( cprop->rowsCol, col );
}

void cprop_undo( CProp *cprop )
{
    Vec_int *vnimpl = cprop->vnimpl;
    Vec_int *vcolimpl = cprop->vcolimpl;
    Vec_double *voldlb = cprop->voldlb;
    Vec_double *voldub = cprop->voldub;
    double *lb = cprop->lb;
    double *ub = cprop->ub;

    assert( vec_int_size( vnimpl ) );

    int n = vec_int_last( vnimpl );

    assert( n >= 1);

    for ( int i=0 ; (i<n) ; ++i )
    {
        int j = vec_int_pop_back( vcolimpl );
        lb[j] = vec_double_pop_back( voldlb );
        ub[j] = vec_double_pop_back( voldub );
    }

    vec_int_pop_back( vnimpl );

    int narcs = vec_int_pop_back( cprop->nNewArcsIG );
    for ( int i=0 ; i<narcs ; ++i )
    {
        const IntPair arc = vec_IntPair_pop_back( cprop->lastArcsIG );
        iset_remove( cprop->implGIn[arc.a], arc.b );
    }

    cprop->feasible = 1;
}

void cprop_add_row_name( CProp *cprop, const char rname[] )
{
    vec_int_push_back( cprop->vrnamest, vec_char_size(cprop->vrnames) );
    int l = strlen( rname );
    for ( int i=0 ; (i<l) ; ++i )
        vec_char_push_back( cprop->vrnames, rname[i] );
    vec_char_push_back( cprop->vrnames, '\0' );
}

const char *cprop_row_name( const CProp *cprop, int i )
{
    return vec_char_getp( cprop->vrnames, vec_int_get( cprop->vrnamest, i ) );
}

char cprop_feasible(const CProp* cprop)
{
    return cprop->feasible;
}

int cprop_n_implications( const CProp *cprop )
{
    return cprop->nimpl;
}

int cprop_implied_var( const CProp *cprop, int i )
{
    if (cprop->nimpl==0)
    {
        fprintf( stderr, "No implications were generated in the last operation.\n" );
        abort();

    }
    return vec_int_get( cprop->vcolimpl, vec_int_size(cprop->vcolimpl)-cprop_n_implications(cprop)+i  );
}

double cprop_get_lb( const CProp *cprop, int j )
{
    assert( j>=0 && j<cprop->cols );
    return cprop->lb[j];
}

double cprop_get_ub( const CProp *cprop, int j )
{
    assert( j>=0 && j<cprop->cols );
    return cprop->ub[j];
}


void cprop_add_arc_impl_g( CProp *cprop, enum IGNType ntSource, int colSource, enum IGNType ntDest, int colDest )
{
    int sIndex = cprop_impl_graph_node_id( cprop, ntSource, colSource );
    int dIndex =  cprop_impl_graph_node_id( cprop, ntDest, colDest );

    if (cprop->implGIn[dIndex]==NULL)
        cprop->implGIn[dIndex] = iset_create( HASH_SIZE_IMPLG );

    iset_add( cprop->implGIn[dIndex], sIndex );

    IntPair arc = { dIndex, sIndex };

    vec_IntPair_push_back( cprop->lastArcsIG, arc );
}

enum IGNType cprop_impl_graph_node_type( const CProp *cprop, int nodeId )
{
    if ( nodeId == cprop->cols*2 )
        return Infeasible;

    if ( nodeId >= cprop->cols )
        return EOne;

    return EZero;
}

int cprop_impl_graph_node_var( const CProp *cprop, int nodeId )
{
    if ( nodeId == cprop->cols*2 )
    {
        fprintf( stderr, "Infeasible node is not related to a variable.\n" );
        abort();
    }

    if ( nodeId >= cprop->cols )
        return nodeId - cprop->cols;

    return nodeId;
}

int cprop_impl_graph_node_id( const CProp *cprop, enum IGNType ntype , int col )
{
    switch (ntype)
    {
        case Infeasible:
            return cprop->cols*2;
            break;
        case EOne:
            return cprop->cols + col;
            break;
        case EZero:
            return col;
            break;
    }

    return col; // keep compilers calm
}

int cprop_impl_graph_in_neigh( const CProp *cprop, int nodeId, int i )
{
    if ( cprop->implGIn[nodeId]==NULL )
    {
        fprintf( stderr, "Node %d has no neighbors in implication graph.\n", nodeId );
        abort();
    }

    return iset_element( cprop->implGIn[nodeId], i );
}

void cprop_set_verbose( CProp *cprop, char verbose )
{
    cprop->verbose = verbose;
}


void cprop_add_msg_inf( CProp *cprop, const char *msg )
{
    if ( cprop->msgInf == NULL )
        cprop->msgInf = vec_char_create();

    int l = strlen(msg);
    for ( int i=0 ; (i<l) ; ++i )
        vec_char_push_back( cprop->msgInf, msg[i] );
    vec_char_push_back(cprop->msgInf, '\0');
}

const char *cprop_inf_msg( const CProp *cprop )
{
    if (cprop->feasible)
    {
        fprintf( stderr, "Infeasibility not detected.");
        abort();
    }

    return vec_char_getp( cprop->msgInf, 0 );
}


char *cprop_constraint_descr( const CProp *cprop, int irow, char *descr )
{
    char strRow[1024] = "";
    int rnz = cprop_nz( cprop, irow );
    const int *idx = cprop_idx( cprop, irow );
    const double *coef = cprop_coef( cprop, irow );
    double rhs = cprop_rhs( cprop, irow );

    int dnz = MIN( rnz, MAX_VAR_ROW_PRINT );
    for ( int j=0 ; (j<dnz) ; ++j )
    {
        char strCol[256];
        sprintf( strCol, "%+g %s ", coef[j], cprop->cname[idx[j]] );
        strcat( strRow, strCol );
    }

    char strMore[16] = "";
    if (dnz<rnz)
        strcpy( strMore, " ..." );

    sprintf( descr, "%s: %s%s <= %g", cprop_row_name(cprop, irow), strRow, strMore, rhs );

    return descr; 
}

char *cprop_constraint_descr_variable_bounds( const CProp *cprop, int irow, char *descr )
{
    int nz = cprop_nz( cprop, irow );
    int printed = 0;
    const int *idx = cprop_idx( cprop, irow );
    const double *lb = cprop->lb;
    const double *ub = cprop->ub;
    const double *olb = cprop->olb;
    const double *oub = cprop->oub;
    const char **cname = (const char**) cprop->cname;

    char strFixed[1024] = "";

    /* first printing variables fixed */
    for ( int j=0 ; (j<nz ) ; ++j )
    {
        if (printed==MAX_VAR_ROW_PRINT)
        {
            strcat( strFixed, " ... " );
            break;
        }

        if (FIXED(lb[idx[j]], ub[idx[j]]))
        {
            char str[256];
            char stri[2] = "";
            if (!FIXED(olb[idx[j]], oub[idx[j]]))
                strcpy( stri, "*" );

            sprintf( str, "%s%s=%g ", stri, cname[idx[j]], lb[idx[j]] );
            strcat( strFixed, str );
            ++printed;
        }
    }
    strcpy( strFixed, descr );

    return descr;
}

void cprop_save_impl_graph( const CProp *cprop, const char *fName )
{
    FILE *f = fopen( fName, "w" );
    assert( f );

    fprintf( f, "digraph ImplicationGraph {\n" );

    const int nNodes = cprop->cols*2+1;
    for ( int i=0 ; i<nNodes ; ++i )
    {
        const ISet *in = cprop->implGIn[i];
        if (in == NULL)
            continue;

        char strDest[256] = "";

        enum IGNType nDestType = cprop_impl_graph_node_type( cprop, i );
        if (nDestType == Infeasible)
            strcpy( strDest, "Infeasible" );
        else
        {
            int ncol = cprop_impl_graph_node_var( cprop, i );

            if ( nDestType == EOne )
                sprintf( strDest, "%sa1", cprop->cname[ncol] );
            else
                sprintf( strDest, "%sa0", cprop->cname[ncol] );
        }

        for ( int j=0 ; (j<iset_n_elements(in)) ; ++j )
        {
            int source = iset_element( in, j );

            char strSource[256] = "";

            enum IGNType nSourceType = cprop_impl_graph_node_type( cprop, source );

            if (nSourceType==Infeasible)
                strcpy( strSource, "Infeasible");
            {
                int ncol = cprop_impl_graph_node_var( cprop, source );

                if ( nSourceType == EOne )
                    sprintf( strSource, "%sa1", cprop->cname[ncol] );
                else
                    sprintf( strSource, "%sa0", cprop->cname[ncol] );
            }

            fprintf( f, "\t%s -> %s;\n", strSource, strDest );
        }
    }

    fprintf( f, "}\n" );

    fclose( f );
}

void cprop_free( CProp **cprop )
{
    free( (*cprop)->lb );
    free( (*cprop)->ub );
    free( (*cprop)->olb );
    free( (*cprop)->oub );
    free( (*cprop)->integer );
    free( (*cprop)->binary );
    if ( (*cprop)->cname )
    {
        free( (*cprop)->cname[0] );
        free( (*cprop)->cname );
    }

    free( (*cprop)->idx );
    free( (*cprop)->coef );
    free( (*cprop)->pv );
    free( (*cprop)->nv );
    vec_char_free( &(*cprop)->ivr );
    vec_int_free( &(*cprop)->stkRows );

    vec_int_free( &(*cprop)->vrstart );
    vec_int_free( &(*cprop)->vrnz );
    vec_int_free( &(*cprop)->vridx );
    vec_double_free( &(*cprop)->vrcoef );
    vec_double_free( &(*cprop)->vrrhs );
    vec_char_free( &(*cprop)->vrtype);
    v2d_free( &(*cprop)->rowsCol );

    vec_int_free( &((*cprop)->vnimpl) );
    vec_int_free( &((*cprop)->vcolimpl) );
    vec_double_free( &((*cprop)->voldlb) );
    vec_double_free( &((*cprop)->voldub) );

    vec_int_free( &((*cprop)->vrnamest) );
    vec_char_free( &((*cprop)->vrnames) );

    vec_IntPair_free( &((*cprop)->lastArcsIG ) );
    vec_int_free( &((*cprop)->nNewArcsIG ) );

    for ( int j=0 ; (j<2*(*cprop)->cols+1) ; ++j )
        if (((*cprop)->implGIn[j]))
            iset_free( &((*cprop)->implGIn[j]) );
    free( (*cprop)->implGIn );

    if ((*cprop)->msgInf)
        vec_char_free( &((*cprop)->msgInf) );

    free( *cprop );
    *cprop = NULL;
}

