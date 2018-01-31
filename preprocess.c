#define _GNU_SOURCE

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <signal.h>
#include "macros.h"
#include "memory.h"
#include "lp.h"
#include "cprop.h"
#include "containers.h"
#include "strutils.h"
#include "cprop_lp.h"
#include "cp_cuts.h"



void conectSignals();
void exitPP( int sig );

/* output information */
LinearProgram *mip = NULL, *ppmip = NULL;
CProp *cprop = NULL;
char probName[256] = "";
char feasibleLP = 0, feasibleCPROP = 0, feasiblePPLP = 0;
double secread = 0, seclpopt = 0, cpsecs = 0, ppsecs = 0, secpplpopt = 0;
double objlp = DBL_MAX, objpplp = DBL_MAX;
char msgError[1024] = "";
int cols = 0, rows = 0, nz = 0;
int ppcols = 0, pprows = 0, ppnz = 0;
int *origCols = NULL;

int totalCutsPC = 0; // cut in probe and cut
double objPC = DBL_MAX;
double secsPC = 0.0;
int roundsPC = 0;

char probeAndCut = False;
char probCOnlyFrac = True;
int probCMaxTries = INT_MAX;
enum DVDir probCDir = DVDBoth;

#define EXEC_AND_STORE_TIME( code, seconds ) ( { clock_t start = clock(); code; clock_t end = clock(); seconds = ( ((double)end-start)/((double)CLOCKS_PER_SEC) ); } )

/* parameters */
char printSol = 0;

#define MAX_DIVE_DEPT 5
char diveVar[MAX_DIVE_DEPT][256];
int diveVal[MAX_DIVE_DEPT];
int nVarsToDive = 0;

char verbose = 1;

void parseParameters( int argc, const char **argv );

int main( int argc, char **argv )
{
    if (argc<2)
    {
        fprintf( stderr, "Enter instance name.\n" );
        exit( EXIT_FAILURE );
    }

    getFileName( probName, argv[1] );
    for ( int i=0 ; (i<MAX_DIVE_DEPT) ; ++i )
    {
        strcpy( diveVar[i], "" );
        diveVal[i] = 1;
    }
    parseParameters( argc, (const char**) argv );
    
    /* now, setting a function to handle errors */
    conectSignals();

    mip = lp_create();

    EXEC_AND_STORE_TIME( lp_read( mip, argv[1] ), secread );
    
    /* debug */
//    {
//        int ridx = lp_row_index( mip, "R4341" );
//        assert( ridx != -1 );
//        int *idx; ALLOCATE_VECTOR( idx, int, lp_cols(mip) );
//        double *coef; ALLOCATE_VECTOR( coef, double, lp_cols(mip) );
//        int rnz = lp_row( mip, ridx, idx, coef );
//        for ( int i=0 ; (i<rnz) ; ++i )
//        {
//            int cidx = idx[i];
//            char cname[256]; lp_col_name( mip, cidx, cname );
//            double clb = lp_col_lb( mip, cidx );
//            double cub = lp_col_ub( mip, cidx );
//            char isint = lp_is_integer( mip, cidx );
//            printf("%s [%g..%g] int %d\n", cname, clb, cub, isint );
//        }
//        
//        free( idx ); free( coef );
//    }

    cols = lp_cols( mip );
    rows = lp_rows( mip );
    nz = lp_nz( mip );
    printf("\n%s Read in %.4f seconds\n", probName, secread );
    printf("Original MIP has dimensions %d/%d/%d (cols/rows/nz)\n", cols, rows, nz );

    int status;
    EXEC_AND_STORE_TIME( status = lp_optimize_as_continuous( mip ), seclpopt );

    if (status!=LP_OPTIMAL)
    {
        feasibleLP = 0;

        if (status!=LP_INFEASIBLE)
        {
            fprintf( stderr, "Error while optimizing. Status %d\n", status );
            exit( EXIT_FAILURE );
        }
    }
    else
    {
        feasibleLP = 1;
        objlp = lp_obj_value( mip );
        printf("Initial LP relaxation solved in %.4f seconds, obj is %g\n", seclpopt, lp_obj_value(mip) );
    } 
        
    EXEC_AND_STORE_TIME( cprop = cprop_create_from_mip( mip, verbose ), cpsecs );    
    
    feasibleCPROP = cprop_feasible( cprop );

    if (!feasibleCPROP)
    {
        printf("Problem is infeasible.\n\t%s\n", cprop_inf_msg(cprop) );
        exitPP(0);
    }
    else
        cprop_enter_relax_sol( cprop, lp_x(mip) );


    if (!cprop)
    {
        strcpy( msgError, "Could not create CPROP.\n" );
        abort();
    }

    printf("CPROP preprocessing took %.4f seconds. Fixed %d columns\n", cpsecs, cprop_n_fixed_at_pre_proc(cprop) );
    
    if ( cprop_n_fixed_at_pre_proc(cprop) )
    {
        const char *fixed = cprop_fixed_at_pre_proc( cprop );
        for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
        {
            if (!fixed[i])
                continue;

            lp_set_col_bounds( mip, i, cprop_get_lb(cprop,i), cprop_get_lb(cprop,i) );
        }
        int status = lp_optimize_as_continuous( mip );
        assert( status == LP_OPTIMAL );
        cprop_enter_relax_sol( cprop, lp_x(mip) );

        if (printSol)
        {
            //printf()
            const double *x = lp_x( mip );
            for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
            {
                if (fabs(x[i])<1e-10)
                    continue;

                char cName[256];
                printf("%40s %g\n", lp_col_name( mip, i, cName ), x[i] );
            } // all variables
        } // print pp sol
        
        if (nVarsToDive)
        {
            printf("\nDiving:\n");
            // checking if diving in some variable should be performed
            for ( int i=0 ; (i<nVarsToDive) ; ++i )
            {
                int idxvtd = lp_col_index( mip, diveVar[i] );
                char cName[256]; 
                assert( strcmp( diveVar[i], lp_col_name( mip, idxvtd, cName ) ) == 0 );
                
                int idxorigvar = lp_col_index( mip, diveVar[i] );
                int res = cprop_update_bound( cprop, idxorigvar, diveVal[i], diveVal[i] );
                printf("\t%s=%d ... ", cName, diveVal[i] );
                switch (res)
                {
                    case -1:
                    {
                        const CPCuts *cp = cprop_cut_pool( cprop );
                        /* printing */
                        printf(" infeasible. %d cuts implied.\n", cpc_n_cuts(cp) );
                        cprop_save_impl_graph( cprop, "impl.dot" );
                        for ( int ic=0 ; (ic<cpc_n_cuts(cp)) ; ++ic )
                        {
                            double lhs = 0;
                            int nz = cpc_nz( cp, ic );
                            const int *idx = cpc_idx( cp, ic );
                            const double *coef = cpc_coef( cp, ic );
                            for ( int ii=0 ; (ii<nz) ; ++ii )
                                lhs += coef[ii] * lp_x(mip)[idx[ii]];

                            if ( lhs - cpc_rhs(cp,ic) < 1e-10 )
                                continue;
                            printf("\t");
                            for ( int ii=0 ; (ii<nz) ; ++ii )
                            {
                                printf("%+g %s ", coef[ii], lp_col_name(mip, idx[ii], cName) );                                
                            }
                            printf("<= %+g (viol: %g)\n", cpc_rhs( cp, ic ), lhs - cpc_rhs(cp,ic) );
                            
                            static int iCut = 0; char cutName[256]; sprintf( cutName, "cut%d", iCut++ );
                            lp_add_row( mip, nz, (int*) idx, (double*) coef,  cutName, 'L', cpc_rhs(cp,ic) );
                        }
                        
                        if (cpc_n_cuts(cp))
                        {
                            int status = lp_optimize_as_continuous( mip );
                            assert( status == LP_OPTIMAL );
                            printf("After cuts obj is %g\n", lp_obj_value(mip) );
                            lp_write_lp( mip, "cuts.lp");
                        }

                        break;
                    }                
                    case 0:
                    {
                        printf("no implications.\n");
                        break;
                    }
                    default:
                    {
                        printf("%d implications.\n\t", cprop_n_implications(cprop) );
                        assert( cprop_n_implications(cprop)>=1 );
                        assert( res==cprop_n_implications(cprop) );
                        for ( int j=0 ; (j<cprop_n_implications(cprop)) ; ++j )
                        {
                            int origIdx = cprop_implied_var(cprop,j);
                            lp_col_name( mip, origIdx, cName );
                            printf("%s=%d ", cName, (int)cprop_get_lb(cprop, origIdx) );
                        }
                        if (cprop_n_implications(cprop))
                            printf("\n");
                    }
                }
            }
        } // if diving should be performed 

        cprop_clear( cprop );
    }
    
    ALLOCATE_VECTOR( origCols, int, cols );
    
    EXEC_AND_STORE_TIME( ppmip = cprop_preprocess( cprop, mip, 1, origCols ), ppsecs );
    
    ppcols = lp_cols( ppmip );
    pprows = lp_rows( ppmip );
    ppnz = lp_nz( ppmip );

    printf("Pre-processed MIP has dimensions %d/%d/%d (cols/rows/nz)\n", ppcols, pprows, ppnz );
    
    EXEC_AND_STORE_TIME( status = lp_optimize_as_continuous( ppmip ), secpplpopt );
    if (status!=LP_OPTIMAL)
    {
        feasiblePPLP = 0;

        if (status!=LP_INFEASIBLE)
        {
            fprintf( stderr, "Error while optimizing. Status %d\n", status );
            exit( EXIT_FAILURE );
        }
    }
    else
    {
        feasiblePPLP = 1;
        objpplp = lp_obj_value( ppmip );
        printf("Initial LP relaxation of Pre-Processed problem solved in %.4f seconds, obj is %g\n", secpplpopt, lp_obj_value(ppmip) );

    }

    if (probeAndCut)
    {
        clock_t startpc = clock();
TRY_TO_CUT:
        {
            const CPCuts *cpCuts = cprop_cut_pool( cprop );
            int nCuts = cpc_n_cuts( cpCuts );
            printf("Starting Probe and Cut, round %d\n", ++roundsPC );
            int newCuts = cprop_probe_and_cut( cprop, lp_x(mip), probCMaxTries, DVCMostFrac, probCOnlyFrac, probCDir  );
            printf("%d new cuts found.\n", newCuts );

            // new cuts
            for ( int i=nCuts ; (i<newCuts) ; ++i )
            {
                int nz = cpc_nz( cpCuts, i );
                int *idx = cpc_idx( cpCuts, i );
                double *coef = cpc_coef( cpCuts, i );
                double rhs = cpc_rhs( cpCuts, i );
                char cutName[256]=""; static int cutId = 0;
                sprintf( cutName, "cut%08d", cutId++ ); 
                lp_add_row( mip, nz, idx, coef, cutName, 'L', rhs );
            }



            if (newCuts>nCuts)
            {
                totalCutsPC += newCuts-nCuts;
                int status = lp_optimize_as_continuous( mip );
                lp_write_lp( mip, "probc.lp");
                
                assert( status == LP_OPTIMAL );
                if ( status == LP_OPTIMAL )
                {
                    objPC = lp_obj_value( mip );
                    goto TRY_TO_CUT;
                }
                else
                {
                    objPC = DBL_MAX;
                }
            }
        }
        

        secsPC = ( (double) clock()-startpc ) / ((double)CLOCKS_PER_SEC);
    }

    // dive propagate and cut
    if (feasiblePPLP)
    {
        /*
        printf("diving into fractional solution\n");
        double *x = lp_x(ppmip);
        Vec_IntPair *vip = vec_IntPair_create();

        for ( int i=0 ; (i<lp_cols(ppmip)) ; ++i )
        {
            double intdist = MIN( ceilf(x[i])-x[i], x[i]-floor(x[i]) );
            IntPair ip = { i, (int)(intdist*100.0) };
            vec_IntPair_push_back( vip, ip );
        }

        vec_IntPair_sort( vip );


        for ( int i=0 ; (i<lp_cols(ppmip)) ; ++i )
        {
            char cname[256];
            int col = vec_IntPair_get( vip, lp_cols(ppmip)-i-1 ).a;
            int v = floor( x[col]+0.5 );
            lp_col_name( ppmip, col, cname );

            printf("fixing %s = %d  (%g) ... ", cname, v, x[col] );

            int idxOrig = lp_col_index( mip, cname );

            if (cprop_fixed(cprop, idxOrig))
            {
                printf("already fixed\n");
            }
            else
            {
                int res = cprop_update_bound( cprop, idxOrig, v, v );
                if ( res == -1 )
                {
                    const CPCuts *cp = cprop_cut_pool( cprop );
                    printf(" infeasible. %d cuts implied.\n", cpc_n_cuts(cp) );
                    break;
                }
                else
                    printf("%d implications.\n", cprop_n_implications(cprop) );
            }

        }

        vec_IntPair_free(&vip); */
    }
    
    exitPP( 0 );
}

void conectSignals()
{
    /* processing stopping due to some error */
    signal( SIGHUP, exitPP ); // user's terminal is disconnected
    signal( SIGABRT, exitPP ); // process detects error and reports by calling abort
    signal( SIGSEGV, exitPP ); // segmentation violation
    signal( SIGSTKFLT, exitPP ); // stack fault
    signal( SIGBUS, exitPP ); // access to an invalid address
    signal( SIGFPE, exitPP ); // floating point exception
    signal( SIGILL, exitPP ); // illegal instruction (usually a corrupted executable)
    signal( SIGSYS, exitPP ); // bad system call
    signal( SIGXFSZ, exitPP ); // file size limit exceeded
    signal( SIGXCPU, exitPP ); // CPU limit exceeded
    signal( SIGPIPE, exitPP ); // broken pipe
    
    /* process stopped by the user */
    signal( SIGQUIT, exitPP ); // terminate process and generate core dump
    signal( SIGINT, exitPP ); // interruptec (ctrl+c)
    signal( SIGTERM, exitPP ); // generated by "kill" command
}

void checkSignal( int signal )
{
    char spc[3] = "";
    if (strlen(msgError))
        strcpy( spc, " " );
    switch (signal)
    {
        case SIGHUP:
            sprintf( msgError, "%suser's terminal is disconnected.", spc );
            break;
        case SIGABRT:
            sprintf( msgError, "%sabort called.", spc );
            break;
        case SIGSEGV:
            sprintf( msgError, "%ssegmentation violation.", spc );
            break;
        case SIGSTKFLT:
            sprintf( msgError, "%sstack fault.", spc );
            break;
        case SIGBUS:
            sprintf( msgError, "%saccess to an invalid address.", spc );
            break;
        case SIGFPE:
            sprintf( msgError, "%sfloating point exception.", spc );
            break;
        case SIGILL:
            sprintf( msgError, "%sillegal instruction (usually a corrupted executable).", spc );
            break;
        case SIGSYS:
            sprintf( msgError, "%sbad system call.", spc );
            break;
        case SIGXFSZ:
            sprintf( msgError, "%sfile size limit exceeded.", spc );
            break;
        case SIGXCPU:
            sprintf( msgError, "%sCPU limit exceeded.", spc );
            break;
        case SIGPIPE:
            sprintf( msgError, "%sbroken pipe.", spc );
            break;
    }

}

void exitPP( int sig )
{
    FILE *flog = NULL;
    
    checkSignal( sig );

    if (strlen(msgError))
    {
        fflush( stdout );  fflush( stderr ); 
        fprintf( stderr, "\n\nERROR: %s\n\n", msgError );
        fflush( stderr ); 
    }


    // checking if this is the first line
    flog = fopen( "summary.csv", "r" );
    char firstLine = (flog == NULL);
    if (!firstLine)
    {
        fclose( flog );
        flog = NULL;
    }
    
    flog = fopen( "summary.csv", "a" );
    assert( flog );

    if (firstLine)    //    1     2    3   4      5          6       7     8          9        10      11     12    13       14        15        16       17
        fprintf( flog, "instance,cols,rows,nz,feasibleLP,secLPOpt,objLP,cpropFeas,cpropTime,ppLPTime,ppCols,ppRows,ppNZ,ppLPOptTime,ppLPObj,ppLPFeasible,roundspc,cutspc,secspc,objpc,error\n" );
                 //  1  2  3  4  5  6    7  8  9    10  11 12 13  14  15 16 18 
    fprintf( flog, "%s,%d,%d,%d,%d,%.4f,%g,%d,%.4f,%.4f,%d,%d,%d,%.4f,%g,%d,%d,%d,%.2f,%g,\"%s\"\n", 
      //   1        2     3    4       5         6        7         8           9       10      11      12     13      14        15           16          17
        probName, cols, rows, nz, feasibleLP, seclpopt, objlp, feasibleCPROP, cpsecs, ppsecs, ppcols, pprows, ppnz, secpplpopt, objpplp, feasiblePPLP, roundsPC, totalCutsPC, secsPC, objPC, msgError );

    
    fclose(flog);

    if (mip)
        lp_free( &mip );
    if (cprop)
        cprop_free( &cprop );
    if (ppmip)
    {
        lp_write_lp( ppmip, "pp.lp" );
        lp_free( &ppmip );
    }
    if (origCols)
        free( origCols );    
    
    exit( strlen( msgError ) ? EXIT_FAILURE : EXIT_SUCCESS );
}

void parseParameters( int argc, const char **argv )
{
    for ( int i=1 ; (i<argc) ; ++i )
    {
        if (strlen(argv[i])<=1)
            continue;
        if (argv[i][0]=='-')
        {
            const char *param = &argv[i][1];
            if (strcasestr( param, "printSol" ))
            {
                printSol = 1;
                continue;
            }
            
            if (strcasestr(param, "diveV"))
            {
                assert( i+1<argc );
                strcpy( diveVar[nVarsToDive++], argv[i+1] );
                ++i;                
                continue;
            }

            if (strcasestr(param, "probeAndCut"))
            {
                probeAndCut = True;
                continue;
            }
            
            if (strcasestr(param, "verbose"))
            {
                assert( i+1<argc );
                verbose = atoi( argv[i+1] );
                ++i;                
                continue;
            }
                
        } // flags
    } // all parameters
}

