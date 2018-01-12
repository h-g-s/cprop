#include <stdio.h>
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

#define EXEC_AND_STORE_TIME( code, seconds ) ( { clock_t start = clock(); code; clock_t end = clock(); seconds = ( ((double)end-start)/((double)CLOCKS_PER_SEC) ); } )


int main( int argc, char **argv )
{
    if (argc<2)
    {
        fprintf( stderr, "Enter instance name.\n" );
        exit( EXIT_FAILURE );
    }

    getFileName( probName, argv[1] );
    
    /* now, setting a function to handle errors */
    conectSignals();

    mip = lp_create();

    EXEC_AND_STORE_TIME( lp_read( mip, argv[1] ), secread );

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
        
    EXEC_AND_STORE_TIME( cprop = cprop_create_from_mip( mip, 0 ), cpsecs );    

    feasibleCPROP = cprop_feasible( cprop );

    if (!cprop)
    {
        strcpy( msgError, "Could not create CPROP.\n" );
        abort();
    }

    printf("CPROP preprocessing took %.4f seconds.\n", cpsecs );
    
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
        fprintf( flog, "instance,cols,rows,nz,feasibleLP,secLPOpt,objLP,cpropFeas,cpropTime,ppLPTime,ppCols,ppRows,ppNZ,ppLPOptTime,ppLPObj,ppLPFeasible,error\n" );
                 //  1  2  3  4  5  6    7  8  9    10  11 12 13  14  15 16 18 
    fprintf( flog, "%s,%d,%d,%d,%d,%.4f,%g,%d,%.4f,%.4f,%d,%d,%d,%.4f,%g,%d,\"%s\"\n", 
      //   1        2     3    4       5         6        7         8           9       10      11      12     13      14        15           16          17
        probName, cols, rows, nz, feasibleLP, seclpopt, objlp, feasibleCPROP, cpsecs, ppsecs, ppcols, pprows, ppnz, secpplpopt, objpplp, feasiblePPLP, msgError );

    
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

