#ifndef MACROS_H_INCLUDED
#define MACROS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define INT_RANDOM( n ) \
   ((int) floor( ((double)(n)) * (((double)rand())/(((double)RAND_MAX)+((double)1.0))) ))

/* generates a random number between 0 and 1 */
#define DBL_RANDOM( ) \
   ((double) rand()) / (((double)RAND_MAX))

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

#define FILE_NAME_SIZE 1024
#define STR_SIZE        256
#define LSTR_SIZE       512
#define LINE_SIZE      2048

#define ALLOCATE( ptr, type ) {\
    ptr = (type*) malloc( sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_INI( ptr, type ) {\
    ptr = calloc( 1, sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_VECTOR( ptr, type, nElements ) {\
    ptr = (type*) malloc( sizeof(type)*(nElements) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

/* allocate filling with zeros */
#define ALLOCATE_VECTOR_INI( ptr, type, nElements ) {\
    ptr = (type*) calloc( (nElements), sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }


#define OPEN_FILE( f, fileName, mode )  \
    f = fopen( fileName, mode ); \
    if (!f) { \
        fflush(stdout); \
        fprintf( stderr, "ERROR: could not open file %s with mode %s. At %s:%d\n", fileName, mode, __FILE__, __LINE__ ); \
        fflush(stderr);  abort(); exit(EXIT_FAILURE); }

#define True  1
#define False 0

/* fills from start until the last element before end */
#define FILL( vector, start, end, value ) { \
    int i; \
    for ( i=start ; (i<end) ; ++i ) vector[i] = value; \
} \

#define EPS 1e-5

#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2)<=EPS )

#endif

