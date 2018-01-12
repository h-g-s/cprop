#ifndef CONTAINERS_H_INCLUDED
#define CONTAINERS_H_INCLUDED

/* set of generic and type safe containers written in ANSI C
 * can be compiled in to modes:
 *    default: fast
 *    with -DDEBUG: includes lots of checkings for easing debug
 *    c 2014, Haroldo Gambini Santos - haroldo@iceb.ufop.br */

#include <limits.h>
#include "macros.h"

typedef struct
{
    int a;
    int b;
} IntPair;

#define NOT_FOUND INT_MAX

/* to be used in qsort */
int intcmp_func( const void *p1, const void *p2 );

/* set of integers */
typedef struct _ISet ISet;
ISet *iset_create( int hashSize );
void iset_ini( ISet *iset, int hashSize );
void iset_cpy( ISet *target, const ISet *source );
void iset_add( ISet *set, int value );
void iset_remove( ISet *set, int value );
void iset_remove_set( ISet *set, const ISet *toRemove );
void iset_clear( ISet *set );
char iset_has( const ISet *set, int value );
int iset_n_elements( const ISet *set );
int iset_element( const ISet *set, int idx );
char iset_equals( const ISet *set1, const ISet *set2 );
void iset_clear_mem( ISet *iset );
void iset_free( ISet **_iset );

/* one dimensional vector */
#define DECLARE_VECTOR_TYPE(type) \
   typedef struct _##Vec_##type{ \
      type *array; \
      int size; \
      int capacity; \
   } Vec_##type; \
Vec_##type *vec_##type##_create( ); \
Vec_##type *vec_##type##_create_cap( const int iniCap ); \
Vec_##type *vec_##type##_create_fill( const int size, const type value ); \
int vec_##type##_size( const Vec_##type *vec ); \
char vec_##type##_equals( const Vec_##type *vec1, const Vec_##type *vec2 ); \
int vec_##type##_capacity( const Vec_##type *vec ); \
void vec_##type##_clear( Vec_##type *vec ); \
void vec_##type##_free( Vec_##type **vec ); \
void vec_##type##_push_back( Vec_##type *vec, const type value ); \
void vec_##type##_push_back_v( Vec_##type *vec, int n, const type v[] ); \
type vec_##type##_pop_back( Vec_##type *vec ); \
/* sets a value into an existing position */\
void vec_##type##_set( Vec_##type *vec, const int pos, const type value ); \
void vec_##type##_swap( Vec_##type *vec, const int pos1, const int pos2 ); \
/* sets a value into a position, if this position does not exists, creates */\
void vec_##type##_fset( Vec_##type *vec, const int pos, const type value ); \
type vec_##type##_first( const Vec_##type *vec ); \
type vec_##type##_last( const Vec_##type *vec ); \
type vec_##type##_get( const Vec_##type *vec, const int pos ); \
/* returns a pointer to a element */ \
type *vec_##type##_getp( Vec_##type *vec, const int pos ); \
/* returns a pointer to the vector */ \
type *vec_##type##_ptr( const Vec_##type *vec ); \
void vec_##type##_cpy( Vec_##type *vecTarget,  const Vec_##type *vecSource ); \
void vec_##type##_resize( Vec_##type *vec, const int size, const type valueNewElements  ); \
int vec_##type##_find( const Vec_##type *vec, const type value ); \
/* removes element with value, changes the order of the vector since the last element will go to the position of the removed element */ \
void vec_##type##_remove(Vec_##type *vec, const type value ); \
void vec_##type##_remove_positions( Vec_##type *vec, int nPos, const int positions[] ); \
void vec_##type##_insert_unique( Vec_##type *vec, const type value );

/* dictionary which maps strings to another type */
#define DECLARE_DICT_TYPE(type) \
    typedef struct { \
        char str[STR_SIZE]; \
        type value; \
        int keyPos; /* in keys vector */ \
    } Dict_Bucket_##type; \
    typedef struct { \
        Dict_Bucket_##type **cont; \
        int *rowSize; \
        int *rowCap; \
        type defValue; \
        unsigned int hashSize; \
        StrV *keys; \
    } Dict_##type; \
    Dict_##type *dict_##type##_create( unsigned int hashSize, type defaultValue ); \
    void dict_##type##_set( Dict_##type *dict, const char *key, const type value ); \
    void dict_##type##_iset( Dict_##type *dict, const int key, const type value ); \
    void dict_##type##_remove( Dict_##type *dict, const char *key ); \
    void dict_##type##_iremove( Dict_##type *dict, int key ); \
    type dict_##type##_get( const Dict_##type *dict, const char *str ); \
    int  dict_##type##_size( const Dict_##type *dict ); \
    const char *dict_##type##_key( const Dict_##type *dict, int pos ); \
    int dict_##type##_ikey( const Dict_##type *dict, int pos ); \
    type dict_##type##_iget( const Dict_##type *dict, const int key ); \
    void dict_##type##_clear( Dict_##type *dict ); \
    void dict_##type##_cpy( Dict_##type *target, const Dict_##type *source ); \
    void dict_##type##_free( Dict_##type **dict );

/* string vector */
typedef struct _StrV StrV;

StrV *strv_create( int strSize );
void strv_resize( StrV *strv, int newSize );
void strv_push_back( StrV *strv, const char *str );
void strv_add_lines( StrV *strv, const StrV *lines );
const char *strv_get( const StrV *strv, int pos );
void strv_set( StrV *strv, int pos, const char *str );
char **strv_ptr( StrV *strv );
int strv_size( const StrV *strv );
void strv_clear( StrV *strv );
int strv_find( const StrV *strv, const char *str );
int strv_find_substr( const StrV *strv, const char *str );
void strv_read_from( StrV *strv, const char *fileName, const char ignoreEmptyLines );
void strv_write_to( StrV *strv, const char *fileName );
void strv_free( StrV **strv );

/* put here all types of vectors you will need */
DECLARE_DICT_TYPE(int)
DECLARE_VECTOR_TYPE(int)
DECLARE_VECTOR_TYPE(char)
DECLARE_VECTOR_TYPE(float)
DECLARE_VECTOR_TYPE(double)
DECLARE_VECTOR_TYPE(IntPair)

/* shuffles a vector of integers, using a vector of intpairs as auxiliary structure */
void vec_int_shuffle( Vec_int *vint, Vec_IntPair *vpair );

/* sorts a vector of int pair by its second field */
void vec_IntPair_sort( Vec_IntPair *vip );

/* two dimensional vector of integers */
typedef Vec_int* Vec_intp;
typedef Vec_IntPair *Vec_IntPairp;
DECLARE_VECTOR_TYPE(Vec_IntPairp)
DECLARE_VECTOR_TYPE(Vec_intp)
typedef Vec_Vec_intp V2D_int;

V2D_int *v2d_create( const int rows );
int v2d_size( const V2D_int *v2d );
void v2d_resize( V2D_int *v2d, const int rows );
void v2d_int_row_push_back( V2D_int *v2d, int row, int value );
void v2d_int_row_insert_unique( V2D_int *v2d, int row, int value );
void v2d_int_row_clear( V2D_int *v2d, int row );
int v2d_int_row_size( const V2D_int *v2d, int row );
int v2d_int_row_get( const V2D_int *v2d, int row, int col );
int *v2d_int_row_ptr( V2D_int *v2d, int row );
char v2d_int_equals( const V2D_int *v2d1, const V2D_int *v2d2 );
void v2d_cpy(  V2D_int *target, const V2D_int *source );
void v2d_free( V2D_int **v2di );

DECLARE_DICT_TYPE( Vec_intp )

#define DEF_INI_CAP 64

#ifdef DEBUG
#define CHECK_NULL_VEC( v ) \
   if ( v==NULL ) { \
      fprintf( stderr, "ERROR: NULL passed as vector in %s:%d\n", __FILE__, __LINE__ ); \
      abort(); \
   }
#else
#define CHECK_NULL_VEC( v ) ;
#endif

#define IMPLEMENT_VECTOR_TYPE(type) \
Vec_##type *vec_##type##_create( ) { \
   return vec_##type##_create_cap( DEF_INI_CAP ); \
} \
Vec_##type *vec_##type##_create_cap( const int iniCap ) { \
   Vec_##type *vec = (Vec_##type *) malloc( sizeof(Vec_##type) ); \
   if (!vec) { \
      fprintf( stderr, "ERROR: out of memory.\n" ); \
      fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
      abort(); \
   } \
   vec->capacity = iniCap; \
   vec->size = 0; \
   vec->array = (type*) malloc( sizeof(type)*vec->capacity ); \
   if (!vec->array) { \
      fprintf( stderr, "ERROR: out of memory.\n" ); \
      fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
      abort(); \
   } \
   return vec; \
} \
Vec_##type *vec_##type##_create_fill( const int size, const type value ) { \
   Vec_##type *vec = vec_##type##_create_cap( size ); \
   int i; \
   for ( i=0 ; (i<size) ; ++i ) \
      vec->array[i] = value; \
   vec->size = size; \
   return vec; \
} \
void vec_##type##_clear( Vec_##type *vec ) { \
   CHECK_NULL_VEC( vec ); \
   vec->size = 0; \
} \
void vec_##type##_free(Vec_##type **vec) { \
   free( (*vec)->array ); \
   free( *vec ); \
   *vec = NULL; \
} \
void vec_##type##_push_back( Vec_##type *vec, const type value ) { \
   if ( vec->size == vec->capacity ) { \
      vec->capacity *= 2; \
      type *p =(type *) realloc( vec->array, sizeof(type)*vec->capacity ); \
      if (!p) { \
         fprintf( stderr, "ERROR: out of memory.\n" ); \
         fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
         abort(); \
      } \
      vec->array = p; \
   } \
   vec->array[vec->size] = value; \
   vec->size++; \
} \
void vec_##type##_push_back_v( Vec_##type *vec, int n, const type v[] ) { \
   if ( vec->size+n > vec->capacity ) { \
      vec->capacity = MAX( vec->capacity*2, vec->size+n+8 ); \
      type *p =(type *) realloc( vec->array, sizeof(type)*vec->capacity ); \
      if (!p) { \
         fprintf( stderr, "ERROR: out of memory.\n" ); \
         fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
         abort(); \
      } \
      vec->array = p; \
   } \
   memcpy( vec->array+vec->size, v, sizeof(type)*n ); \
   vec->size += n; \
} \
type vec_##type##_pop_back( Vec_##type *vec ) { \
   if (vec->size==0) { \
      fprintf( stderr, "ERROR: vec_type_pop_back - trying to remove element from empty vector.\n" ); \
      abort(); \
   } \
   vec->size--; \
   return vec->array[vec->size]; \
} \
void vec_##type##_set( Vec_##type *vec, const int pos, const type value) { \
   if ( (pos>=vec->size) || (pos<0) ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_set.\n"); \
      fprintf(stderr, "\tpos: %d array size: %d array capacity: %d\n", pos, vec->size, vec->capacity ); \
      abort(); \
   } \
   vec->array[pos] = value; \
} \
void vec_##type##_swap( Vec_##type *vec, const int pos1, const int pos2 ) { \
   if ( (pos1>=vec->size) || (pos1<0) ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_swap.\n"); \
      fprintf(stderr, "\tpos1: %d array size: %d array capacity: %d\n", pos1, vec->size, vec->capacity ); \
      abort(); \
   } \
   if ( (pos2>=vec->size) || (pos2<0) ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_swap.\n"); \
      fprintf(stderr, "\tpos2: %d array size: %d array capacity: %d\n", pos2, vec->size, vec->capacity ); \
      abort(); \
   } \
   type aux = vec->array[pos1]; \
   vec->array[pos1] = vec->array[pos2]; \
   vec->array[pos2] = aux; \
}; \
/* resizes without filling new elements */ \
void vec_##type##_resize_nf( Vec_##type *vec, const int size ) { \
    if ( vec->capacity>=size ) { \
        vec->size = size; \
        return; \
    } \
    const int oldSize = vec->size; \
    vec->size = vec->capacity = size; \
    type *p = (type*) realloc( vec->array, sizeof(type)*vec->capacity ); \
    if (!p) { \
        fprintf( stderr, "ERROR: out of memory resizing vector from %d to %d elements.\n", oldSize, size ); \
        fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
        abort(); exit( EXIT_FAILURE ); \
    } \
    vec->array = p; \
}\
void vec_##type##_fset( Vec_##type *vec, const int pos, const type value ) { \
   if ( pos<0 ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_set.\n"); \
      fprintf(stderr, "\tpos: %d array size: %d array capacity: %d\n", pos, vec->size, vec->capacity ); \
      abort(); \
   } \
   if ( pos >= vec->capacity ) { \
      vec->capacity = MAX( 2*vec->capacity, pos+1 ); \
      type *p =(type *) realloc( vec->array, sizeof(type)*vec->capacity ); \
      if (!p) { \
         fprintf( stderr, "ERROR: out of memory.\n" ); \
         fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
         abort(); \
      } \
      vec->array = p; \
   } \
   vec->array[pos] = value; \
   vec->size = MAX( vec->size, pos+1 ); \
} \
type vec_##type##_first( const Vec_##type *vec ) { \
    assert( vec->size ); return vec->array[0]; \
} \
type vec_##type##_last( const Vec_##type *vec ) { \
    assert( vec->size ); return vec->array[vec->size-1]; \
} \
type vec_##type##_get( const Vec_##type *vec, const int pos ) { \
   if ( (pos>=vec->size) || (pos<0) ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_get.\n"); \
      fprintf(stderr, "\tpos: %d array size: %d array capacity: %d\n", pos, vec->size, vec->capacity ); \
      abort(); \
   } \
   return vec->array[pos]; \
} \
type *vec_##type##_getp( Vec_##type *vec, const int pos ) { \
   if ( (pos>=vec->size) || (pos<0) ) { \
      fprintf(stderr, "ERROR: accessing array out of bounds at Vec_##type##_getp.\n"); \
      fprintf(stderr, "\tpos: %d array size: %d array capacity: %d\n", pos, vec->size, vec->capacity ); \
      abort(); \
   } \
   return &(vec->array[pos]); \
} \
void vec_##type##_cpy( Vec_##type *vecTarget, const Vec_##type *vecSource ) { \
   if ( vecTarget->capacity < vecSource->size ) { \
      vecTarget->capacity = vecSource->size; \
      type *p = (type *)realloc( vecTarget->array, sizeof(type)*vecTarget->capacity ); \
      if (!p) { \
         fprintf( stderr, "ERROR: out of memory.\n" ); \
         fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
         abort(); \
      } \
      vecTarget->array = p; \
   } \
   vecTarget->size = vecSource->size; \
   memcpy( vecTarget->array, vecSource->array, vecSource->size*sizeof(type) ); \
} \
int vec_##type##_size( const Vec_##type *vec ) { \
   return vec->size; \
} \
char vec_##type##_equals( const Vec_##type *vec1, const Vec_##type *vec2 ) { \
    if ( vec_##type##_size(vec1) != vec_##type##_size(vec2) ) return 0; \
    int i; const int sz = vec_##type##_size(vec1); \
    for ( i=0 ; (i<sz) ; ++i ) \
        if ( !EQUAL_##type( vec_##type##_get(vec1,i), vec_##type##_get(vec2,i) ) ) \
            return 0; \
    return 1; \
}\
int vec_##type##_capacity( const Vec_##type *vec ) { \
   return vec->capacity; \
} \
type *vec_##type##_ptr( const Vec_##type *vec ) { \
   return vec->array; \
}\
void vec_##type##_resize( Vec_##type *vec, const int size, const type valueNewElements  ) { \
    if ( vec->capacity>=size ) { \
        int i; \
        for ( i=vec->size ; (i<size) ; ++i ) \
            vec->array[i] = valueNewElements; \
        vec->size = size; \
        return; \
    } \
    const int oldSize = vec->size; \
    vec->size = vec->capacity = size; \
    type *p = (type *) realloc( vec->array, sizeof(type)*vec->capacity ); \
    if (!p) { \
        fprintf( stderr, "ERROR: out of memory resizing vector from %d to %d elements.\n", oldSize, size ); \
        fprintf( stderr, "\t%s:%d\n", __FILE__, __LINE__ ); \
        abort(); exit( EXIT_FAILURE ); \
    } \
    vec->array = p; \
    int i; \
    for ( i=oldSize ; (i<vec->size) ; i++ ) \
        vec_##type##_set( vec, i, valueNewElements ); \
}\
int vec_##type##_find( const Vec_##type *vec, const type value ) { \
    int i; \
    for ( i=0 ; (i<vec_##type##_size(vec)) ; ++i ) \
        if (EQUAL_##type( vec->array[i], value )) return i; \
    return NOT_FOUND; \
} \
void vec_##type##_remove(Vec_##type *vec, const type value ) { \
    int pos = vec_##type##_find( vec, value ); \
    if ( pos == NOT_FOUND ) return; \
    if ( pos == vec_##type##_size(vec)-1 ) { vec->size--; return; } \
    vec_##type##_swap( vec, pos, vec_##type##_size(vec)-1 ); \
    vec->size--; \
} \
void vec_##type##_remove_positions( Vec_##type *vec, int nPos, const int positions[] ) { \
    assert( nPos <= vec->size && nPos >= 1 ); \
    assert( vec!=NULL && positions!=NULL ); \
    while ( (--nPos) >=0 ) { \
        int p = positions[nPos]; \
        assert( p>=0 && p<vec->size ); \
        if ( nPos && p<=positions[nPos-1] ) { \
            fprintf( stderr, "Positions should be informed in ascending order.\n" ); \
            abort(); \
        } \
        if ( p == vec->size-1 )  \
            vec->size--; \
        else { \
            vec_##type##_swap( vec, p, vec->size-1 ); \
            vec->size--; \
        } \
    } \
} \
void vec_##type##_insert_unique( Vec_##type *vec, const type value ) { \
    if (vec_##type##_find( vec, value )==NOT_FOUND) \
        vec_##type##_push_back( vec, value ); \
}

#define IMPLEMENT_DICT_TYPE(type) \
Dict_##type *dict_##type##_create( unsigned int hashSize, type defaultValue ) \
{ \
    Dict_##type *res; \
    ALLOCATE( res, Dict_##type ); \
    ALLOCATE_VECTOR_INI( res->cont, Dict_Bucket_##type *, hashSize ); \
    ALLOCATE_VECTOR_INI( res->rowSize, int, hashSize ); \
    ALLOCATE_VECTOR_INI( res->rowCap, int, hashSize ); \
    res->defValue = defaultValue; \
    res->keys = strv_create( STR_SIZE ); \
    res->hashSize = hashSize; \
    return res; \
} \
void dict_##type##_set( Dict_##type *dict, const char *key, const type value ) \
{ \
    unsigned int hashPos = str_hash( key, dict->hashSize ); \
    int i; \
    char found = False; \
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) { \
        if ( strcmp( key, dict->cont[hashPos][i].str) == 0 ) { \
            found = True; \
            break; \
        } \
    } \
    if (found) { \
        dict->cont[hashPos][i].value = value; \
        return; \
    } \
    if ( dict->rowSize[hashPos]+1 > dict->rowCap[hashPos] ) \
    { \
        if (dict->rowSize[hashPos]==0) { \
            ALLOCATE_VECTOR( dict->cont[hashPos], Dict_Bucket_##type, DEF_INI_CAP ); \
            dict->rowCap[hashPos] = DEF_INI_CAP; \
        } else { \
           dict->rowCap[hashPos]*=2; \
           Dict_Bucket_##type *bigger =(Dict_Bucket_##type *) realloc( dict->cont[hashPos], sizeof(Dict_Bucket_##type)*dict->rowCap[hashPos] ); \
           if (!bigger) { \
               fprintf( stderr, "ERROR: no more memory available" ); abort(); exit(EXIT_FAILURE); \
           } \
           dict->cont[hashPos] = bigger; \
        } \
    } \
    strncpy( dict->cont[hashPos][dict->rowSize[hashPos]].str, key, STR_SIZE ); \
    dict->cont[hashPos][dict->rowSize[hashPos]].value = value; \
    dict->cont[hashPos][dict->rowSize[hashPos]].keyPos = strv_size( dict->keys ); \
    dict->rowSize[hashPos]++; \
    strv_push_back( dict->keys, key ); \
}\
void dict_##type##_remove( Dict_##type *dict, const char *key ) { \
    /* finding bucket for key */ \
    unsigned int hashPos = str_hash( key, dict->hashSize ); \
    int i, pBucket = -1; \
    Dict_Bucket_##type *bucket = NULL; \
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) { \
        if ( strcmp( key, dict->cont[hashPos][i].str) == 0 ) { \
            pBucket = i; \
            bucket = &dict->cont[hashPos][i]; \
            break; \
        } \
    } \
    /* if not found getting out */ \
    if ( pBucket == -1 ) return; \
    /* removing key from keys vector */ \
    if ( bucket->keyPos != (strv_size(dict->keys)-1)  ) { \
        /* moving last key in str vector to pos */ \
        const char *keyToUpdate = strv_get(dict->keys, strv_size(dict->keys)-1); \
        strv_set( dict->keys, bucket->keyPos, keyToUpdate ); \
        { \
            unsigned int hashPosK = str_hash( keyToUpdate, dict->hashSize ); \
            Dict_Bucket_##type *bucketK = NULL; \
            int ik; \
            for ( ik=0 ; (ik<dict->rowSize[hashPosK]) ; ik++ ) { \
                if ( strcmp(  keyToUpdate, dict->cont[hashPosK][ik].str) == 0 ) { \
                    bucketK = dict->cont[hashPosK] + ik; \
                    break; \
                } \
            } \
            assert( bucketK!=NULL ); \
            bucketK->keyPos = bucket->keyPos; \
        } \
    } \
    strv_resize( dict->keys, strv_size(dict->keys)-1 ); \
    /* now updating row in hash matrix */ \
    if ( pBucket != dict->rowSize[hashPos]-1 ) { \
        dict->cont[hashPos][pBucket] = dict->cont[hashPos][dict->rowSize[hashPos]-1]; \
    } \
    dict->rowSize[hashPos]--; \
} \
void dict_##type##_iset( Dict_##type *dict, const int key, const type value ) { \
    char str[32]; sprintf( str, "%d", key ); \
    dict_##type##_set( dict, str, value ); \
} \
type dict_##type##_iget( const Dict_##type *dict, const int key ) { \
    char str[32]; sprintf( str, "%d", key ); \
    return dict_##type##_get( dict, str); \
} \
type dict_##type##_get( const Dict_##type *dict, const char *str ) \
{ \
    unsigned int hashPos = str_hash( str, dict->hashSize ); \
    int i; \
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) { \
        if ( strcmp( str, dict->cont[hashPos][i].str) == 0 ) { \
            return dict->cont[hashPos][i].value; \
            break; \
        } \
    } \
    return dict->defValue; \
} \
void dict_##type##_iremove( Dict_##type *dict, int key ) { \
    char strKey[STR_SIZE]; sprintf( strKey, "%d", key ); dict_##type##_remove(dict, strKey); \
} \
int  dict_##type##_size( const Dict_##type *dict ) { \
    return strv_size( dict->keys ); \
} \
const char *dict_##type##_key( const Dict_##type *dict, int pos ) { \
    return strv_get( dict->keys, pos ); \
} \
int dict_##type##_ikey( const Dict_##type *dict, int pos ) { \
    return atoi(strv_get(dict->keys, pos)); \
}\
void dict_##type##_clear( Dict_##type *dict ) { \
    memset( dict->rowSize, 0, dict->hashSize*sizeof(int) );  \
    strv_resize( dict->keys, 0 ); \
} \
void dict_##type##_cpy( Dict_##type *target, const Dict_##type *source ) { \
    dict_##type##_clear( target ); \
    int i; \
    for ( i=0 ; (i<dict_##type##_size(source)) ; ++i ) \
        dict_##type##_set( target, dict_##type##_key(source,i), dict_##type##_get(source, dict_##type##_key(source,i))  ); \
} \
void dict_##type##_free( Dict_##type **dict ) \
{ \
    Dict_##type *d = *dict; \
    unsigned int i; \
    for (i=0;(i<d->hashSize);++i) \
        if (d->rowCap[i]) \
            free( d->cont[i] ); \
    free( (*dict)->cont ); \
    free( (*dict)->rowSize ); \
    free( (*dict)->rowCap ); \
    strv_free( &(*dict)->keys ); \
    free( *dict ); \
}

void strv_sort( StrV *strv );

#endif

