#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include "containers.h"
#include "memory.h"
#include "strutils.h"
#include "macros.h"

static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157, 821, 257, 263, 389, 397, \
                                        457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 9, 53, 59  };

static const unsigned int nHashvalues = sizeof(hashval)/sizeof(int);

int cmp_string( const void *p1, const void *p2 );

unsigned int str_hash( const char *str, const unsigned int hashSize )
{
    const unsigned int len = (unsigned int) MIN( nHashvalues , strlen(str) );

    unsigned int sum = 0, i;

    for ( i=0; (i<len) ; i++ )
        sum += (unsigned int) str[i]*hashval[i];

    return sum%hashSize;
}

/* including code for primitive types */
#define EQUAL_int( x,y ) ( x == y )
IMPLEMENT_VECTOR_TYPE(int)
#define EQUAL_char( x,y ) (x==y)
IMPLEMENT_VECTOR_TYPE(char)
#define EQUAL_double( x,y ) (x==y)
IMPLEMENT_VECTOR_TYPE(double)
#define EQUAL_float( x,y ) (x==y)
IMPLEMENT_VECTOR_TYPE(float)
#define EQUAL_IntPair( x,y ) ( (x.a==y.a) && (x.b==y.b) )
IMPLEMENT_VECTOR_TYPE( IntPair )
#define EQUAL_Vec_intp( x,y ) ( (x==y) )
IMPLEMENT_VECTOR_TYPE( Vec_intp )

#define EQUAL_Vec_IntPairp( x,y ) ( (x==y) )
IMPLEMENT_VECTOR_TYPE( Vec_IntPairp )

IMPLEMENT_DICT_TYPE( int )
IMPLEMENT_DICT_TYPE( Vec_intp )

struct _StrV
{
   int capacity;
   int size;

   char **sv;
   char *s;

   int strSize;
};

#ifdef DEBUG
void strv_check_valid_strv( const StrV *strv );
#endif

void strv_increase_capacity_to( StrV *strv, int newCapacity );

StrV *strv_create( int strSize )
{
#ifdef DEBUG
   assert( strSize > 0 );
#endif
   StrV *res = xmalloc( sizeof(StrV) );

   res->capacity = DEF_INI_CAP;
   res->size = 0;
   res->strSize = strSize;
   res->sv = xmalloc( sizeof(char*)*res->capacity );

   res->s = xcalloc( res->capacity*strSize, sizeof(char) );

   res->sv[0] = res->s;
   int i;
   for ( i=1 ; (i<res->capacity) ; i++ )
      res->sv[i] = res->sv[i-1]+strSize;

#ifdef DEBUG
   strv_check_valid_strv( res );
#endif

   return res;
}

int strv_find( const StrV *strv, const char *str )
{
    int i;
    for ( i=0; (i<strv_size(strv)) ; i++ )
        if ( strcmp( str, strv_get(strv,i) ) == 0 )
            return i;

    return NOT_FOUND;
}

int strv_find_substr( const StrV *strv, const char *str )
{
    int i;
    for ( i=0; (i<strv_size(strv)) ; i++ )
        if ( strstr( strv_get(strv,i), str ) == 0 )
            return i;

    return NOT_FOUND;
}

void strv_write_to( StrV *strv, const char *fileName )
{
    FILE *f = fopen( fileName, "w" );
    int i;
    for ( i=0; (i<strv_size(strv)) ; i++ )
        fprintf( f, "%s\n", strv_get( strv,i) );
    fclose( f );
}

void strv_read_from( StrV *strv, const char *fileName, const char ignoreEmptyLines )
{
    char line[LINE_SIZE], *s;

    FILE *f = fopen( fileName, "r" );
    if (!f)
    {
        fprintf( stderr, "could not open %s.\n", fileName );
        abort(); exit( EXIT_FAILURE );
    }

    while ( (s=fgets(line,FILE_NAME_SIZE,f)) )
    {
        str_remove_sps_eol( s );
        if ( ignoreEmptyLines && (!strlen(s)) )
            continue;

        strv_push_back( strv, s );
    }

    fclose(f);
}

int strv_size( const StrV *strv )
{
    return strv->size;
}

void strv_add_lines( StrV *strv, const StrV *lines )
{
    int i;
    for ( i=0 ; i<strv_size(lines) ; ++i )
        strv_push_back( strv, strv_get( lines, i ) );
}

void strv_push_back( StrV *strv, const char *str )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
   assert( str );
#endif
   if ( strv->size+1 > strv->capacity )
      strv_increase_capacity_to( strv, strv->capacity*2 );

   strncpy( strv->sv[strv->size], str, strv->strSize );
   strv->size++;
}


const char *strv_get( const StrV *strv, int pos )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
   assert( pos >= 0 );
   assert( pos < strv->size );
#endif

   return strv->sv[pos];
}

void strv_resize( StrV *strv, int newSize )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
   assert( newSize>=0 );
#endif
   if ( newSize > strv->capacity )
      strv_increase_capacity_to( strv, newSize );
   strv->size = newSize;
}

#ifdef DEBUG
void strv_check_valid_strv( const StrV *strv )
{
   assert( strv );
   assert( strv->s );
   assert( strv->sv );
   assert( strv->sv[0] == strv->s );
   assert( strv->capacity > 0 );
   assert( strv->size >= 0 );
   assert( strv->strSize > 0 );
   assert( strv->capacity >= strv->size );
}
#endif

void strv_increase_capacity_to( StrV *strv, int newCapacity )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
   assert( newCapacity > strv->capacity );
#endif

   char **sv = xmalloc( sizeof(char*)*newCapacity );
   char *s = xcalloc( newCapacity*strv->strSize, sizeof(char) );
   sv[0] = s;
   int i;
   for ( i=1; (i<newCapacity) ; i++ )
      sv[i] = sv[i-1] + strv->strSize;

   memcpy( s, strv->s, sizeof(char)*strv->size*strv->strSize );

   free( strv->sv );
   free( strv->s );
   strv->sv = sv;
   strv->s = s;

   strv->capacity = newCapacity;
}

void strv_set( StrV *strv, int pos, const char *str )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
   assert( pos>=0 );
   assert( pos<strv->size );
   assert( str );
#endif
   strncpy( strv->sv[pos], str, strv->strSize );
}

char **strv_ptr( StrV *strv )
{
#ifdef DEBUG
   strv_check_valid_strv( strv );
#endif
   return(strv->sv);
}

void strv_free( StrV **strv )
{
#ifdef DEBUG
   strv_check_valid_strv( *strv );
#endif

   free( (*strv)->s );
   free( (*strv)->sv );
   free( (*strv) );
   (*strv) = NULL;
}

V2D_int *v2d_create( const int rows )
{
    V2D_int *v2di = vec_Vec_intp_create_fill( rows, NULL );
    int i;
    for ( i=0; (i<rows) ; i++ )
        vec_Vec_intp_set( v2di, i, vec_int_create() );

    return v2di;
}

void v2d_resize( V2D_int *v2d, const int rows )
{
    if (rows<vec_Vec_intp_capacity(v2d))
    {
        vec_Vec_intp_resize( v2d, rows, NULL );
        return;
    }

    const int oldSize = vec_Vec_intp_size( v2d );
    vec_Vec_intp_resize( v2d, rows, NULL );
    int i;
    for ( i=oldSize ; (i<vec_Vec_intp_size( v2d )) ; ++i )
        if (vec_Vec_intp_get(v2d,i)==NULL)
            vec_Vec_intp_set(v2d,i,vec_int_create());
}

int v2d_int_row_size( V2D_int *v2d, int row )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    return vec_int_size( vi );
}

int v2d_int_row_get( V2D_int *v2d, int row, int col )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    return vec_int_get( vi, col );
}

int *v2d_int_row_ptr( V2D_int *v2d, int row )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    return vec_int_ptr( vi );
}

void v2d_int_row_clear( V2D_int *v2d, int row, int value )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    vec_int_clear( vi );
}

void v2d_int_row_push_back( V2D_int *v2d, int row, int value )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    vec_int_push_back( vi, value );
}

void v2d_int_row_insert_unique( V2D_int *v2d, int row, int value )
{
    Vec_int *vi = vec_Vec_intp_get( v2d, row );
    vec_int_insert_unique( vi, value );
}

int v2d_size( const V2D_int *v2d )
{
    return vec_Vec_intp_size( v2d );
}

int intcmp_func( const void *p1, const void *p2 )
{
    return  (*((const int *)p1)) - (*((const int *)p2));
}

void v2d_free( V2D_int **v2di )
{
    Vec_Vec_intp *v = *v2di;
    const int n = vec_Vec_intp_size(v);
    int i;
    for ( i=0; (i<n) ; i++ )
    {
        Vec_int *vi = vec_Vec_intp_get( v, i );
        vec_int_free( &vi );
    }

    vec_Vec_intp_free( v2di );
}


struct _ISet
{
    int hashSize;

    Vec_Vec_IntPairp *rows;

    Vec_int *elements;
};

ISet *iset_create( int hashSize )
{
    ISet *iset;
    ALLOCATE( iset, ISet );
    iset->hashSize = hashSize;
    iset->rows = vec_Vec_IntPairp_create_fill( hashSize, NULL );
    iset->elements = vec_int_create();

    return iset;
}

void iset_cpy( ISet *target, const ISet *source )
{
    iset_clear( target );
    int i;
    for ( i=0 ; (i<iset_n_elements(source)) ; ++i )
        iset_add( target, iset_element(source,i));
}

void iset_clear( ISet *set )
{
    vec_int_clear( set->elements );
    int i;
    for ( i=0; (i<vec_Vec_IntPairp_size(set->rows)) ; i++ )
    {
        Vec_IntPair *v = vec_Vec_IntPairp_get(set->rows, i);
        if (!v)
            continue;

        vec_IntPair_clear( v );
    }
}

void iset_add( ISet *set, int value )
{
    int pos = value % set->hashSize;

    /* getting vector of bucket */
    Vec_IntPair *vip = vec_Vec_IntPairp_get( set->rows, pos );
    if (!vip)
    {
        vip = vec_IntPair_create();
        vec_Vec_IntPairp_set( set->rows, pos, vip );
    }

    /* checking if it is already inserted */
    int posFound = NOT_FOUND;
    {
        int i;
        for ( i=0; (i<vec_IntPair_size(vip)) ; i++ )
        {
            if ( vec_IntPair_get( vip, i ).a == value )
            {
                posFound = i;
                break;
            }
        }
    }

    if (posFound != NOT_FOUND )
        return;

    vec_int_push_back( set->elements, value );
    vec_IntPair_push_back( vip, (IntPair) { value, (vec_int_size( set->elements )-1) } );
}

IntPair *iset_get_pair( ISet *set, int value )
{
    int pos = value % set->hashSize;

    /* getting vector of bucket */
    Vec_IntPair *vip = vec_Vec_IntPairp_get( set->rows, pos );
    if (!vip)
        return NULL;

    /* checking if it is already inserted */
    int i;
    for ( i=0; (i<vec_IntPair_size(vip)) ; i++ )
        if ( vec_IntPair_get( vip, i ).a == value )
            return vec_IntPair_getp( vip, i );

    return NULL;
}

void iset_remove( ISet *set, int value )
{
    if ( vec_int_size(set->elements) == 0 )
        return;

    int pos = value % set->hashSize;

    /* getting vector of bucket */
    Vec_IntPair *vip = vec_Vec_IntPairp_get( set->rows, pos );
    if (!vip)
    {
        vip = vec_IntPair_create();
        vec_Vec_IntPairp_set( set->rows, pos, vip );
    }

    /* checking if it is already inserted */
    int posFound = NOT_FOUND;
    {
        int i;
        for ( i=0; (i<vec_IntPair_size(vip)) ; i++ )
        {
            if ( vec_IntPair_get( vip, i ).a == value )
            {
                posFound = i;
                break;
            }
        }
    }

    if (posFound == NOT_FOUND )
        return;

    /* position in the elements vector */
    int posValue = vec_IntPair_get( vip, posFound ).b;

    if ( vec_IntPair_size( vip ) == 1 )
    {
        /* just one element */
        vec_IntPair_clear( vip );
    }
    else
    {
        /* more than one */
        /* not last, swapping */
        if ( posFound != vec_IntPair_size( vip )-1 )
            vec_IntPair_swap( vip, posFound, vec_IntPair_size( vip )-1 );

        vec_IntPair_resize( vip, vec_IntPair_size(vip)-1, (IntPair) { 0, 0 } );
    }

    if ( (iset_n_elements(set)>=1) && (posValue!=vec_int_size( set->elements )-1) )
    {
        int currLast = vec_int_last( set->elements );
        vec_int_swap( set->elements, vec_int_size( set->elements )-1, posValue );

        /* the one which was the last element must have its position updated */
        IntPair *ipairLast = iset_get_pair( set, currLast );
        assert( ipairLast );
        ipairLast->b = posValue;
    }

    vec_int_resize( set->elements, vec_int_size(set->elements)-1, 0 );
}

int iset_n_elements( const ISet *set )
{
    return vec_int_size( set->elements );
}

int iset_element( const ISet *set, int idx )
{
    return vec_int_get( set->elements, idx );
}

char iset_has( const ISet *set, int value )
{
    int pos = value % set->hashSize;

    /* getting vector of bucket */
    Vec_IntPair *vip = vec_Vec_IntPairp_get( set->rows, pos );
    if (!vip)
    {
        vip = vec_IntPair_create();
        vec_Vec_IntPairp_set( set->rows, pos, vip );
    }

    /* checking if it is already inserted */
    int posFound = NOT_FOUND;
    {
        int i;
        for ( i=0; (i<vec_IntPair_size(vip)) ; i++ )
        {
            if ( vec_IntPair_get( vip, i ).a == value )
            {
                posFound = i;
                break;
            }
        }
    }

    return (posFound!=NOT_FOUND);
 }

void iset_remove_set( ISet *set, const ISet *toRemove )
{
    int i;
    for ( i=0 ; (i<iset_n_elements(toRemove)) ; ++i )
        iset_remove( set, iset_element(toRemove,i) );
}

void strv_clear( StrV *strv )
{
    strv_resize( strv, 0 );
}

void iset_free( ISet **_iset )
{
    ISet *iset = *_iset;
    int i;
    for ( i=0; (i<vec_Vec_IntPairp_size(iset->rows)) ; i++ )
    {
        Vec_IntPair *v = vec_Vec_IntPairp_get( iset->rows, i );
        if (v)
            vec_IntPair_free( &v );
    }
    vec_Vec_IntPairp_free( &iset->rows );
    vec_int_free( &iset->elements );
    free( iset );
    *_iset = NULL;
}

void strv_sort( StrV *strv )
{
    qsort( strv->sv, strv->size, sizeof(char *), cmp_string );
}

int cmp_string( const void *p1, const void *p2 )
{
    const char **s1 = (const char **) p1;
    const char **s2 = (const char **) p2;

    return strcasecmp( *s1, *s2 );
}

int cmp_IntPair_b( const void *v1, const void *v2 );

void vec_int_shuffle( Vec_int *vint, Vec_IntPair *vpair )
{
    vec_IntPair_clear( vpair );
    int i;
    for ( i=0 ; (i<vec_int_size(vint)) ; ++i )
    {
        const IntPair ipair = { vec_int_get(vint,i), rand() };
        vec_IntPair_push_back( vpair, ipair );
    }

    vec_IntPair_sort(vpair);
    vec_int_clear(vint);
    for ( i=0 ; (i<vec_IntPair_size(vpair)) ; ++i )
        vec_int_push_back( vint, vec_IntPair_getp(vpair,i)->a );
}

int cmp_IntPair_b( const void *v1, const void *v2 )
{
    const IntPair *ip1 = v1;
    const IntPair *ip2 = v2;

    return ip1->b - ip2->b;
}

void vec_IntPair_sort( Vec_IntPair *vip )
{
    qsort( vip->array, vip->size, sizeof(IntPair), cmp_IntPair_b );
}

