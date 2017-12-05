/*
 * memory.c
 * Developed by Haroldo Gambini Santos
 * hsantos@ic.uff.br
 */

#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

void *xmalloc( const size_t size )
{
   void *result = malloc( size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      exit(1);
   }

   return result;
}

void *xcalloc( const size_t elements, const size_t size )
{
   void *result = calloc( elements, size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to callocte %zu bytes.", size);
      exit(1);
   }

   return result;
}

void *xrealloc( void *ptr, const size_t size )
{
   void *result = realloc( ptr, size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      exit(1);
   }
   return result;
}
