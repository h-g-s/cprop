/*
 * memory.h
 * Developed by Haroldo Gambini Santos
 * haroldo@iceb.ufop.br
 */

#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>
#include <sys/types.h>

void *xmalloc( const size_t size );

void *xcalloc( const size_t elements, const size_t size );

void *xrealloc( void *ptr, const size_t size );

#endif /* ifndef MEMORY_H */

