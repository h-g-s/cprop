#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "macros.h"

char char_is_invisible( const char c )
{
    return ( (c==' ') || (c=='\t') || (c=='\r') || (c=='\n') );
}

char* applyInversion(char* str) {
   int head, tail;
   char t;
   for ( head=0,tail=(strlen(str)-1) ; head<tail ; head++,tail-- ) {
      t = str[head];
      str[head] = str[tail];
      str[tail] = t;
   }
   return str;
}

char *getFileName(char *destiny, const char *fileWithPath) {
   int i,pos;
   /* Returning till found a slash */
   for ( i=(strlen(fileWithPath)-1),pos=0 ; (i>=0) ; i--,pos++ ) {
      if (fileWithPath[i]=='/') break;
      destiny[pos] = fileWithPath[i];
   }
   destiny[pos]='\0';
   applyInversion(destiny);
   /* Now, removing the . */
   const int endIdx = strlen(destiny)-1;
   for ( i=endIdx ; (i>=0) ; i-- )
   {
      if (destiny[i]=='.')
         break;
   }

   return destiny;
}

char *str_remove_sps_sol( char *dest, const char *str )
{
   char *startDest = dest;
   const char *send = str + strlen( str );
   const char *s = str;
   while ( (char_is_invisible(*s)) && (s<send) )
      s++;
   while ( (s<send) )
   {
      *dest = *s;
      ++dest;
      ++s;
   }
   *dest = '\0';

   return startDest;
}


char* getParamName(char* target, const char* str) {
   unsigned int i;
   unsigned int size;
   unsigned int pos;
   pos = 0;
   size = strlen(str);
   for ( i=0 ; i<size ; i++ ) {
      if (str[i] != '-') {
         if (str[i] != '=') {
            target[pos++] = str[i];
         } else {
            break;
         }
      }
   }
   target[pos] = '\0';
   return target;
}

char* getParamValue(char* target, const char* str) {
   unsigned int i;
   unsigned int size;
   size = strlen(str);
   unsigned int destSize = 0;
   for ( i=0 ; i<size ; i++) {
      if (str[i]=='=') break;
   }
   i++;
   for ( ; i<size ; i++) {
      target[destSize] = str[i];
      destSize++;
   }
   target[destSize] = 0;
   return target;
}

char *str_remove_sps_eol( char *str )
{
   char *p;

   p = str + strlen(str)-1;
   while ( p>=str )
   {
      if (char_is_invisible(*p))
         *p = '\0';
      else
         break;

      --p;
   }

   return str;
}

void str_all_to_upper( char *str )
{
    int l = strlen( str );
    {
        int i;
        for ( i=0 ; (i<l) ; ++i )
            str[i] = toupper(str[i]);
    }

}

void str_fill_spaces_left( char* dest, const char* str, int n )
{
    int len = strlen( str );
    int toFill = n-len;
    toFill = MAX( 0, toFill );
    if (toFill)
        memset( dest ,' ', sizeof(char)*toFill );
    strcpy( dest+toFill, str );
}

void str_fill_spaces_right(char* dest, const char* str, int n)
{
    int len = strlen( str );
    int toFill = n-len;
    toFill = MAX( 0, toFill );
    strcpy( dest, str );
    if (toFill)
        memset( dest+len ,' ', sizeof(char)*toFill );
    *(dest+len+toFill) = '\0';
}

void str_fill_spaces_both( char *dest, const char *str, int n )
{
    str_fill_spaces_left( dest, str, n/2);
    str_fill_spaces_right( dest, str, n/2);
}

char *str_clear_spaces(char *input) {
    int i,j;
    char *output=input;
    for (i = 0, j = 0; i<strlen(input); i++,j++)
    {
        if (input[i]!=' '&&input[i]!='\n')
            output[j]=input[i];
        else
            j--;
    }
    output[j]='\0';
    return output;
}

char *str_remove_dbl_spaces( char *str )
{
    /* removing spaces at extremes */
    char dest[LINE_SIZE];
    str_remove_sps_sol( dest, str );
    str_remove_sps_eol( dest );

    strcpy( str, dest );
    if (strlen(str)<=3)   /* no double spaces possible  A A */
        return str;

    char *s, *v, *prev;

    /* transforming all intermediate
    invisible chars in spaces */
    for ( s=str+1; (*s!='\0') ; ++s )
        if (char_is_invisible(*s))
            *s = ' ';

    v    = str+2,
    s    = str+2,
    prev = str+1;

SKIP_CHAR:
    /* jump just one if valid is a repeated space or tab */
    if ( *prev==' ' )
        while ( *v==' ' )
            ++v;

    if ( *v == '\0' )
        goto CLOSE_STR;

    prev = v;

    if ( v!=s )
        *s = *v;

    ++v; ++s;

    if ( *v != '\0' )
        goto SKIP_CHAR;
CLOSE_STR:
    *s = '\0';

    return str;
}
