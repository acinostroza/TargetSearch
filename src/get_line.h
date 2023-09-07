#ifndef GET_LINE_H
#define GET_LINE_H

#define EMPTY_CHAR -1
#define READSIZE 4096

#include <R.h>
#ifndef USING_R
	#include <stdio.h>
#endif

/**
 * read an entire line from stream
 *
 * This is a simple re-implementation of getline(3). The function
 * supports LF, CRLF, and CR line breaks which are the standard
 * line breaks supported in Linux, Windows and Mac operating systems.
 *
 * The line pointer `lineptr` must be set to NULL, and the initial `size`
 * set to zero on the first call. The function automatically allocates
 * necessary memory to store the line, internally calling R_realloc
 * (this function checks automatically for allocation errors).
 *
 * The pointer `next` is used to store the next character in the stream
 * so it needs to be set to EMPTY_CHAR on the first call.
 *
 * The function returns the number of character read, or zero if
 * it reaches the end of the stream.
 *
 * If there is an allocation error, then returns a negative value,
 * though it should not happen in an R environment
 */
int get_line(char **lineptr, int *size, int *next, FILE *fp);

#endif
