/**
 * string utilities for parsing TXT files
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "strutils.h"

int ascii(const char *line, int n)
{
	for(int i = 0; i < n; i++) {
		int c = (int) line[i];
		if((c >= ' ' && c <= '~') || c == '\t' || c == '\n' || c == '\r')
			continue;
		return 0;
	}
	return 1;
}

int stod(const char *s, double *pd)
{
        char *end = NULL;
        errno = 0;
        double d = strtod(s, &end);
        if(errno != 0 || s == end)
                return 0;
        /* check for white space till end */
        while(*end)
                if(!isspace(*(end++)))
                        return 0;
        *pd = d;
        return 1;
}

char * tokenize(char *str, char sep)
{
	while(*str != '\0') {
		if(*str == sep) {
			*str = '\0';
			return ++str;
		}
		++str;
	}
	return NULL;
}

void untokenize(char *str, int len, char sep)
{
	for(int i = 0; i < len; i++)
		if(str[i] == '\0')
			str[i] = sep;
}

int rstrip(char * s)
{
        int i = strlen(s) - 1;
        while(i >= 0) {
                if(isspace(s[i]))
                        s[i--] = '\0';
                else
                        break;
        }
	return i + 1;
}

int get_col_index(const char *line, const char *column, char sep)
{
	if(line == NULL || column == NULL)
		return -1;

	int n = strlen(column), field = 0;

	while(*line) {
		if(strncmp(line, column, n) == 0 &&
				(line[n] == '\0' || line[n] == sep))
			return field;
		field++;
		while(*line != sep && *line != '\0')
			line++;
		if(*line == sep)
			line++;
	}
	return -1;
}
