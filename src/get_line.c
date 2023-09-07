
#include "get_line.h"
#include <Rdefines.h>

#define ALLOC_ERROR -1

static inline int _getc(FILE *fp, int *next)
{
	if(*next == EMPTY_CHAR)
		return fgetc(fp);
	int c = *next;
	*next = EMPTY_CHAR;
	return c;
}

static inline int _pushc(char **line, int *size, int c, int k)
{
	if(k >= *size) {
		*size = (*size) ? (*size) * 2 : READSIZE;
		char * tmp = R_Realloc(*line, *size, char);
		if(tmp == NULL)
			return ALLOC_ERROR;
		*line = tmp;
	}
	(*line)[k] = (char) c;
	return k + 1;
}

int get_line(char **line, int *size, int *next, FILE *fp)
{
	int c, len = 0, ret;

	while(1) {
		c = _getc(fp, next);

		if(c == EOF) {
			if(_pushc(line, size, '\0', len) == ALLOC_ERROR)
				return ALLOC_ERROR;
			return len;
		}

		if((ret = _pushc(line, size, c, len)) == ALLOC_ERROR)
			 return ret;
		len = ret;

		if(c == '\n')
			break;

		if(c == '\r') {
			int c2 = _getc(fp, next);
			if(c2 == '\n') {
				if((ret = _pushc(line, size, c2, len - 1)) == ALLOC_ERROR)
					return ret;
				len = ret;
			} else {
				(*line)[len - 1] = '\n';
				*next = c2;
			}
			break;
		}
	}
	if(_pushc(line, size, '\0', len) == ALLOC_ERROR)
		return ALLOC_ERROR;
	return len;
}
