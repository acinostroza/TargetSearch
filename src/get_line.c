#include <R.h>
#include <Rdefines.h>
#include "get_line.h"

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
		*line = R_Realloc(*line, *size, char);
	}
	(*line)[k] = (char) c;
	return k + 1;
}

int get_line(char **line, int *size, int *next, FILE *fp)
{
	int c, len = 0;

	while(1) {
		c = _getc(fp, next);

		if(c == EOF) {
			_pushc(line, size, '\0', len);
			return len;
		}

		len = _pushc(line, size, c, len);

		if(c == '\n')
			break;

		if(c == '\r') {
			int c2 = _getc(fp, next);
			if(c2 == '\n') {
				len = _pushc(line, size, c2, len - 1);
			} else {
				(*line)[len - 1] = '\n';
				*next = c2;
			}
			break;
		}
	}
	_pushc(line, size, '\0', len);
	return len;
}
