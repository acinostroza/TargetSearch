#include <R.h>
#include <Rdefines.h>
#include "getLine.h"

int getLine(char **line, int *n, FILE *fp) {

	int c, linepos = 0;

	if(*line == NULL) {
		*n = BUFSIZE;
		*line = (char *) R_chk_calloc(*n, sizeof(char));
	}

	while(1) {
		c = getc(fp);
		linepos++;
		if(linepos == *n) {
			*n += BUFSIZE;
			*line = (char *) R_chk_realloc(*line, *n * sizeof(char));
		}

		if(c == EOF) {
			(*line)[linepos-1] = '\0';
			return -1;
		}

		(*line)[linepos-1] = (char) c;

		if(c == '\n') {
			break;
		}
	}
	(*line)[linepos] = '\0';

	return linepos;
}
