/**
 * string utilities for parsing TXT files
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
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


/**
 * returns the number of utf8 bytes or zero if it's not a UTF-8 sequence.
 *
 * The function checks whether a sequence a chars follows the specification of
 * UTF-8 encoding (ref: wikipedia) and return the number of bytes it uses.
 * It does not check whether the character is actually valid, only that the
 * masks are correct.
 */
static inline int isutf8(char *s, int n)
{
#define MASK(x) (((x) & 0xc0) == 0x80)
#define BYTE(k) ((k) < n ? s[(k)] : 0x00)
	char b0 = BYTE(0), b1 = BYTE(1), b2 = BYTE(2), b3 = BYTE(3);
	/* 0xxxxxxx */
	if((b0 & 0x80) == 0x00)
		return 1;
	/* 110xxxxx 10xxxxxx */
	if((b0 & 0xe0) == 0xc0 && MASK(b1))
		return 2;
	/* 1110xxxx 10xxxxxx 10xxxxxx */
	if((b0 & 0xf0) == 0xe0 && MASK(b1) && MASK(b2))
		return 3;
	/* 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx */
	if((b0 & 0xf8) == 0xf0 && MASK(b1) && MASK(b2) && MASK(b3))
		return 4;
	return 0;
#undef MASK
#undef BYTE
}

/* checks for printable chars */
static inline int is_ascii(char c) {
	return (c >= ' ' && c <= '~') || c == '\t' || c == '\n' || c == '\r';
}

/**
 * guess the file type of a file
 *
 * open the file and read the first READSIZE bytes and check whether is
 * a binary file or a text (UTF-8) file. The check is quite simple so
 * there could be false positives (it will guess `binary` by default)
 *
 * @param fname string. path to the file name
 * @return and integer. zero if binary or empty; one if text (utf-8),
 *        negative (minus one) if file cannot be open;
 */

#define READSIZE 256
int file_type(const char * fname)
{
	FILE * fp = fopen(fname, "rb");
	if(fp == NULL)
		return -1;

	char buf[READSIZE];
	int len = fread(buf, 1, sizeof(buf), fp);
	fclose(fp);

	if(len == 0)
		return 0;

	for(int i = 0, j = 0; i < len; i += j) {
		if((j = isutf8(buf + i, len - i)) == 0)
			return 0;
		if(j == 1 && !is_ascii(buf[i]))
			return 0;
	}
	return 1;
}

/**
 * simple routine to get the endianness at runtime
 *
 * @returns 0 for little endian, 1 for big endian
 */
int endianness(void) {
	uint16_t x = 0x1122;
	uint8_t * y = (uint8_t *) &x;
	return *y == 0x11;
}
