#pragma once

int rstrip(char * s);
char * tokenize(char *str, char sep);
void untokenize(char *str, int len, char sep);
int stod(const char *s, double *pd);
int ascii(const char *line, int n);
int get_col_index(const char *, const char *, char);
