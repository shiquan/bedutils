#ifndef NUMBER_HEADER
#define NUMBER_HEADER
#include <stdio.h>
#include <stdlib.h>

extern int get_numbase(const char *s);
extern int is_ieee_magic_val(const char *val);
extern double nondec2num(char *str, int length);
extern int check_num_likely(const char *str);
extern double force2num(char *str);
extern int str2int(char *str);
#endif
