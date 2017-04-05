/*
 * Here is some common interface for use in particalur implementation.
 */

#ifndef LOADBAL_H
#define LOADBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned int	uint;
typedef unsigned long	ulong;

#define PRERROR(s, err)                                                 \
	do {                                                            \
		fprintf(stderr, "%s%s\n", s, strerror(err));		\
		exit(EXIT_FAILURE);                                     \
	} while(0);

#define max(x, y) ({                            \
	typeof(x) _max1 = (x);                  \
	typeof(y) _max2 = (y);                  \
	(void) (&_max1 == &_max2);              \
	_max1 > _max2 ? _max1 : _max2; })

#define min(x, y) ({                            \
	typeof(x) _min1 = (x);                  \
	typeof(y) _min2 = (y);                  \
	(void) (&_min1 == &_min2);              \
	_min1 < _min2 ? _min1 : _min2; })

#endif // LOADBAL_H
