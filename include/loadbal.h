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
typedef unsigned short	ushort;

typedef ushort		bool;

#define true		1
#define false		0

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

#define swap(a, b)				\
	do { typeof(a) __tmp = (a); (a) = (b); (b) = __tmp; } while (0)

#endif // LOADBAL_H
