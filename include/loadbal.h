/*
 * Here is some common interface for use in particalur implementation.
 */

#ifndef LOADBAL_H
#define LOADBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

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

#define WLOAD_TAG	0xFEED
#define JOBPREP_TAG	0xDEAD
#define JOBEXCH_TAG	0xBEEF
#define JOBTRANS_TAG	0xBEAD
#define READY_TAG	0xFAAC

#define NSECPERSEC	1000000000

/* Here are some functions for managing timespec structures
 */
static inline struct timespec tssub(struct timespec lhs, struct timespec rhs)
{
	struct timespec res;

	res.tv_sec = lhs.tv_sec - rhs.tv_sec;
	res.tv_nsec = lhs.tv_nsec - rhs.tv_nsec;
	if (res.tv_nsec < 0) {
		res.tv_nsec = NSECPERSEC + res.tv_nsec;
		res.tv_sec--;
	}

	return res;
}

static inline ulong ts_to_ns(struct timespec ts)
{
	return (ulong)((long)(ts.tv_sec * NSECPERSEC) + ts.tv_nsec);
}

#endif // LOADBAL_H
