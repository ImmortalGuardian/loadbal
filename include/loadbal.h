/*
 * Here is some common interface for use in particalur implementation.
 */

#ifndef LOADBAL_H
#define LOADBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PRERROR(s, err)                                                 \
	do {                                                            \
		fprintf(stderr, "%s%s\n", s, strerror(err));		\
		exit(EXIT_FAILURE);                                     \
	} while(0);

#endif
