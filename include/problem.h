#ifndef PROBLEM_H
#define PROBLEM_H

#include <loadbal.h>
#include <errno.h>

#define scale(num)	(16 - log2base(num) + 1)

/*
 *        Calculation area
 *        ****************
 *     Y .
 *      /|\
 * Ytop  |_______________
 *       |               |
 *       |               |
 *       |               |
 *       |               |
 *       |               |    X
 * Ylow _|_______________|____\
 *       |                    /
 *        Xlft           Xryt
 */
double Xlft;
double Xryt;
double Ylow;
double Ytop;

double Xlen;
double Ylen;

unsigned int xgridsize;
unsigned int ygridsize;

/*
 * Lottka--Volterra model parameters
 */
double alpha;
double beta;
double gama;
double delta;
double Dn;
double Dm;

/*
 * Calculation method parameters
 */
double h;
double dt;

/*
 * Elementary job unit.
 * It is to be distributed/transfered among processes.
 */
typedef struct {
	double **data[2];
	int num;
	int top;
	int low;
	int lft;
	int ryt;
} job_t;

/*
 * log2base - calculates number of bits up to the most significant non-zero bit
 * (from right to left)
 */
static inline unsigned int log2base(unsigned int n)
{
	if (!n)
		PRERROR("log2base: really want calculate log(0)? Think twice!",
			ERANGE);
	unsigned int base = 1;
	unsigned int index = 0;

	do {
		base <<= 1;
		index++;
	} while (base < n);

	return index;
}

void assist_init(void);
job_t *form_jobs(unsigned int np, unsigned int *jobsnum);

#endif // PROBLEM_H
