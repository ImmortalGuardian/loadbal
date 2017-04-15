#ifndef PROBLEM_H
#define PROBLEM_H

#include <loadbal.h>
#include <errno.h>

#define scale(num)	(16 - log2base(num) + 1)

#define NOCELL	(-1)	/* Indicates that there's no cell in the neighbourhood */
#define NOPROC	(-1)	/* Indicates that there's no proc in the neighbourhood */

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

uint xgridsize;
uint ygridsize;

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
 * Info about neighbours of the cell
 * @nbr_rank: corresponding process' rank
 * @nbr_cell: cell number
 */
typedef struct {
	int nbr_rank;
	int nbr_cell;
} border_t;

typedef struct {
	border_t top, low, lft, ryt;
} bord_info_t;

/*
 * Elementary job unit.
 * It is to be distributed/transfered among processes.
 */
enum layers {old = 0, pred, new};
typedef struct {
	double **N[3];
	double **M[3];
	bord_info_t brds;
	int xnodes;
	int ynodes;
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
static inline uint log2base(uint n)
{
	if (!n)
		PRERROR("log2base: really want calculate log(0)? Think twice!",
			ERANGE);
	uint base = 1;
	uint index = 0;

	do {
		base <<= 1;
		index++;
	} while (base < n);

	return index;
}

void assist_init(void);
job_t *form_jobs(uint np, uint *jobsnum);

double M0(double x, double y);	/* Initial	*/
double N0(double x, double y);	/* conditions	*/

#endif // PROBLEM_H
