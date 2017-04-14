#include <stdlib.h>
#include <problem.h>
#include <loadbal.h>
#include <assert.h>
#include <math.h>

extern double Xlft;
extern double Xryt;
extern double Ylow;
extern double Ytop;

extern double Xlen;
extern double Ylen;

extern uint xgridsize;
extern uint ygridsize;

extern double alpha;
extern double beta;
extern double gama;
extern double delta;
extern double Dn;
extern double Dm;

extern double h;
extern double dt;

const double bord_val = 15.0;
const double amplitude = 5.0;

void assist_init(void)
{
	Xlft = 0.0;
	Xryt = 1.0;
	Ylow = 0.0;
	Ytop = 1.0;

	alpha = 1.0;
	beta = 1.05;
	gama = 1.03;
	delta = 1.04;
	Dn = 0.005;
	Dm = 0.005;

	h = 0.01;
	dt = 0.001;

	Xlen = Xryt - Xlft;
	Ylen = Ytop - Ylow;

	xgridsize = (uint)(Xlen / h) + 1;
	ygridsize = (uint)(Ylen / h) + 1;
}

job_t *form_jobs(uint np, uint *jobsnum)
{
	uint horjobsnum = min(np, scale(np));
	uint vertjobsnum = max(np, scale(np));
	uint xgridperjob = xgridsize / horjobsnum;
	uint ygridperjob = ygridsize / vertjobsnum;
	uint xgridrem = xgridsize % horjobsnum;
	uint ygridrem = ygridsize % vertjobsnum;
	uint i, j;

	*jobsnum = horjobsnum * vertjobsnum;

	uint ind;
	job_t *jobs = calloc(*jobsnum, sizeof(job_t));
	if (!jobs)
		PRERROR("form_jobs: cannot allocate memory: ", ENOMEM);
	for (i = 0; i < vertjobsnum; i++)
		for (j = 0; j < horjobsnum; j++) {
			/* Here we distribute nodes among jobs evenly, getting
			 * smth like this:
			 *
			 * ################
			 * #___________   #
			 * #000|000|000|  #
			 * #000|000|000|  #
			 * #-----------|  #
			 * #000|000|000|  #
			 * #000|000|000|  #
			 * ################
			 */
			ind = i * horjobsnum + j;
			jobs[ind].num = ind;
			jobs[ind].lft = xgridperjob * j;
			jobs[ind].low = ygridperjob * i;
			jobs[ind].ryt = jobs[ind].lft + xgridperjob - 1;
			jobs[ind].top = jobs[ind].low + ygridperjob - 1;

			/* And here we redistribute remaining nodes to fit the
			 * external borders:
			 *
			 * ################
			 * #0000|0000|000|#
			 * #0000|0000|000|#
			 * #-------------|#
			 * #0000|0000|000|#
			 * #0000|0000|000|#
			 * #0000|0000|000|#
			 * ################
			 */
			jobs[ind].lft += min(j, xgridrem);
			jobs[ind].ryt += min(j, xgridrem);
			jobs[ind].low += min(i, ygridrem);
			jobs[ind].top += min(i, ygridrem);
			if (j < xgridrem)
				jobs[ind].ryt++;
			if (i < ygridrem)
				jobs[ind].top++;

			jobs[ind].xnodes = jobs[ind].ryt - jobs[ind].lft + 1;
			jobs[ind].ynodes = jobs[ind].top - jobs[ind].low + 1;
		}
	/* Just to be sure we did everything right; the user isn't supposed
	 * to see theese asserts.
	 */
	assert(jobs[*jobsnum - 1].ryt == xgridsize - 1);
	assert(jobs[*jobsnum - 1].top == ygridsize - 1);

	return jobs;
}


void alloc_memory(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, j, k;
	int num;

	/* Increase grid nodes ammount by two in both directions to keep
	 * boundary values. Allocate space for 3 time layers, as computational
	 * method dictates it.
	 */
	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		for (j = 0; j < 3; j++) {
			alljobs[num].N[j] = calloc(alljobs[num].ynodes +2,
				sizeof(double*));
			alljobs[num].M[j] = calloc(alljobs[num].ynodes +2,
				sizeof(double*));
			if (!(alljobs[num].N[j]) || !(alljobs[num].M[j]))
				PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);

			for (k = 0; k < alljobs[num].ynodes + 2; k++) {
				alljobs[num].N[j][k] = calloc(alljobs[num].xnodes + 2,
					sizeof(double));
				alljobs[num].M[j][k] = calloc(alljobs[num].xnodes + 2,
					sizeof(double));
				if (!(alljobs[num].N[j][k]) || !(alljobs[num].M[j][k]))
					PRERROR("alloc_memory: cannot allocate memory: ",
						ENOMEM);
			}
		}
	}
}

void free_resources(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, j, k;
	int num;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		for (j = 0; j < 3; j++) {
			for (k = 0; k < alljobs[num].ynodes + 2; k++) {
				free(alljobs[num].N[j][k]);
				free(alljobs[num].M[j][k]);
			}
			free(alljobs[num].N[j]);
			free(alljobs[num].M[j]);
		}
	}

	free(activejobs);
	free(alljobs);
}

/* Set testing periodical initial conditions
 */

/* Here we have sin wave */
double M0(double x, double y)
{
	double xperiod = (M_PI * 2) * 5;
	double yperiod = (M_PI * 2) * 5;

	return bord_val + amplitude * sin((x / Xlen) * xperiod) *
		sin((y / Ylen) / yperiod);
}

/* And here - triangle wave */
double N0(double x, double y)
{
	double xperiod = Xlen / 5;
	double yperiod = Ylen / 5;
	double xbord = Xlft;
	double ybord = Ylow;
	double xlean = amplitude / (xperiod / 4);
	double ylean = amplitude / (yperiod / 4);
	double xmul, ymul;

	while (x > xbord + xperiod)
		xbord += xperiod;
	while (y > ybord + yperiod)
		ybord += yperiod;
	x -= xbord;
	y -= ybord;

	if ((x > xperiod / 4) && (x < 3 * xperiod / 4))
		xmul = -xlean * x + 2 * amplitude + bord_val;
	else if (x < xperiod / 4)
		xmul = xlean * x;
	else
		xmul = xlean * x - 4 * amplitude + bord_val;

	if ((y > yperiod / 4) && (y < 3 * yperiod / 4))
		ymul = -ylean * y + 2 * amplitude + bord_val;
	else if (y < yperiod / 4)
		ymul = ylean * y;
	else
		ymul = ylean * y - 4 * amplitude + bord_val;

	return xmul * ymul;
}

void set_init_cond(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, j, k;
	int num;
	double x, y;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];

		/* Here the boundaries (from neighboring cells) are being filled
		 * in. Even those, which are outside of the calculation area,
		 * since they will not distract us anyhow.
		 */
		for (j = 0; j < alljobs[num].ynodes + 2; j++)
			for (k = 0; k < alljobs[num].xnodes + 2; k++) {
				x = Xlft + (alljobs[num].lft + k - 1) * h;
				y = Ylow + (alljobs[num].low + j - 1) * h;

				alljobs[num].N[old][j][k] = N0(x, y);
				alljobs[num].M[old][j][k] = M0(x, y);

				/* Just in case */
				assert((x <= Xryt + h) && (y <= Ytop + h));
			}
	}
}
