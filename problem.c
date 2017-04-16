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

			if (i == 0) 
				jobs[ind].brds.low.nbr_cell = NOCELL;
			else
				jobs[ind].brds.low.nbr_cell = ind - horjobsnum;

			if (i == vertjobsnum - 1)
				jobs[ind].brds.top.nbr_cell = NOCELL;
			else
				jobs[ind].brds.top.nbr_cell = ind + horjobsnum;

			if (j == 0)
				jobs[ind].brds.lft.nbr_cell = NOCELL;
			else
				jobs[ind].brds.lft.nbr_cell = ind - 1;

			if (j == horjobsnum - 1)
				jobs[ind].brds.ryt.nbr_cell = NOCELL;
			else
				jobs[ind].brds.ryt.nbr_cell = ind + 1;
		}
	/* Just to be sure we did everything right; the user isn't supposed
	 * to see theese asserts.
	 */
	assert(jobs[*jobsnum - 1].ryt == xgridsize - 1);
	assert(jobs[*jobsnum - 1].top == ygridsize - 1);

	return jobs;
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

/* This set of simple functions is to determine whether the cell has shared
 * borders.
 */
inline bool comm_lft_bord(job_t *job) {
	if (job->brds.lft.nbr_rank != NOPROC && job->brds.lft.nbr_rank != job->rank)
		return true;
	else
		return false;
}

inline bool comm_ryt_bord(job_t *job) {
	if (job->brds.ryt.nbr_rank != NOPROC && job->brds.ryt.nbr_rank != job->rank)
		return true;
	else
		return false;
}

inline bool comm_top_bord(job_t *job) {
	if (job->brds.top.nbr_rank != NOPROC && job->brds.top.nbr_rank != job->rank)
		return true;
	else
		return false;
}

inline bool comm_low_bord(job_t *job) {
	if (job->brds.low.nbr_rank != NOPROC && job->brds.low.nbr_rank != job->rank)
		return true;
	else
		return false;
}

inline bool is_sharing_bord(job_t *job) {
	return comm_lft_bord(job) || comm_ryt_bord(job) ||
		comm_top_bord(job) || comm_low_bord(job);
}


static inline bool border_on_lft(job_t *job)
{
	if (job->lft == 0)
		return true;
	else
		return false;
}

static inline bool border_on_ryt(job_t *job)
{
	if (job->ryt == xgridsize - 1)
		return true;
	else
		return false;
}

static inline bool border_on_low(job_t *job)
{
	if (job->low == 0)
		return true;
	else
		return false;
}

static inline bool border_on_top(job_t *job)
{
	if (job->top == ygridsize - 1)
		return true;
	else
		return false;
}

static inline bool is_bound_cell(job_t *job)
{
	return (border_on_low(job) || border_on_lft(job) || border_on_top(job) ||
			border_on_ryt(job));
}

static void calc_predict(job_t *job)
{
	int j, k;
	int lft_indnt, ryt_indnt, low_indnt, top_indnt;
	double sigma = dt / (2 * h * h);
	double **Nold, **Mold, **Npred, **Mpred, **Nnew, **Mnew;

	lft_indnt = ryt_indnt = low_indnt = top_indnt = 0;

	if (is_bound_cell(job)) {
		if (border_on_lft(job))
			lft_indnt = 1;
		if (border_on_ryt(job))
			ryt_indnt = 1;
		if (border_on_low(job))
			low_indnt = 1;
		if (border_on_top(job))
			top_indnt = 1;
	}

	Nold = job->N[old];
	Mold = job->M[old];
	Npred = job->N[pred];
	Mpred = job->M[pred];
	Nnew = job->N[new];
	Mnew = job->M[new];
	for (j = 1 + low_indnt; j <= job->ynodes - top_indnt; j++)
	for (k = 1 + lft_indnt; k <= job->xnodes - ryt_indnt; k++) {
		Npred[j][k] = Nold[j][k] + (Dn * sigma) * ((Nold[j][k+1] - 2 * Nold[j][k] +
			Nold[j][k-1]) + (Nold[j+1][k] - 2 * Nold[j][k] + Nold[j-1][k]));
		Mpred[j][k] = Mold[j][k] + (Dm * sigma) * ((Mold[j][k+1] - 2 * Mold[j][k] +
			Mold[j][k-1]) + (Mold[j+1][k] - 2 * Mold[j][k] + Mold[j-1][k]));
	}

	if (is_bound_cell(job)) {
		if (border_on_lft(job)) {
			k = 1;
			for (j = 1; j <= job->ynodes; j++) {
				Npred[j][k] = Nnew[j][k] = bord_val;
				Mpred[j][k] = Mnew[j][k] = bord_val;
			}
		}
		if (border_on_ryt(job)) {
			k = job->xnodes;
			for (j = 1; j <= job->ynodes; j++) {
				Npred[j][k] = Nnew[j][k] = bord_val;
				Mpred[j][k] = Mnew[j][k] = bord_val;
			}
		}
		if (border_on_low(job)) {
			j = 1;
			for (k = 1; k <= job->xnodes; k++) {
				Npred[j][k] = Nnew[j][k] = bord_val;
				Mpred[j][k] = Mnew[j][k] = bord_val;
			}
		}
		if (border_on_top(job)) {
			j = job->ynodes;
			for (k = 1; k <= job->xnodes; k++) {
				Npred[j][k] = Nnew[j][k] = bord_val;
				Mpred[j][k] = Mnew[j][k] = bord_val;
			}
		}
	}
}

static inline double fn(double N, double M)
{
	return (alpha - gama * M) * N;
}

static inline double fm(double N, double M)
{
	return (-beta + delta * N) * M;
}

enum {n = 0, m};
static double (*f[2]) (double N, double M) = {fn, fm};

static void runge_kutta(job_t *job, double dt)
{
	double k1[2], k2[2], k3[2], k4[2];
	int j, k;
	double **Npred, **Mpred, **Nnew, **Mnew;
	double N, M;

	Npred = job->N[pred];
	Mpred = job->M[pred];
	Nnew = job->N[new];
	Mnew = job->M[new];

	for (j = 1; j <= job->ynodes; j++)
		for (k = 1; k <= job->xnodes; k++) {
			N = Npred[j][k];
			M = Mpred[j][k];

			k1[n] = f[n](N, M);
			k1[m] = f[m](N, M);

			k2[n] = f[n](N + k1[n] * dt / 2, M + k1[m] * dt / 2);
			k2[m] = f[m](N + k1[n] * dt / 2, M + k1[m] * dt / 2);

			k3[n] = f[n](N + k2[n] * dt / 2, M + k2[m] * dt / 2);
			k3[m] = f[m](N + k2[n] * dt / 2, M + k2[m] * dt / 2);

			k4[n] = f[n](N + k3[n] * dt, M + k3[m] * dt);
			k4[m] = f[m](N + k2[n] * dt, M + k2[m] * dt);

			Nnew[j][k] = Npred[j][k] + (dt / 6) * (k1[n] + 2 * k2[n] +
				2 * k3[n] + k4[n]);
			Mnew[j][k] = Mpred[j][k] + (dt / 6) * (k1[m] + 2 * k2[m] +
				2 * k3[m] + k4[m]);
		}
}

void make_timestep(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i;
	int num;
	job_t *job;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);
		calc_predict(job);
		runge_kutta(job, dt);

		swap(job->N[new], job->N[old]);
		swap(job->M[new], job->M[old]);
	}
}
