#include <assert.h>
#include <assist.h>
#include <problem.h>
#include <limits.h>

/*int which_proc(int cellnum, int np, uint alljobsnum) {
		return cellnum % np;
}
*/

void map_nbrs(int rank, int np, job_t *alljobs, uint alljobsnum, uint *activejobs,
		uint actjobsnum, int *jobsmap)
{
	int i, j, num;
	int *ranks, *loc_ranks, *sizes;
	uint max, *inds, *loc_inds;
	job_t *job;

	sizes = calloc(np, sizeof(uint));
	sizes[rank] = actjobsnum;
	MPI_Allgather(&(sizes[rank]), 1, MPI_UNSIGNED, sizes, 1, MPI_UNSIGNED,
			MPI_COMM_WORLD);
	max = 0;
	for (i = 0; i < np; i++)
		if (sizes[i] > max)
			max = sizes[i];
	free(sizes);

	inds = calloc(np*max, sizeof(uint));
	ranks = calloc(np*max, sizeof(int));
	loc_inds = calloc(max, sizeof(uint));
	loc_ranks = calloc(max, sizeof(int));
	for (i = 0; i < max; i++)
		if (i >= actjobsnum) {
			loc_inds[i] = -1;
			loc_ranks[i] = -1;
		}
		else {
			loc_inds[i] = activejobs[i];
			loc_ranks[i] = rank;
		}
	MPI_Allgather(loc_inds, max, MPI_UNSIGNED, inds, max, MPI_UNSIGNED,
			MPI_COMM_WORLD);
	MPI_Allgather(loc_ranks, max, MPI_INT, ranks, max, MPI_INT, MPI_COMM_WORLD);
	free(loc_ranks);
	free(loc_inds);

	j = 0;
	for (i = 0; i < np*max; i++)
		if (inds[i] != -1) {
			jobsmap[inds[i]] = ranks[i];
			j++;
		}
		else
			assert(ranks[i] == -1);
	assert(j == alljobsnum);

	free(inds);
	free(ranks);

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);

		job->rank = rank;

		if (job->brds.low.nbr_cell == NOCELL)
			job->brds.low.nbr_rank = NOPROC;
		else
			job->brds.low.nbr_rank = jobsmap[job->brds.low.nbr_cell];

		if (job->brds.top.nbr_cell == NOCELL)
			job->brds.top.nbr_rank = NOPROC;
		else
			job->brds.top.nbr_rank = jobsmap[job->brds.top.nbr_cell];

		if (job->brds.lft.nbr_cell == NOCELL)
			job->brds.lft.nbr_rank = NOPROC;
		else
			job->brds.lft.nbr_rank = jobsmap[job->brds.lft.nbr_cell];

		if (job->brds.ryt.nbr_cell == NOCELL)
			job->brds.ryt.nbr_rank = NOPROC;
		else
			job->brds.ryt.nbr_rank = jobsmap[job->brds.ryt.nbr_cell];
	}
}

uint
*distr_jobs(int rank, int np, job_t *alljobs, uint alljobsnum, uint *actjobsnum,
		int *jobsmap)
{
	uint rem, i, currnum, jobsperproc;
	uint *activejobs;

	jobsperproc = alljobsnum / np;
	*actjobsnum = jobsperproc;
	rem = alljobsnum % np;
	if (rank < rem)
		*actjobsnum += 1;
	activejobs = (uint *)calloc(*actjobsnum, sizeof(int));

	currnum = rank;
	for (i = 0; i < *actjobsnum; i++) {
		activejobs[i] = currnum;
		currnum += np;
	}

	map_nbrs(rank, np, alljobs, alljobsnum, activejobs, *actjobsnum, jobsmap);

	return activejobs;
}

static inline uint get_jnum(int i, int j, int xsize)
{
	return (i-1) * xsize + (j-1);
}
uint
*distr_jobs1(int rank, int np, job_t *alljobs, uint alljobsnum, uint *actjobsnum,
		int *jobsmap)
{
	uint others, hlf, perproc, rem;
	uint sum, prev, *jobnums, *activejobs;
	int i, j, k, ysize, xsize, dir;

	jobnums = calloc(np, sizeof(uint));

	hlf = alljobsnum / 2;
	jobnums[0] = hlf;
	others = alljobsnum - hlf;
	perproc = others / (np-1);
	rem = others - perproc*(np-1);

	for (i = 1; i < np; i++)
		jobnums[i] = perproc;
	for (i = 1; i < np; i++) {
		if (rem == 0)
			break;
		jobnums[i]++;
		rem--;
	}

	sum = 0;
	for (i = 0; i < np; i++)
		sum += jobnums[i];
	assert(sum == alljobsnum);

	*actjobsnum = jobnums[rank];
	activejobs = (uint *)calloc(*actjobsnum, sizeof(int));

	prev = 0;
	for (i = 0; i < rank; i++)
		prev += jobnums[i];
	ysize = alljobsnum / np;
	xsize = np;
	assert(ysize*xsize == alljobsnum);

	j = prev / ysize + 1;
	if (j%2 == 1)
		i = prev % ysize + 1;
	else
		i = ysize - prev % ysize;
	/* decide wether we're going up or down */
	if (j%2 == 0)
		dir = 3;
	else
		dir = 1;
	if (i == 1 || i == ysize) {
		if (prev % ysize == 0)
			dir = 2;
		else
			if (i == 1)
				dir = 3;
			else
				dir = 1;
	}
	for (k = 0; k < *actjobsnum; k++) {
		activejobs[k] = get_jnum(i, j, xsize);
		if (dir == 1) {
			i++;
			if (i > ysize) {
				i--;
				dir = 2;
				j++;
			}
		}
		else if (dir == 2) {
			if (i == 1) {
				i++;
				dir = 1;
			}
			else {
				i--;
				dir = 3;
			}
		}
		else if (dir == 3) {
			i--;
			if (i == 0) {
				i++;
				dir = 2;
				j++;
			}
		}
	}
	/*
	for (i = 0; i < *actjobsnum; i++)
		printf("%d, ", activejobs[i]);
	printf("\n");
	fflush(stdout);
	*/

	free(jobnums);
	map_nbrs(rank, np, alljobs, alljobsnum, activejobs, *actjobsnum, jobsmap);

	return activejobs;
}

static void alloc_stor(job_t *job)
{
	int j, k;

	/* Increase grid nodes ammount by two in both directions to keep
	 * boundary values. Allocate space for 3 time layers, as computational
	 * method dictates it.
	 */
	for (j = 0; j < 3; j++) {
		job->N[j] = calloc(job->ynodes + 2, sizeof(double*));
		job->M[j] = calloc(job->ynodes + 2, sizeof(double*));
		if (!(job->N[j]) || !(job->M[j]))
			PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);

		for (k = 0; k < job->ynodes + 2; k++) {
			job->N[j][k] = calloc(job->xnodes + 2, sizeof(double));
			job->M[j][k] = calloc(job->xnodes + 2, sizeof(double));
			if (!(job->N[j][k]) || !(job->M[j][k]))
				PRERROR("alloc_memory: cannot allocate memory: ",
					ENOMEM);
		}
	}
}

static void alloc_bords(job_t *job)
{

	if (comm_lft_bord(job)) {
		job->clms_snd[0] = calloc(job->ynodes * 2, sizeof(double));
		job->clms_rcv[0] = calloc(job->ynodes * 2, sizeof(double));
		if (!(job->clms_snd[0]) || !(job->clms_rcv[0]))
			PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);
	}

	if (comm_ryt_bord(job)) {
		job->clms_snd[1] = calloc(job->ynodes * 2, sizeof(double));
		job->clms_rcv[1] = calloc(job->ynodes * 2, sizeof(double));
		if (!(job->clms_snd[1]) || !(job->clms_rcv[1]))
			PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);
	}

	if (comm_low_bord(job)) {
		job->rows_snd[0] = calloc(job->xnodes * 2, sizeof(double));
		job->rows_rcv[0] = calloc(job->xnodes * 2, sizeof(double));
		if (!(job->rows_snd[0]) || !(job->rows_rcv[0]))
			PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);
	}

	if (comm_top_bord(job)) {
		job->rows_snd[1] = calloc(job->xnodes * 2, sizeof(double));
		job->rows_rcv[1] = calloc(job->xnodes * 2, sizeof(double));
		if (!(job->rows_snd[1]) || !(job->rows_rcv[1]))
			PRERROR("alloc_memory: cannot allocate memory: ", ENOMEM);
	}
}

void alloc_cell(job_t *job)
{
	alloc_stor(job);
	alloc_bords(job);
}

enum {n = 0, m};
void alloc_memory(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, num;
	job_t *job;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);
		alloc_cell(job);
	}
}

/* Count number of borders (in this case, boundary segments) the process shares
 * with other processes. It's used then while exchanging boundary values.
 */
uint count_nbredges(job_t *alljobs, uint *activejobs, uint actjobsnum, int rank) {
	uint cnt = 0;
	int i, num;
	job_t *job;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);

		if (job->brds.top.nbr_rank != NOPROC)
			if (job->brds.top.nbr_rank != rank)
				cnt++;
		if (job->brds.low.nbr_rank != NOPROC)
			if (job->brds.low.nbr_rank != rank)
				cnt++;
		if (job->brds.lft.nbr_rank != NOPROC)
			if (job->brds.lft.nbr_rank != rank)
				cnt++;
		if (job->brds.ryt.nbr_rank != NOPROC)
			if (job->brds.ryt.nbr_rank != rank)
				cnt++;
	}

	return cnt;
}

MPI_Request *prep_shrreqs(uint nbredgenum) {
	return calloc(nbredgenum * 2, sizeof(MPI_Request));
}

MPI_Request *prep_wldreqs(uint nbrsnum)
{
	return calloc(nbrsnum * 2, sizeof(MPI_Request));
}

int *get_nbrs(int rank, int np, job_t *alljobs, uint *activejobs,
		uint actjobsnum, uint *nbrsnum)
{
	uint cnt = 0;
	int i, j, num, *nbrs;
	job_t *job;
	ushort *procmap;

	procmap = calloc(np, sizeof(ushort));
	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);
		
		if (job->brds.top.nbr_rank != NOPROC &&
			job->brds.top.nbr_rank != rank)
			if (procmap[job->brds.top.nbr_rank] == 0) {
				procmap[job->brds.top.nbr_rank] = 1;
				cnt++;
			}
		if (job->brds.low.nbr_rank != NOPROC &&
			job->brds.low.nbr_rank != rank)
			if (procmap[job->brds.low.nbr_rank] == 0) {
				procmap[job->brds.low.nbr_rank] = 1;
				cnt++;
			}
		if (job->brds.lft.nbr_rank != NOPROC &&
			job->brds.lft.nbr_rank != rank)
			if (procmap[job->brds.lft.nbr_rank] == 0) {
				procmap[job->brds.lft.nbr_rank] = 1;
				cnt++;
			}
		if (job->brds.ryt.nbr_rank != NOPROC &&
			job->brds.ryt.nbr_rank != rank)
			if (procmap[job->brds.ryt.nbr_rank] == 0) {
				procmap[job->brds.ryt.nbr_rank] = 1;
				cnt++;
			}
	}

	nbrs = calloc(cnt, sizeof(int));
	*nbrsnum = cnt;
	cnt = 0;
	for (j = 0; j < np; j++)
		if (procmap[j] == 1) {
			nbrs[cnt] = j;
			cnt++;
		}
	assert(cnt == *nbrsnum);

	free(procmap);
	return nbrs;
}

static void norm_wloads(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	job_t *job;
	int num, i;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);

		job->avg_ctime = job->ctime / job->iternum;
		job->ctime = 0;
		job->iternum = 0;
	}
}

static ulong sum_avg_wloads(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	job_t *job;
	int i, num;
	ulong res = 0;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);
		res += job->avg_ctime;
	}

	return res;
}

static ulong *share_wloads(int rank, job_t *alljobs, uint *activejobs,
		uint actjobsnum, int *nbrs, uint nbrsnum, MPI_Request *wloadreqs,
		uint total_wload)
{
	void *buf;
	int cnt, proc, tag;
	int i;
	uint reqnum;
	MPI_Request *req;
	ulong *wloads, tot;

	if (nbrsnum == 0)
		return NULL;

	tot = total_wload;
	wloads = calloc(nbrsnum, sizeof(ulong));
	reqnum = nbrsnum * 2;
	//printf("%d: %lu\n", rank, total_wload);
	//fflush(stdout);

	for (i = 0; i < nbrsnum; i++) {
		buf = (void *)(&(tot));
		cnt = 1;
		proc = nbrs[i];
		tag = WLOAD_TAG;
		req = &(wloadreqs[i*2]);
		MPI_Isend(buf, cnt, MPI_UNSIGNED_LONG, proc, tag,
				MPI_COMM_WORLD, req);

		buf = (void *)(&(wloads[i]));
		req = &(wloadreqs[i*2 + 1]);
		MPI_Irecv(buf, cnt, MPI_UNSIGNED_LONG, proc, tag,
				MPI_COMM_WORLD, req);
	}

	MPI_Waitall(reqnum, wloadreqs, MPI_STATUSES_IGNORE);

	return wloads;
}

static int sort_wl_diffs(int *nbrs, uint nbrsnum, long *wl_diffs,
		long **wldiffs_srt, int **nbrs_srt, int rank)
{
	int i, j, num;
	int *buf_nbrs;
	long *buf_wl;

	num = 0;
	for (i = 0; i < nbrsnum; i++)
		if (wl_diffs[i] > 0)
			num++;
	if (num == 0) {
		*wldiffs_srt = NULL;
		*nbrs_srt = NULL;
		return 0;
	}
	*wldiffs_srt = calloc(num, sizeof(long));
	*nbrs_srt = calloc(num, sizeof(int));

	buf_wl = calloc(nbrsnum, sizeof(long));
	buf_nbrs = calloc(nbrsnum, sizeof(int));
	memcpy(buf_wl, wl_diffs, nbrsnum * sizeof(long));
	memcpy(buf_nbrs, nbrs, nbrsnum * sizeof(int));

	for (i = 0; i < nbrsnum; i++)
		for (j = nbrsnum - 1; j > i; j--)
			if (buf_wl[j] > buf_wl[j-1]) {
				swap(buf_wl[j], buf_wl[j-1]);
				swap(buf_nbrs[j], buf_nbrs[j-1]);
			}
	for (i = 0; i < num; i++) {
		(*wldiffs_srt)[i] = buf_wl[i];
		(*nbrs_srt)[i] = buf_nbrs[i];
	}

	free(buf_wl);
	free(buf_nbrs);

	return num;
}

static int get_srtjobs(job_t *alljobs, uint *activejobs, uint actjobsnum,
		uint **jobs_srt, int nbr_rank)
{
	int i, j, num, cnt;
	ulong ctime1, ctime2;
	job_t *job;

	cnt = 0;
	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);

		if (job->brds.ryt.nbr_rank == nbr_rank ||
			job->brds.lft.nbr_rank == nbr_rank ||
			job->brds.top.nbr_rank == nbr_rank ||
			job->brds.low.nbr_rank == nbr_rank)
			cnt++;
	}
	if (cnt == 0) {
		*jobs_srt = NULL;
		return 0;
	}

	*jobs_srt = calloc(cnt, sizeof(uint));
	j = 0;
	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);

		if (job->brds.ryt.nbr_rank == nbr_rank ||
			job->brds.lft.nbr_rank == nbr_rank ||
			job->brds.top.nbr_rank == nbr_rank ||
			job->brds.low.nbr_rank == nbr_rank) {
			(*jobs_srt)[j] = num;
			j++;
		}
	}

	for (i = 0; i < cnt; i++)
		for (j = cnt-1; j > i; j--) {
			num = (*jobs_srt)[j];
			job = &(alljobs[num]);
			ctime1 = job->avg_ctime;

			num = (*jobs_srt)[j-1];
			job = &(alljobs[num]);
			ctime2 = job->avg_ctime;

			if (ctime1 > ctime2)
				swap((*jobs_srt)[j], (*jobs_srt)[j-1]);
		}

	return cnt;
}

static int choose_to_send(job_t *alljobs, uint *jobs_srt, int jobsrtnum,
		uint nbrsnum, uint **jobs_to_snd, long wl_diff, ushort *jbusymap)
{
	int i, j, num, cnt;
	bool *curmap, *optmap;
	ulong *job_wloads;
	double target, thrhold;
	long optsum, cursum;
	job_t *job;

	target = (double)wl_diff / (nbrsnum + 1);
	thrhold = (1.0 + THRHOLD_PERC / 100.0) * target;

	job_wloads = calloc(jobsrtnum, sizeof(ulong));
	curmap = calloc(jobsrtnum, sizeof(bool));
	optmap = calloc(jobsrtnum, sizeof(bool));
	for (i = 0; i < jobsrtnum; i++) {
		num = jobs_srt[i];
		job = &(alljobs[num]);
		job_wloads[i] = job->avg_ctime;
	}

	optsum = 0;
	for (i = 0; i < jobsrtnum; i++) {
		cursum = 0;
		memset((void *)curmap, '\0', jobsrtnum*sizeof(bool));
		for (j = i; j < jobsrtnum; j++) {
			num = jobs_srt[j];
			if ((jbusymap[num] == 0) &&
				((double)(cursum + job_wloads[j]) < thrhold)) {
				curmap[j] = true;
				cursum += job_wloads[j];
			}
		}
		if (cursum > optsum) {
			optsum = cursum;
			memcpy(optmap, curmap, jobsrtnum*sizeof(bool));
		}
	}

	cnt = 0;
	for (i = 0; i < jobsrtnum; i++)
		if (optmap[i]) {
			num = jobs_srt[i];
			jbusymap[num] = 1;
			cnt++;
		}
	*jobs_to_snd = calloc(cnt, sizeof(uint));
	j = 0;
	for (i = 0; i < jobsrtnum; i++)
		if (optmap[i]) {
			num = jobs_srt[i];
			(*jobs_to_snd)[j] = num;
			j++;
		}
	assert(j == cnt);

	free(job_wloads);
	free(curmap);
	free(optmap);

	return cnt;
}

static void get_to_rcv(int *nbrs, uint nbrsnum, int *nbrs_srt, int nbrsrtnum,
	int *jobsndnum, int *jobrcvnum, uint **jobs_to_snd, uint **jobs_to_rcv,
	int *numtosnd, int *numtorcv, int rank)
{
	int i, j, k, *num_to_snd;
	int *buf, count, dest, src, tag;
	uint *ubuf;
	MPI_Request *sndreqs, *rcvreqs, *req;

	sndreqs = calloc(nbrsnum, sizeof(MPI_Request));
	rcvreqs = calloc(nbrsnum, sizeof(MPI_Request));
	num_to_snd = calloc(nbrsnum, sizeof(int));

	for (i = 0; i < nbrsnum; i++) {
		for (j = 0; j < nbrsrtnum; j++)
			if (nbrs[i] == nbrs_srt[j])
				break;
		if (j == nbrsrtnum)
			num_to_snd[i] = 0;
		else
			num_to_snd[i] = jobsndnum[j];

		buf = &(num_to_snd[i]);
		count = 1;
		dest = nbrs[i];
		tag = JOBPREP_TAG;
		req = &(sndreqs[i]);
		MPI_Isend(buf, count, MPI_INT, dest, tag, MPI_COMM_WORLD, req);

		buf = &(jobrcvnum[i]);
		src = nbrs[i];
		req = &(rcvreqs[i]);
		MPI_Irecv(buf, count, MPI_INT, src, tag, MPI_COMM_WORLD, req);
	}

	MPI_Waitall(nbrsnum, rcvreqs, MPI_STATUSES_IGNORE);
	MPI_Waitall(nbrsnum, sndreqs, MPI_STATUSES_IGNORE);

	for (i = 0; i < nbrsnum; i++)
		jobs_to_rcv[i] = calloc(jobrcvnum[i], sizeof(uint));

	k = 0;
	for (i = 0; i < nbrsrtnum; i++) {
		if (jobsndnum[i] > 0) {
			ubuf = jobs_to_snd[i];
			count = jobsndnum[i];
			dest = nbrs_srt[i];
			tag = JOBEXCH_TAG;
			req = &(sndreqs[k++]);
			MPI_Isend(ubuf, count, MPI_UNSIGNED, dest, tag,
				MPI_COMM_WORLD, req);
		}
	}
	j = 0;
	for (i = 0; i < nbrsnum; i++) {
		if (jobrcvnum[i] > 0) {
			ubuf = jobs_to_rcv[i];
			count = jobrcvnum[i];
			src = nbrs[i];
			tag = JOBEXCH_TAG;
			req = &(rcvreqs[j++]);
			MPI_Irecv(ubuf, count, MPI_UNSIGNED, src, tag,
				MPI_COMM_WORLD, req);
		}
	}
	if (j > 0)
		MPI_Waitall(j, rcvreqs, MPI_STATUSES_IGNORE);
	if (k > 0)
		MPI_Waitall(k, sndreqs, MPI_STATUSES_IGNORE);

	*numtosnd = 0;
	*numtorcv = 0;
	for (i = 0; i < nbrsrtnum; i++)
		*numtosnd += jobsndnum[i];
	for (i = 0; i < nbrsnum; i++)
		*numtorcv += jobrcvnum[i];

	free(num_to_snd);
	free(sndreqs);
	free(rcvreqs);
}

void send_job(job_t *alljobs, int num, int dest, double **buf, MPI_Request *req,
		int rank)
{
	int i, linelen;
	int count, tag;
	double *dstptr;
	job_t *job;

	job = &(alljobs[num]);
	count = (job->xnodes + 2) * (job->ynodes + 2);
	count *= 2;
	tag = num;
	*buf = calloc(count, sizeof(double));

	linelen = job->xnodes + 2;
	for (i = 0; i < job->ynodes+2; i++) {
		dstptr = *buf + linelen*i;
		memcpy(dstptr, job->N[old][i], linelen*sizeof(double));
		dstptr = *buf + (count/2) + linelen*i;
		memcpy(dstptr, job->M[old][i], linelen*sizeof(double));
	}

	MPI_Isend(*buf, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);

	release_cell(job, rank);
}

void rcv_job(job_t *alljobs, int num, int src, double **buf, MPI_Request *req)
{
	int count, tag;
	job_t *job;

	job = &(alljobs[num]);
	count = (job->xnodes + 2) * (job->ynodes + 2);
	count *= 2;
	tag = num;
	*buf = calloc(count, sizeof(double));
	MPI_Irecv(*buf, count, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, req);
}

static void collect_incoming(job_t *alljobs, int num, double *buf)
{
	int i, linelen, count;
	double *srcptr;
	job_t *job;

	job = &(alljobs[num]);
	alloc_stor(job);
	linelen = job->xnodes + 2;
	count = (job->xnodes + 2) * (job->ynodes + 2);

	for (i = 0; i < job->ynodes+2; i++) {
		srcptr = buf + linelen*i;
		memcpy(job->N[old][i], srcptr, linelen*sizeof(double));
		srcptr = buf + count + linelen*i;
		memcpy(job->M[old][i], srcptr, linelen*sizeof(double));
	}
}

static void end_alloc(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, num;
	job_t *job;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		job = &(alljobs[num]);
		alloc_bords(job);
	}
}

static void transfer_jobs(job_t *alljobs, int *nbrs, uint nbrsnum, int *nbrs_srt,
		int nbrsrtnum, int *jobsndnum, uint **jobs_to_snd, int *jobrcvnum,
		uint **jobs_to_rcv, int rank)
{
	int i, j, k, m, num;
	int numtosnd, numtorcv;
	int src, dest;
	double **sndbuf, **rcvbuf;
	MPI_Request *sndreqs, *rcvreqs;

	numtosnd = 0;
	numtorcv = 0;
	for (i = 0; i < nbrsrtnum; i++)
		numtosnd += jobsndnum[i];
	for (i = 0; i < nbrsnum; i++)
		numtorcv += jobrcvnum[i];
	sndreqs = calloc(numtosnd, sizeof(MPI_Request));
	rcvreqs = calloc(numtorcv, sizeof(MPI_Request));
	sndbuf = calloc(numtosnd, sizeof(double *));
	rcvbuf = calloc(numtorcv, sizeof(double *));

	k = 0;
	m = 0;
	for (i = 0; i < nbrsrtnum; i++)
		for (j = 0; j < jobsndnum[i]; j++) {
			num = jobs_to_snd[i][j];
			dest = nbrs_srt[i];
			send_job(alljobs, num, dest, &(sndbuf[k]), &(sndreqs[k]),
					rank);
			k++;
		}
	for (i = 0; i < nbrsnum; i++)
		for (j = 0; j < jobrcvnum[i]; j++) {
			num = jobs_to_rcv[i][j];
			src = nbrs[i];
			rcv_job(alljobs, num, src, &(rcvbuf[m]), &(rcvreqs[m]));
			m++;
		}

	assert(k == numtosnd);
	assert(m == numtorcv);
	MPI_Waitall(numtosnd, sndreqs, MPI_STATUSES_IGNORE);
	for (i = 0; i < numtosnd; i++)
		free(sndbuf[i]);
	MPI_Waitall(numtorcv, rcvreqs, MPI_STATUSES_IGNORE);
	m = 0;
	/*
	if (rank == 3) {
		printf("3th rcving: ");
		for (i = 0; i < nbrsnum; i++)
			for (j = 0; j < jobrcvnum[i]; j++)
				printf("%d, ", jobs_to_rcv[i][j]);
		printf("; and snding: ");
		for (i = 0; i < nbrsrtnum; i++)
			for (j = 0; j < jobsndnum[i]; j++)
				printf("%d, ", jobs_to_snd[i][j]);
		printf("\n");
		fflush(stdout);
	}*/
	for (i = 0; i < nbrsnum; i++)
		for (j = 0; j < jobrcvnum[i]; j++) {
			num = jobs_to_rcv[i][j];
			collect_incoming(alljobs, num, rcvbuf[m]);
			free(rcvbuf[m]);
			m++;
		}

	free(sndreqs);
	free(rcvreqs);
	free(sndbuf);
	free(rcvbuf);
}

static void renew_actjobs(job_t *alljobs, uint **activejobs, uint *actjobsnum,
		int *jobsndnum, uint **jobs_to_snd, int *jobrcvnum,
		uint **jobs_to_rcv, int nbrsrtnum, int nbrsnum,
		int numtosnd, int numtorcv, int rank)
{
	int i, j, cnt, k, num;

	cnt = 0;
	for (i = 0; i < nbrsrtnum; i++)
		for (j = 0; j < jobsndnum[i]; j++) {
			num = jobs_to_snd[i][j];
			for (k = 0; k < *actjobsnum; k++)
				if ((*activejobs)[k] == num)
					break;
			assert(k < *actjobsnum);
			(*activejobs)[k] = UINT_MAX;
			cnt++;
		}
	assert(cnt == numtosnd);

	for (i = 0; i < *actjobsnum-numtosnd; i++)
		if((*activejobs)[i] == UINT_MAX)
			for (j = i+1; j < *actjobsnum; j++)
				if ((*activejobs)[j] != UINT_MAX)
					swap((*activejobs)[i], (*activejobs)[j]);
	for (j = i; j < *actjobsnum; j++)
		assert((*activejobs)[j] == UINT_MAX);

	k = *actjobsnum - numtosnd;
	*actjobsnum = *actjobsnum - numtosnd + numtorcv;
	*activejobs = realloc(*activejobs, (*actjobsnum)*sizeof(uint));

	cnt = 0;
	for (i = 0; i < nbrsnum; i++)
		for (j = 0; j < jobrcvnum[i]; j++) {
			num = jobs_to_rcv[i][j];
			(*activejobs)[k] = num;
			k++;
		}
	assert(k == *actjobsnum);
}

void rebalance(int rank, int np, job_t *alljobs, uint alljobsnum, uint **activejobs,
		uint *actjobsnum, int *nbrs, uint nbrsnum, int *jobsmap,
		MPI_Request *wloadreqs)
{
	int i, *nbrs_srt, nbrsrtnum, nbr_rank, jobsrtnum, *jobsndnum,
		*jobrcvnum, numtosnd, numtorcv;
	ulong total_wload, *wloads;
	ushort *jbusymap;
	long *wl_diffs, *wldiffs_srt;
	uint *jobs_srt, **jobs_to_snd, **jobs_to_rcv;

	norm_wloads(alljobs, *activejobs, *actjobsnum);
	total_wload = sum_avg_wloads(alljobs, *activejobs, *actjobsnum);
	wloads = share_wloads(rank, alljobs, *activejobs, *actjobsnum,
			nbrs, nbrsnum, wloadreqs, total_wload);
	wl_diffs = calloc(nbrsnum, sizeof(long));

	/*
	if (rank == 0) {
		for (i = 0; i < nbrsnum; i++)
			printf("%lu, ", wloads[i]);
		printf("%lu\n", total_wload);
		fflush(stdout);
	}
	*/
	/*
	sum = total_wload;
	for (i = 0; i < nbrsnum; i++)
		sum += wloads[i];
	avg_distr = (double)sum / (nbrsnum + 1);
	deviation = (double)total_wload - avg_distr;
	*/

	for (i = 0; i < nbrsnum; i++)
		wl_diffs[i] = (long)total_wload - (long)wloads[i];
	free(wloads);
	nbrsrtnum = sort_wl_diffs(nbrs, nbrsnum, wl_diffs,
			&wldiffs_srt, &nbrs_srt, rank);
	free(wl_diffs);

	jobs_to_snd = calloc(nbrsrtnum, sizeof(uint *));
	jobsndnum = calloc(nbrsrtnum, sizeof(int));
	jbusymap = calloc(alljobsnum, sizeof(ushort));
	for (i = 0; i < nbrsrtnum; i++) {
		nbr_rank = nbrs_srt[i];
		jobsrtnum = get_srtjobs(alljobs, *activejobs, *actjobsnum,
				&jobs_srt, nbr_rank);
		jobsndnum[i] = choose_to_send(alljobs, jobs_srt, jobsrtnum, nbrsnum,
				&(jobs_to_snd[i]), wldiffs_srt[i], jbusymap);
		free(jobs_srt);
	}
	free(jbusymap);
	free(wldiffs_srt);

	jobs_to_rcv = calloc(nbrsnum, sizeof(uint *));
	jobrcvnum = calloc(nbrsnum, sizeof(int));
	get_to_rcv(nbrs, nbrsnum, nbrs_srt, nbrsrtnum, jobsndnum, jobrcvnum,
			jobs_to_snd, jobs_to_rcv, &numtosnd, &numtorcv, rank);
	if (numtosnd > 0 || numtorcv > 0) {
		transfer_jobs(alljobs, nbrs, nbrsnum, nbrs_srt, nbrsrtnum, jobsndnum,
				jobs_to_snd, jobrcvnum, jobs_to_rcv, rank);
		renew_actjobs(alljobs, activejobs, actjobsnum, jobsndnum, jobs_to_snd,
				jobrcvnum, jobs_to_rcv, nbrsrtnum, nbrsnum, numtosnd,
				numtorcv, rank);
	}
	map_nbrs(rank, np, alljobs, alljobsnum, *activejobs, *actjobsnum,
			jobsmap);
	end_alloc(alljobs, *activejobs, *actjobsnum);

	for (i = 0; i < nbrsrtnum; i++)
		free(jobs_to_snd[i]);
	for (i = 0; i < nbrsnum; i++)
		free(jobs_to_rcv[i]);

	free(nbrs_srt);
	free(jobs_to_snd);
	free(jobsndnum);
	free(jobs_to_rcv);
	free(jobrcvnum);
}

void renew_resources(int rank, int np, job_t *alljobs, uint *activejobs,
		uint actjobsnum, MPI_Request **sharereqs, MPI_Request **wloadreqs,
		uint *nbredgenum, int **nbrs, uint *nbrsnum)
{
	*nbredgenum = count_nbredges(alljobs, activejobs, actjobsnum, rank);
	free(*sharereqs);
	*sharereqs = prep_shrreqs(*nbredgenum);
	*nbrs = get_nbrs(rank, np, alljobs, activejobs, actjobsnum, nbrsnum);
	free(*wloadreqs);
	*wloadreqs = prep_wldreqs(*nbrsnum);
}

void release_cell(job_t *job, int rank)
{
	int j, k;
	/*
	if (job->num == 43 && rank == 3) {
		printf("HEYYA!! %p, %p, %p\n", job->N[0], job->N[1], job->N[2]);
		fflush(stdout);
	}*/

	assert(job->N != NULL);
	assert(job->M != NULL);

	for (j = 0; j < 3; j++) {
		assert(job->N[j] != NULL);
		assert(job->M[j] != NULL);
		for (k = 0; k < job->ynodes + 2; k++) {
			assert(job->N[j][k] != NULL);
			assert(job->M[j][k] != NULL);
			free(job->N[j][k]);
			free(job->M[j][k]);
		}
		assert(job->N[j] != NULL);
		assert(job->M[j] != NULL);
		free(job->N[j]);
		free(job->M[j]);
	}

	if (comm_lft_bord(job)) {
		assert(job->clms_snd[0] != NULL);
		assert(job->clms_rcv[0] != NULL);
		free(job->clms_snd[0]);
		free(job->clms_rcv[0]);
	}
	if (comm_ryt_bord(job)) {
		assert(job->clms_snd[1] != NULL);
		assert(job->clms_rcv[1] != NULL);
		free(job->clms_snd[1]);
		free(job->clms_rcv[1]);
	}
	if (comm_low_bord(job)) {
		assert(job->rows_snd[0] != NULL);
		assert(job->rows_rcv[0] != NULL);
		free(job->rows_snd[0]);
		free(job->rows_rcv[0]);
	}
	if (comm_top_bord(job)) {
		assert(job->rows_snd[1] != NULL);
		assert(job->rows_rcv[1] != NULL);
		free(job->rows_snd[1]);
		free(job->rows_rcv[1]);
	}
}

void free_resources(job_t *alljobs, uint *activejobs, uint actjobsnum,
		MPI_Request *sharereqs, MPI_Request *wloadreqs, int *nbrs,
		int *jobsmap, int rank)
{
	int i, num;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		assert(&(alljobs[num]) != NULL);
		release_cell(&(alljobs[num]), rank);
	}

	assert(activejobs != NULL);
	assert(alljobs != NULL);
	assert(sharereqs != NULL);
	assert(wloadreqs != NULL);
	assert(nbrs != NULL);
	assert(jobsmap != NULL);

	free(activejobs);
	free(alljobs);
	free(sharereqs);
	free(wloadreqs);
	free(nbrs);
	free(jobsmap);
}
