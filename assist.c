#include <assert.h>
#include <assist.h>
#include <problem.h>

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

void release_cell(job_t *job) {
	int j, k;

	for (j = 0; j < 3; j++) {
		for (k = 0; k < job->ynodes + 2; k++) {
			free(job->N[j][k]);
			free(job->M[j][k]);
		}
		free(job->N[j]);
		free(job->M[j]);
	}
}

void free_resources(job_t *alljobs, uint *activejobs, uint actjobsnum)
{
	int i, num;

	for (i = 0; i < actjobsnum; i++) {
		num = activejobs[i];
		release_cell(&(alljobs[num]));
	}

	free(activejobs);
	free(alljobs);
}
