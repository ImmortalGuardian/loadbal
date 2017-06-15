#include <mpi.h>
#include <loadbal.h>
#include <assist.h>

int main (int argc, char *argv[])
{
	int rank, np, *nbrs, j, *jobsmap;
	uint alljobsnum, actjobsnum;
	uint nbredgenum, nbrsnum;
	job_t *alljobs;
	uint *activejobs;
	MPI_Request *sharereqs, *wloadreqs;
	int errnum;

	errnum = MPI_Init(&argc, &argv);
	if (errnum != MPI_SUCCESS)
		PRERROR("main: cannot initialize: ", errnum);
	errnum = MPI_Comm_size(MPI_COMM_WORLD, &np);
	if (errnum != MPI_SUCCESS)
		PRERROR("main: cannot get size: ", errnum);
	errnum = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (errnum != MPI_SUCCESS)
		PRERROR("main: cannot get process rank: ", errnum);

	assist_init();
	alljobs = form_jobs(np, &alljobsnum);
	jobsmap = calloc(alljobsnum, sizeof(int));
	activejobs = distr_jobs1(rank, np, alljobs, alljobsnum, &actjobsnum,
			jobsmap);
	/*
	printf("%d(%d): ", rank, actjobsnum);
	for (j = 0; j < actjobsnum; j++)
		printf("%d, ", activejobs[j]);
	printf("\n");
	fflush(stdout);
	*/

	alloc_memory(alljobs, activejobs, actjobsnum);
	set_init_cond(alljobs, activejobs, actjobsnum);
	nbredgenum = count_nbredges(alljobs, activejobs, actjobsnum, rank);
	sharereqs = prep_shrreqs(nbredgenum);
	nbrs = get_nbrs(rank, np, alljobs, activejobs, actjobsnum, &nbrsnum);
	wloadreqs = prep_wldreqs(nbrsnum);

	for (j = 1; j <= 10; j++)
		make_timestep(alljobs, activejobs, actjobsnum, sharereqs,
				nbredgenum, rank);
	rebalance(rank, np, alljobs, alljobsnum, &activejobs, &actjobsnum, nbrs,
			nbrsnum, jobsmap, wloadreqs);
	renew_resources(rank, np, alljobs, activejobs, actjobsnum, &sharereqs,
			&wloadreqs, &nbredgenum, &nbrs, &nbrsnum);
	for (j = 1; j <= 10; j++)
		make_timestep(alljobs, activejobs, actjobsnum, sharereqs,
				nbredgenum, rank);
	rebalance(rank, np, alljobs, alljobsnum, &activejobs, &actjobsnum, nbrs,
			nbrsnum, jobsmap, wloadreqs);
	renew_resources(rank, np, alljobs, activejobs, actjobsnum, &sharereqs,
			&wloadreqs, &nbredgenum, &nbrs, &nbrsnum);

	free_resources(alljobs, activejobs, actjobsnum, sharereqs, wloadreqs,
			nbrs, jobsmap, rank);
	errnum =  MPI_Finalize();
	if (errnum != MPI_SUCCESS)
		PRERROR("main: cannot finalize: ", errnum);
	return EXIT_SUCCESS;
}
