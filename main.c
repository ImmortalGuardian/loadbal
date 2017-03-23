#include <mpi.h>
#include <loadbal.h>
#include <assist.h>

int main (int argc, char *argv[])
{
	int rank, np;
	unsigned int alljobsnum, actjobsnum;
	job_t *alljobs;
	unsigned int *activejobs;
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
	activejobs = distr_jobs(rank, np, alljobsnum, &actjobsnum);

	errnum =  MPI_Finalize();
	if (errnum != MPI_SUCCESS)
		PRERROR("main: cannot finalize: ", errnum);
	return EXIT_SUCCESS;
}
