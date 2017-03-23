#include <assist.h>
#include <problem.h>

unsigned int
*distr_jobs(int rank, int np, unsigned int alljobsnum, unsigned int *actjobsnum)
{
	unsigned int rem, i, currnum, jobsperproc;
	unsigned int *activejobs;

	jobsperproc = alljobsnum / np;
	*actjobsnum = jobsperproc;
	rem = alljobsnum % np;
	if (rank < rem)
		*actjobsnum += 1;
	activejobs = (unsigned int *)calloc(*actjobsnum, sizeof(int));

	currnum = rank;
	for (i = 0; i < *actjobsnum; i++) {
		activejobs[i] = currnum;
		currnum += np;
	}

	return activejobs;
}
