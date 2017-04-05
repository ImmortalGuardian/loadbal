#include <assist.h>
#include <problem.h>

uint
*distr_jobs(int rank, int np, uint alljobsnum, uint *actjobsnum)
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

	return activejobs;
}
