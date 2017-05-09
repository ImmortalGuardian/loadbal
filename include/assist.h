#ifndef ASSIST_H
#define ASSIST_H

#include <problem.h>
#include <loadbal.h>

extern void assist_init(void);
extern job_t *form_jobs(uint np, uint *alljobsnum);

extern void set_init_cond(job_t *alljobs, uint *activejobs, uint actjobsnum);
extern void make_timestep(job_t *alljobs, uint *activejobs, uint actjobsnum,
		MPI_Request *sharereqs, uint nbredgenum);

// int which_proc(int cellnum, int np, uint alljobsnum);
void map_nbrs(int rank, int np, job_t *alljobs, uint alljobsnum, uint *activejobs,
		uint actjobsnum, int *jobsmap);
uint *distr_jobs(int rank, int np, job_t *alljobs, uint alljobsnum,
		uint *actjobsnum, int *jobsmap);
uint count_nbredges(job_t *alljobs, uint *activejobs, uint actjobsnum, int rank);
int *get_nbrs(int rank, int np, job_t *alljobs, uint *activejobs,
		uint actjobsnum, uint *nbrsnum);
MPI_Request *prep_shrreqs(uint nbredgenum);
MPI_Request *prep_wldreqs(uint nbrsnum);
void alloc_memory(job_t *alljobs, uint *activejobs, uint actjobsnum);
void free_resources(job_t *alljobs, uint *activejobs, uint actjobsnum,
		MPI_Request *sharereqs, MPI_Request *wloadreqs);

#endif // ASSIST_H
