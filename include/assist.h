#ifndef ASSIST_H
#define ASSIST_H

#include <problem.h>
#include <loadbal.h>

extern void assist_init(void);
extern job_t *form_jobs(uint np, uint *alljobsnum);
uint *distr_jobs(int rank, int np, uint alljobsnum,
		uint *actjobsnum);
extern void alloc_memory(job_t *alljobs, uint *activejobs, uint actjobsnum);
extern void set_init_cond(job_t *alljobs, uint *activejobs, uint actjobsnum);

#endif // ASSIST_H
