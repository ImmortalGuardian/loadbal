#ifndef ASSIST_H
#define ASSIST_H

#include <problem.h>

extern void assist_init(void);
extern job_t *form_jobs(unsigned int np, unsigned int *alljobsnum);
unsigned int *distr_jobs(int rank, int np, unsigned int alljobsnum,
		unsigned int *actjobsnum);

#endif // ASSIST_H
