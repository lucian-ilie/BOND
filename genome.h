#ifndef _GENOME_
#define _GENOME_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdint.h>
#include <omp.h>

typedef struct {
	int nbytes; // amount of bytes for coding genome
	int gSize; // size of genome
	uint8_t * gCoded;
} genome;
extern unsigned int PRIME;
extern unsigned int PRIME_IDEN;
extern unsigned int FRIME;
extern int TEST_PARAM;
extern int POS_SIZE;
extern int BLOCK;
extern int KANE_SIMILARITY;
extern int KANE_IDENTITY;
extern int OLIG_LENGTH;
extern int DELTATM;
extern int GC_CONT_MIN;
extern int GC_CONT_MAX;
extern int geneno;
extern int m;
extern uint16_t preNOR[65536];
extern int MAX_NUM_THREAD;
void pre_Process();
void gc_Content(genome *, uint8_t *);
void similarity_Check(genome *, uint8_t *, const char *);
void identity_Check(genome *, uint8_t *);
void melting_Temp(genome *, uint8_t *, int *, float *);
unsigned long print_Oligo_S8W9(genome *, uint8_t *, unsigned char *, int *, FILE *, int, float *);
void destroy_Genome(genome *);
int gene_find(int , int *);
void print_Time();
#endif
