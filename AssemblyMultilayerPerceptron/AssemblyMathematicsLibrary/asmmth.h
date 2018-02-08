#pragma once

#include "asmclib.h" // include c library helper functions

/*
* x64 Assembly math routines
*/

typedef unsigned __int64	uint64;
typedef unsigned __int32	uint32;
typedef unsigned __int16	uint16;
typedef unsigned __int8		uint8;

typedef __int64	int64;
typedef __int32	int32;
typedef __int16	int16;
typedef __int8	int8;


#ifndef VECTYPES
	#define VECTYPES
struct vec {
	float* data;
	uint64 size;
};
typedef struct vec vec;

struct mat {
	float** data;
	uint64 rows;
	uint64 cols;
};
typedef struct mat mat;

#endif

extern void* malloc_asm(uint64 size);
extern void* calloc_asm(uint64 count, uint64 size);
extern void free_asm(void* block_ptr);

/**
* Basic uninitialized vec memory allocation
* @param
* @return
*/
extern vec vecalloc_asm(/* result pointer (hidden),*/ uint64 size); // Completed
extern void vecfree_asm(vec* u); // Completed

/**
* Basic uninitialized mat memory allocation
* @param
* @return
*/
extern mat matalloc_asm(/* result pointer (hidden),*/ uint64 rows, uint64 cols); // Completed
extern void matfree_asm(mat* m); // Completed

/**
* Make zero-initialized array of given size
* @param
* @return
*/
extern vec veczeros_asm(/* result pointer (hidden),*/ uint64 size); // Completed

/**
* Make one-initialized array of given size
* @param
* @return
*/
extern vec vecones_asm(/* result pointer (hidden),*/ uint64 size); // Completed

/**
* Make zero-initialized 2D-array
* @param
* @return
*/
extern mat matzeros_asm(/* result pointer (hidden),*/ uint64 rows, uint64 cols); // Completed

/**
* Make one-initialized 2D-array
* @param
* @return
*/
extern mat matones_asm(/* result pointer (hidden),*/ uint64 rows, uint64 cols); // Completed

/**
* Make idendity-initialized 2D-array
* @param
* @return
*/
extern mat matidentity_asm(/* result pointer (hidden),*/ uint64 rows, uint64 cols); // Completed

extern vec vecrand_asm(/* result pointer (hidden),*/ uint64 size); // Completed

extern mat matrand_asm(/* result pointer (hidden),*/ uint64 rows, uint64 cols); // Completed

extern vec* veccopy_asm(vec* fromvec, vec* tovec); // Completed
extern mat* matcopy_asm(mat* frommat, mat* tomat); // Completed

/**
* Generate an array of linearly spaced values from an interval.
* @param
* @return
*/
// TODO: Handle 0 steps
extern vec linspace_asm(/* result pointer (hidden),*/ float start, float end, uint64 steps); // Completed


extern float vecdotvec_asm(vec* u, vec* v); // Completed


extern vec* matdotvec_asm(mat* m, vec* v, vec* result); // Completed, Buggy =( 
extern mat* matdotmat_asm(mat* m, mat* n, mat* result); // Completed

// Unsafe, but faster, copyless functions fail if v and result points to the same or overlapping location
extern vec* matdotvec_unsafe_asm(mat* m, vec* v, vec* result); // Completed
extern vec* matTdotvec_unsafe_asm(mat* m, vec* v, vec* result); // Completed

extern vec* matTdotvec_asm(mat* m, vec* v, vec* result); // Completed

extern mat* outervec_asm(vec* u, vec* v, mat* result);  // Completed

extern vec* addvec_asm(vec* u, vec* v, vec* result);    // Completed
extern vec* subvec_asm(vec* u, vec* v, vec* result);	// Completed
extern vec* mulvec_asm(vec* u, vec* v, vec* result);	// Completed
extern vec* divvec_asm(vec* u, vec* v, vec* result);	// Completed
extern vec* scalevec_asm(float a, vec* v, vec* result); // Completed

extern vec* expvec_asm(vec* u, vec* result);			// Completed

extern uint64 equalsvec_asm(vec* u, vec* v);			// Completed

extern mat* addmat_asm(mat* m, mat* n, mat* result);    // Completed
extern mat* submat_asm(mat* m, mat* n, mat* result);	// Completed
extern mat* mulmat_asm(mat* m, mat* n, mat* result);	// Completed
extern mat* divmat_asm(mat* m, mat* n, mat* result);	// Completed
extern mat* scalemat_asm(float a, mat* m, mat* result);	// Completed

extern mat* expmat_asm(mat* u, mat* result);			// Completed

extern uint64 equalsmat_asm(mat* m, mat* n);			// Completed
extern void transposemat_asm(mat* m, mat* result);		// Completed