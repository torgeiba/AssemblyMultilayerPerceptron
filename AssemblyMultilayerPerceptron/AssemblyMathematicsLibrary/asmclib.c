#include "asmclib.h"

#include "stdlib.h" // malloc, calloc, 

#include "windows.h"

// Slow approximation to exp function
#if 0
float exp_asm(float x) {
	int iterations = 50;

	// Compute iteger exponent part
	/*
	double xlog_2e = x * 1.4426950408889634;
	double abs_xlog_2e = xlog_2e >= 0 ? xlog_2e : -xlog_2e;
	double integer_part = 1.;
	if (xlog_2e > 0.)		{ integer_part = (double)(1ULL << (uint64)abs_xlog_2e);		}
	else if (xlog_2e < 0.)  { integer_part = 1. / ((double)(1ULL << (uint64)abs_xlog_2e)); }
	*/

	// Taylor series approximation. TODO: Cosider computing integer part of exponent as outlined above
	double f = x;// -(int)x;
	double fractional_part = 1.;
	for (int n = iterations; n > 0; n--)
	{
		fractional_part = 1. + (f / (double)n) * fractional_part;
	}
	double res = fractional_part;
	return (float)(res);
}

void* malloc_asm(uint64 size) {
	return HeapAlloc(GetProcessHeap(), 0, size);
}

void* calloc_asm(uint64 count, uint64 size) {
	return HeapAlloc(GetProcessHeap(), 8, size*count); // HEAP_ZERO_MEMORY flag = 0x00000008
}

void free_asm(void* block_ptr) {
	HeapFree(GetProcessHeap(), 0, block_ptr); 
}

#endif


//int rand_asm() { return rand(); }
/*
int rand_asm()
{ 
	static unsigned int rnd_state = 7777;
	rnd_state ^= rnd_state << 13;
	rnd_state ^= rnd_state >> 17;
	rnd_state ^= rnd_state << 5;
	return (int)(rnd_state & 0x7fff);
}
*/