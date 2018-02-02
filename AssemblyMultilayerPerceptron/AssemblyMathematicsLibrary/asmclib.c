#include "asmclib.h"

#include "math.h"   // expf, sinf, cosf
#include "stdlib.h" // malloc, calloc, 

float exp_asm(float x) { 
	return expf(x);
}
void* malloc_asm(uint64 size) { return malloc(size); }
void* calloc_asm(uint64 count, uint64 size) { return calloc(count, size); }
void free_asm(void* block_ptr) { free(block_ptr); }



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