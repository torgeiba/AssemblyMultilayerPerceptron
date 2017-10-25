#include "asmclib.h"

#include "math.h"   // expf, sinf, cosf
#include "stdlib.h" // malloc, calloc, 

float exp_asm(float x) { return expf(x); }
void* malloc_asm(uint64 size) { return malloc(size); }
void* calloc_asm(uint64 count, uint64 size) { return calloc(count, size); }
void free_asm(void* block_ptr) { free(block_ptr); }
int rand_asm() { return rand(); }
