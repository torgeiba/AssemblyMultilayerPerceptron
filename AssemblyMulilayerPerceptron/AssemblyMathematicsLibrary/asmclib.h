#pragma once

// Helper header for asmmth.h for providing implementations of standard c library functions

typedef unsigned long long uint64;

float exp_asm(float x);
void* malloc_asm(uint64 size);
void* calloc_asm(uint64 count, uint64 size);
void free_asm(void* block_ptr);
int rand_asm();
