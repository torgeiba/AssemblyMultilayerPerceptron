#pragma once

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
#endif // VECTYPES

vec vecalloc(uint64 size);
void vecfree(vec u);

mat matalloc(uint64 rows, uint64 cols);
void matfree(mat m);

vec veczeros(uint64 size);

vec vecones(uint64 size);

mat matzeros(uint64 rows, uint64 cols);

mat matones(uint64 rows, uint64 cols);

mat matidendity(uint64 rows, uint64 cols);

vec vecrand(uint64 size);

mat matrand(uint64 rows, uint64 cols);

vec veccopy(vec fromvec, vec tovec);
mat matcopy(mat m);

vec linspace(float start, float end, uint64 steps);


float vecdotvec(vec u, vec v);

vec matdotvec(mat m, vec v, vec result);
mat matdotmat(mat m, mat n, mat result);

vec matTdotvec(mat m, vec v, vec result);

mat outervec(vec u, vec v, mat result);

vec addvec(vec u, vec v, vec result);
vec subvec(vec u, vec v, vec result);
vec mulvec(vec u, vec v, vec result);
vec divvec(vec u, vec v, vec result);
vec scalevec(float a, vec v, vec result);

uint64 equalsvec(vec u, vec v);

mat addmat(mat m, mat n, mat result);
mat submat(mat m, mat n, mat result);
mat mulmat(mat m, mat n, mat result);
mat divmat(mat m, mat n, mat result);
mat scalemat(float a, mat m, mat result);

uint64 equalsmat(mat m, mat n);
//mat transposemat(mat m, mat result);
