#include "procmth.h"

#include "stdlib.h"

#ifndef max
#define max(a, b) ((a)>(b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b) ((a)<(b) ? (a) : (b))
#endif


vec vecalloc(uint64 size)
{
	vec result;
	result.size = size;
	result.data = (float*) malloc(size * sizeof(float));
	return result;
}
void vecfree(vec u) { u.size = 0; free(u.data); }

mat matalloc(uint64 rows, uint64 cols)
{
	mat result;
	result.rows = rows;
	result.cols = cols;
	result.data = (float**) malloc(rows * sizeof(float*));
	result.data[0] = (float*) malloc(rows * cols * sizeof(float));
	uint64 row = 1;
	while (row < rows) {
		result.data[row] = result.data[row - 1] + cols;
		row++;
	}
	return result;
}
void matfree(mat m)
{
	m.rows = 0;
	m.cols = 0;
	free(m.data[0]);
	free(m.data);
}

vec veczeros( uint64 size)
{
	vec result;
	result.size = size;
	result.data = (float*)calloc(size, sizeof(float));
	return result;
}


vec vecones(uint64 size)
{
	vec result = vecalloc(size);
	uint64 index = 0;
	while (index < result.size)
	{
		result.data[index] = 1.f;
		index++;
	}
	return result;
}

mat matzeros(uint64 rows, uint64 cols)
{
	mat result;
	result.rows = rows;
	result.cols = cols;
	result.data = (float**)calloc(rows, sizeof(float*));
	result.data[0] = (float*)calloc(rows * cols, sizeof(float));
	uint64 row = 1;
	while (row < rows) {
		result.data[row] = result.data[row - 1] + cols;
		row++;
	}
	return result;
}

mat matones(uint64 rows, uint64 cols)
{
	mat result = matalloc(rows, cols);
	uint64 r = 0;
	while (r < result.rows)
	{
		uint64 c = 0;
		while (c < result.cols)
		{
			result.data[r][c] = 1.f;
			c++;
		}
		r++;
	}
	return result;
}

vec linspace(float start, float end, uint64 steps)
{
	vec result = vecalloc(steps);
	if (steps <= 0ULL) { return result; }
	else if (steps == 1ULL) { 
		result.data[0] = start; return result;
	}
	double stepsize = ((double)end - (double)start) / ((double)steps - 1.);
	uint64 index = 0;
	while (index < steps) {
		result.data[index] = (float)(((double)index * stepsize) + (double)start);
		index++;
	}
	return result;
}

float vecdotvec(vec u, vec v)
{
	float result = 0.f;
	uint64 size = min(u.size, v.size);// u.size < v.size ? u.size : v.size;
	uint64 index = 0;
	while (index < size)
	{
		result += u.data[index] * v.data[index];
		index++;
	}
	return result;
}

vec matdotvec(mat m, vec v, vec result)
{
	//vec result = vecalloc(m.rows);
	vec row;
	row.size = m.cols;
	uint64 index = 0ULL;
	vec tempvec = vecalloc(m.rows);
	while(index < m.rows)
	{
		row.data = m.data[index];
		tempvec.data[index] = vecdotvec(row, v);
		index++;
	}
	veccopy(tempvec, result);
	return result;
}

vec matTdotvec(mat m, vec v, vec result)
{
	//vec result = veczeros(m.cols);
	uint64 size = min(v.size, m.rows); //v.size < m.rows ? v.size : m.rows;
	uint64 c = 0ULL;
	vec tempvec = veczeros(result.size);
	while (c < m.cols)
	{
		uint64 r = 0ULL;
		while (r < size)
		{
			tempvec.data[c] += m.data[r][c] * v.data[r];
			r++;
		}
		c++;
	}
	veccopy(tempvec, result);
	return result;
}

mat matdotmat(mat m, mat n, mat result)
{
	//mat result = matzeros(m.rows, n.cols);
	uint64 size = min(m.cols, n.rows); // m.cols < n.rows ? m.cols : n.rows;
	uint64 r = 0ULL;
	while (r < result.rows)
	{
		uint64 c = 0ULL;
		while (c < result.cols)
		{
			uint64 i = 0ULL;
			while (i < size)
			{
				result.data[r][c] += m.data[r][i] *  n.data[i][c];
				i++;
			}
			c++;
		}
		r++;
	}
	return result;
}

mat outervec(vec u, vec v, mat result)
{
	//mat result = matzeros(u.size, v.size);
	uint64 r = 0ULL;
	vec outrow;
	outrow.size = result.cols;
	while (r < result.rows)
	{
		outrow.data = result.data[r];
		scalevec(u.data[r], v, outrow);
		/*
		uint64 c = 0ULL;
		while (c < result.cols)
		{
			result.data[r][c] = u.data[r] * v.data[c];
			c++;
		}*/
		r++;
	}
	return result;
}

vec addvec(vec u, vec v, vec result)
{
	uint64 size = min(u.size, v.size); //u.size < v.size ? u.size : v.size;
	//vec result = vecalloc(u.size);
	uint64 index = 0ULL;
	while (index < size)
	{
		result.data[index] = u.data[index] + v.data[index];
		index++;
	}
	return result;
}
vec subvec(vec u, vec v, vec result)
{
	uint64 size = min(u.size, v.size); //u.size < v.size ? u.size : v.size;
	//vec result = vecalloc(u.size);
	uint64 index = 0ULL;
	while (index < size)
	{
		result.data[index] = u.data[index] - v.data[index];
		index++;
	}
	return result;
}
vec mulvec(vec u, vec v, vec result)
{
	uint64 size = min(u.size, v.size); //u.size < v.size ? u.size : v.size;
	//vec result = vecalloc(u.size);
	uint64 index = 0ULL;
	while (index < size)
	{
		result.data[index] = u.data[index] * v.data[index];
		index++;
	}
	return result;
}
vec divvec(vec u, vec v, vec result)
{
	uint64 size = min(u.size, v.size); //u.size < v.size ? u.size : v.size;
	//vec result = vecalloc(u.size);
	uint64 index = 0ULL;
	while (index < size)
	{
		result.data[index] = u.data[index] / v.data[index];
		index++;
	}
	return result;
}
vec scalevec(float a, vec u, vec result)
{
	//vec result = vecalloc(u.size);
	uint64 index = 0ULL;
	while (index < u.size)
	{
		result.data[index] = a * u.data[index];
		index++;
	}
	return result;
}

mat addmat(mat m, mat n, mat result)
{
	uint64 rows = min(m.rows, n.rows); // m.rows < n.rows ? m.rows: n.rows;
	uint64 cols = min(m.cols, n.cols); // m.cols < n.cols ? m.cols : n.cols;
	//mat result = matalloc(rows, cols);
	uint64 r = 0ULL;
	while (r < rows)
	{
		uint64 c = 0ULL;
		while (c < cols)
		{
			result.data[r][c] = m.data[r][c] + n.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}
mat submat(mat m, mat n, mat result)
{
	uint64 rows = min(m.rows, n.rows); // m.rows < n.rows ? m.rows : n.rows;
	uint64 cols = min(m.cols, n.cols); // m.cols < n.cols ? m.cols : n.cols;
	//mat result = matalloc(rows, cols);
	uint64 r = 0ULL;
	while (r < rows)
	{
		uint64 c = 0ULL;
		while (c < cols)
		{
			result.data[r][c] = m.data[r][c] - n.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}
mat mulmat(mat m, mat n, mat result)
{
	uint64 rows = min(m.rows, n.rows); // m.rows < n.rows ? m.rows : n.rows;
	uint64 cols = min(m.cols, n.cols); // m.cols < n.cols ? m.cols : n.cols;
	//mat result = matalloc(rows, cols);
	uint64 r = 0ULL;
	while (r < rows)
	{
		uint64 c = 0ULL;
		while (c < cols)
		{
			result.data[r][c] = m.data[r][c] * n.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}
mat divmat(mat m, mat n, mat result)
{
	uint64 rows = min(m.rows, n.rows); // m.rows < n.rows ? m.rows : n.rows;
	uint64 cols = min(m.cols, n.cols); // m.cols < n.cols ? m.cols : n.cols;
	//mat result = matalloc(rows, cols);
	uint64 r = 0ULL;
	while (r < rows)
	{
		uint64 c = 0ULL;
		while (c < cols)
		{
			result.data[r][c] = m.data[r][c] / n.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}
mat scalemat(float a, mat m, mat result)
{
	//mat result = matalloc(m.rows, m.cols);
	uint64 r = 0ULL;
	while (r < m.rows)
	{
		uint64 c = 0ULL;
		while (c < m.cols)
		{
			result.data[r][c] = a * m.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}

//mat transposemat(mat m, mat result)
//{
//	//mat result = matalloc(m.cols, m.rows);
//	uint64 i = 0ULL;
//	while (i < m.rows)
//	{
//		uint64 j = 0ULL;
//		while (j < m.cols)
//		{
//			// Does not work for m =  result!
//			result.data[j][i] = m.data[i][j];
//			j++;
//		}
//		i++;
//	}
//	return result;
//}

float frand()
{
	return (float)((double)rand() / (double)RAND_MAX);
}

vec vecrand(uint64 size)
{
	vec result = vecalloc(size);
	uint64 index = 0;
	while (index < result.size)
	{
		result.data[index] = frand();
		index++;
	}
	return result;
}

mat matrand(uint64 rows, uint64 cols)
{
	mat result = matalloc(rows, cols);
	uint64 r = 0;
	while (r < result.rows)
	{
		uint64 c = 0;
		while (c < result.cols)
		{
			result.data[r][c] = frand();
			c++;
		}
		r++;
	}
	return result;
}

mat matidendity(uint64 rows, uint64 cols)
{
	mat result = matzeros(rows, cols);
	uint64 size = min(rows, cols);
	uint64 i = 0;
	while (i < size)
	{
		result.data[i][i] = 1.f;
		i++;
	}
	return result;
}

vec veccopy(vec fromvec, vec tovec)
{
	uint64 size = min(fromvec.size, tovec.size);
	uint64 index = 0;
	while (index < size)
	{
		tovec.data[index] = fromvec.data[index];
 		index++;
	}
	return tovec;
}
mat matcopy(mat m)
{
	mat result = matalloc(m.rows, m.cols);
	uint64 r = 0;
	uint64 c = 0;
	while (r < m.rows)
	{
		c = 0;
		while (c < m.cols)
		{
			result.data[r][c] = m.data[r][c];
			c++;
		}
		r++;
	}
	return result;
}

uint64 equalsvec(vec u, vec v)
{
	if (u.size != v.size) return 0ULL;
	
	uint64 index = 0;
	while (index < u.size)
	{
		if (u.data[index] != v.data[index]) return 0ULL;
		index++;
	}
	return 1ULL;
}

uint64 equalsmat(mat m, mat n)
{
	if ((m.rows != n.rows) || (m.cols != n.cols)) return 0ULL;

	uint64 r = 0;
	uint64 c = 0;
	while (r < m.rows)
	{
		while (c < m.cols)
		{
			if (m.data[r][c] != n.data[r][c]) return 0ULL;
			c++;
		}
		r++;
	}
	return 1ULL;
}