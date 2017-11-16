#include "procmthutils.h"

#include "stdio.h"

void printrow(vec u)
{
	printf("[");
	for (int i = 0; i + 1 < u.size; i++)
	{
		printf("%f, ", u.data[i]);
	}
	if (u.size > 0) printf("%f", u.data[u.size - 1]);
	printf("]");
}

void printcol(vec u)
{
	printf("[");
	for (int i = 0; i + 1 < u.size; i++)
	{
		printf("%f,\n", u.data[i]);
	}
	if (u.size > 0) printf("%f", u.data[u.size - 1]);
	printf("]");
}

void printvec(vec u) { printcol(u); printf("\n"); }



void printmat(mat m)
{
	vec row;
	row.size = m.cols;
	printf("[\n");
	for (int i = 0; i + 1 < m.rows; i++)
	{
		row.data = m.data[i];
		printrow(row);
		printf(",\n");
	}
	if (m.rows > 0) {
		row.data = m.data[m.rows - 1];
		printrow(row);
	}
	printf("\n]\n");
}