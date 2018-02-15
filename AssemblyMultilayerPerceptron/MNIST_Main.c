#include "ProceduralMultilayerPerceptron\procmlp.h"
#include "AssemblyMultilayerPerceptron/asmmlp.h"
#include "AssemblyMathematicsLibrary\asmclib.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"

char* image_path = "train-images.idx3-ubyte";
char* label_path = "train-labels.idx1-ubyte";

#pragma warning(disable:4996)

void decimal_to_vec(vec* v, uint8 L)
{
	for (uint32 i = 0; i < 10; i++)
	{
		v->data[i] = 0.f;
		if (i == (uint32)L) v->data[i] = 1.0;
	}
}

uint32 vec_to_decimal(vec* v)
{
	uint32 current_result = 256;
	float current_absdiff = 2.f;
	for (uint32 i = 0; i < 10; i++)
	{
		float f = v->data[i];
		float absdiff = (float)fabs(1.f - f);
		if (absdiff < current_absdiff)
		{
			current_absdiff = absdiff;
			current_result = i;
		}
	}
	return current_result;
}

mat mnist_read_labels(char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp) {
		fseek(fp, 0, SEEK_END);				// Jump to the end of the file
		uint64 filesize = ftell(fp);		// Get the current byte offset in the file
		uint8 *buffer = (uint8*)malloc(filesize);
		rewind(fp);							// Jump back to the beginning of the file
		fread(buffer, filesize, 1, fp);
		fclose(fp);

		uint8 magic[4] = { buffer[3], buffer[2], buffer[1], buffer[0] };
		uint8 numItems[4] = { buffer[7], buffer[6], buffer[5], buffer[4] };
		uint64 size = (uint64)(*(uint32*)numItems);
		mat labels = matalloc_asm(size, 10);

		uint8 *L = &buffer[8];
		vec v = vecalloc_asm(10);
		for (uint64 i = 0; i < size; i++)
		{
			decimal_to_vec(&v, L[i]);
			vec d;
			d.size = 10;
			d.data = labels.data[i];
			veccopy_asm(&v, &d);
		}
		vecfree_asm(&v);
		free(buffer);
		return labels;
	}
	mat errmat = { 0 };
	return errmat;
}

mat mnist_read_images(char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp) {
		fseek(fp, 0, SEEK_END);          // Jump to the end of the file
		uint64 filesize = ftell(fp); // Get the current byte offset in the file
		uint8 *buffer = (uint8*)malloc(filesize);
		rewind(fp);                      // Jump back to the beginning of the file
		fread(buffer, filesize, 1, fp);
		fclose(fp);

		uint8 magic[4] = { buffer[4], buffer[3], buffer[2], buffer[1] };
		uint8 numItems[4] = { buffer[7], buffer[6], buffer[5], buffer[4] };
		uint8 rowsize_bytes[4] = { buffer[11], buffer[10], buffer[9], buffer[8] };
		uint8 colsize_bytes[4] = { buffer[15], buffer[14], buffer[13], buffer[12] };
		uint64 itemcnt = (uint64)(*(uint32*)numItems);
		uint64 rowsize = (uint64)(*(uint32*)rowsize_bytes);
		uint64 colsize = (uint64)(*(uint32*)colsize_bytes);
		uint64 imgsize = rowsize * colsize;
		mat images = matalloc_asm(itemcnt, imgsize);

		uint8 *L = &buffer[16];
		for (uint64 i = 0; i < itemcnt; i++)
		{
			for(uint64 j = 0; j < imgsize; j++)
			{
				images.data[i][j] = ((float)(*L)) / 255.f;
				L++;
			}
		}

		return images;
	}
	mat errmat = { 0 };
	return errmat;
}

#pragma warning(enable:4996)



int main(int argc, char** argv)
{
	uint64 sizes[3] = { 784, 100, 10 };
	mlp* nn = makemlp_asm(3, sizes, 20, 0.01f); // note argument switch for the two last arguments!

	mat inputs = mnist_read_images(image_path); //matalloc_asm(trainingsize, inputsize);
	mat outputs = mnist_read_labels(label_path); //matalloc_asm(trainingsize, outputsize);

	train_asm(nn, &inputs, &outputs);


	// Verify training and print results :

	uint32 passed = 0;
	uint32 failed = 0;
	uint32 total  = 0;

	vec v;
	v.size = 28 * 28;
	vec w;
	w.size = 10;
	for (uint32 i = 0; i < inputs.rows; i++)
	{
		v.data = inputs.data[i];
		w.data = outputs.data[i];
		fwd_asm(nn, &v);
		uint32 result = vec_to_decimal(&nn->layers[2]);
		uint32 answer = vec_to_decimal(&w);
		const char* outcome_ok = "OK";
		const char* outcome_miss = "Miss";
		uint32 outcome = result == answer;
		printf("Test\t %u: \tPredicted: %u | Actual: %u |", i, result, answer);
		if (outcome)
		{
			printf(outcome_ok); printf("\n");
			passed++;
		}
		else
		{
			printf(outcome_miss); printf("\n");
			failed++;
		}
		total++;
	}
	printf("\n");
	printf("Passed: %u\n", passed);
	printf("Failed: %u\n", failed);
	printf("Total: %u\n", total);
	printf("\n");
	printf("Percentage passed: %f%%", (100.f * (float)passed)/((float)total));
	printf("\n");

	matfree_asm(&inputs);
	matfree_asm(&outputs);
	freemlp_asm(nn);
	getchar();

	return 0;
}
