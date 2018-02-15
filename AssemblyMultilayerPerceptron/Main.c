#include "ProceduralMultilayerPerceptron\procmlp.h"
#include "AssemblyMultilayerPerceptron/asmmlp.h"
#include "AssemblyMathematicsLibrary\asmclib.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"

float target_function(float x)
{
	return sinf(x) + 0.f * 0.025f * cosf(28 * x);
	// return sinf(x) + 0.025 * cosf(28 * x);
	//return x < 0.f ? -1.f : 1.f;
}
#if 0
int main(int argc, char** argv)
{
	char inseed[10];
	printf("Input random seed: ");
	//gets_s(inseed, 10);
	int randseed = 7777;
	//int randseed = abs(atoi(inseed));

	srand(randseed);
	printf("Using random seed %d.\n\n", randseed);

	uint64 sizes[4] = { 1, 5, 10, 1 };
	//mlp* nn = makemlp(4, sizes, 0.2f, 800);
	mlp* nn = makemlp_asm(4, sizes, 800, 0.2f); // note argument switch for the two last arguments!

	printmlp(nn);

	uint64 trainingsize = 100;
	uint64 inputsize  = 1;
	uint64 outputsize = 1;

	mat inputs  = matalloc_asm(trainingsize, inputsize);
	mat outputs = matalloc_asm(trainingsize, outputsize);

	float start = -3.1415f;
	float stop  =  3.1415f;
	vec x = linspace_asm(start, stop, trainingsize);

	// Create inputs / outputs:
	uint64 i = 0;
	while (i < trainingsize)
	{
		inputs.data[i][0] = x.data[i];
		outputs.data[i][0] = target_function(x.data[i]);
		i++;
	}
	////////////////////
	
	train_asm(nn, &inputs, &outputs);

	printmlp(nn);
	printoutputs(nn, inputs);
	printerrors(nn, inputs, outputs);
	vecfree_asm(&x);      
	matfree_asm(&inputs);
	matfree_asm(&outputs);
	freemlp_asm(nn);
	
	/* // Test exp_asm
	uint64 testsize = 100;
	vec t = linspace_asm(-10., 10., testsize);
	vec y = vecalloc_asm(testsize);
	expvec_asm(&t, &y);
	printvec(y);

	vecfree_asm(&t);
	vecfree_asm(&y);
	*/
	getchar();

	return 0;
}
#endif 