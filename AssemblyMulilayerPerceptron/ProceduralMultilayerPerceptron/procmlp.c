#include "procmlp.h"

#include "../ProceduralMathematicsLibrary/procmth.h"
#include "../ProceduralMathematicsLibrary/procmthutils.h"

#include "../AssemblyMathematicsLibrary/asmmth.h"

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

void test(int input)
{
	if (!input) {
		printf("Error!");
	}
}

mlp* makemlp(uint64 numlayers, uint64* layersizes, float learningrate, uint64 numepochs)
{
	mlp* result = malloc(sizeof(mlp));

	result->numlayers = numlayers;
	result->learningrate = learningrate;
	result->numepochs = numepochs;

	uint64 numweights = numlayers - 1;

	result->layers = (vec*)malloc(numlayers * sizeof(vec));
	result->layersErr = (vec*)malloc(numlayers * sizeof(vec));

	result->weights = (mat*)malloc(numweights * sizeof(mat));
	result->weightsErr = (mat*)malloc(numweights * sizeof(mat));

	result->biases = (vec*)malloc(numweights * sizeof(vec));
	result->biasesErr = (vec*)malloc(numweights * sizeof(vec));

	uint64 l = 0;
	while (l < numlayers)
	{
		result->layers[l] = veczeros_asm(layersizes[l]);
		result->layersErr[l] = veczeros_asm(layersizes[l]);
		l++;
	}

	uint64 w = 0;
	while (w < numweights)
	{
		mat m = matrand_asm(layersizes[w + 1], layersizes[w]);
		mat m1 = matones_asm(m.rows, m.cols);
		scalemat_asm(.5f, &m1, &m1);
		submat_asm(&m, &m1, &m);
		result->weights[w] = m;
		result->weightsErr[w] = matzeros_asm(layersizes[w + 1], layersizes[w]);
		matfree_asm(&m1);

		vec v = vecrand_asm(layersizes[w + 1]);
		vec v1 = vecones_asm(v.size);
		scalevec_asm(.5f, &v1, &v1);
		subvec_asm(&v, &v1, &v);
		result->biases[w] = v;
		result->biasesErr[w] = veczeros_asm(layersizes[w + 1]);
		vecfree_asm(&v1);

		w++;
	}

	return result;
}
void freemlp(mlp* net)
{
	uint64 numweights = net->numlayers - 1;

	uint64 l = 0;
	while (l < net->numlayers)
	{
		vecfree_asm(&net->layers[l]);
		vecfree_asm(&net->layersErr[l]);
		l++;
	}

	uint64 w = 0;
	while (w < numweights)
	{
		matfree_asm(&net->weights[w]);
		matfree_asm(&net->weightsErr[w]);
		vecfree_asm(&net->biases[w]);
		vecfree_asm(&net->biasesErr[w]);
		w++;
	}

	free(net->layers);
	free(net->layersErr);
	free(net->weights);
	free(net->weightsErr);
	free(net->biases);
	free(net->biasesErr);
}

void sigmoid_activation(vec v, vec result)
{
	uint64 index = 0;
	uint64 size = min(v.size, result.size);
	while (index < size)
	{
		result.data[index] = 1.f / (1.f + expf(-v.data[index]));
		index++;
	}
}

void transfer(mlp* net, uint64 l, int bActivate)
{
	test(net->weights[l - 1].cols == net->layers[l - 1].size);
	test(net->weights[l - 1].rows == net->layers[l].size);
	//matdotvec_unsafe_asm(&net->weights[l - 1], &net->layers[l - 1], &net->layers[l]);
	matdotvec_asm(&net->weights[l - 1], &net->layers[l - 1], &net->layers[l]);

	test(net->layers[l].size == net->biases[l - 1].size);
	addvec_asm(&net->layers[l], &net->biases[l - 1], &net->layers[l]); // Add bias term

	//  - Apply activation function
	if (bActivate) { sigmoid_activation(net->layers[l], net->layers[l]); }
}

void fwd(mlp* net, vec input)
{
	// Insert input into first layer
	veccopy_asm(&input, &net->layers[0]);

	uint64 L = net->numlayers - 1ULL;
	// For each layer
	uint64 l = 1;
	while (l < net->numlayers)
	{
		// Skip activation function for last layer
		int bActivate = l < L;

		//  - transfer to next layer ( from l-1 to l )
		transfer(net, l, bActivate);
		l++;
	}
}

void derivative_sigmoid_activation(vec v, vec result)
{
	uint64 index = 0;
	uint64 size = min(v.size, result.size);
	while (index < size)
	{
		float sigma = 1.f / (1.f + expf(-v.data[index]));
		result.data[index] = sigma * (1.f - sigma);
		index++;
	}
}

void bwd(mlp* net, vec output)
{
	uint64 L = net->numlayers - 1ULL;
	float rate = net->learningrate;

	// calculate output error gradient
	subvec_asm(&net->layers[L], &output, &net->layersErr[L]); // dE/dx_L = eps_L = x_L - y

	uint64 l = L;
	while (l > 0) // For all weights L-1, ..., 0
	{
		// Go to next backwards layer
		l--;

		// Propagate backwards
		// - calc propagation factor 
		// rho_l = eps_{l+1} * f'(x^{l+1})
		vec rho = vecalloc(net->layers[l + 1].size);
		if (l == L - 1) {
			veccopy_asm(&net->layersErr[l + 1], &rho);
		}
		else {
			//derivative_sigmoid_activation(net->layers[l + 1], rho); // incorrect, activation applied twice
			//veccopy_asm(&net->layers[l], &rho);
			//matdotvec(net->weights[l], rho, rho);
			matdotvec_asm(&net->weights[l], &net->layers[l], &rho);
			addvec_asm(&net->biases[l], &rho, &rho);
			derivative_sigmoid_activation(rho, rho);
			mulvec_asm(&rho, &net->layersErr[l + 1], &rho);
		}

		// - calc weight error
		// dE/dW_l = rho_l x_l^T
		outervec_asm(&rho, &net->layers[l], &net->weightsErr[l]);
		scalemat_asm(rate, &net->weightsErr[l], &net->weightsErr[l]); // Use Err as temp storage

																// - calc bias error
																// dE/db_l = rho_l
		veccopy_asm(&rho, &net->biasesErr[l]);
		scalevec_asm(rate, &net->biasesErr[l], &net->biasesErr[l]); // Use Err as temp storage

															  // Skip input error calc
															  // - calc next backwards layer error
		if (l > 0)
		{
			// eps_l = W_l^T rho_l
			matTdotvec_asm(&net->weights[l], &rho, &net->layersErr[l]);
		}

		// update weights by learning rate
		subvec_asm(&net->biases[l], &net->biasesErr[l], &net->biases[l]);
		submat_asm(&net->weights[l], &net->weightsErr[l], &net->weights[l]);
		// free allocated vectors
		vecfree_asm(&rho);
	}
}

void learn(mlp* net, vec input, vec output)
{
	fwd(net, input);
	bwd(net, output);
	// system("cls");
	// printmlp(net);
}

void train(mlp* net, mat inputs, mat outputs)
{
	uint64 size = min(inputs.rows, outputs.rows);
	uint64 e = 0;
	vec input;
	vec output;
	input.size = inputs.cols;
	output.size = outputs.cols;
	while (e < net->numepochs)
	{
		//printf("Epoch: %I64d\n", e);
		uint64 i = 0;
		while (i < size)
		{
			int idx = rand() % size;
			input.data = inputs.data[idx];
			output.data = outputs.data[idx];
			learn(net, input, output);
			//printvec(net->layers[net->numlayers - 1]);
			//printf(", ");
			i++;
		}
		e++;
	}
}

void printoutputs(mlp* net, mat inputs)
{
	uint64 size = inputs.rows;
	vec input;
	input.size = inputs.cols;
	printf("outputs = ");
	printf("np.array([");

	uint64 i = 0;
	while (i < size)
	{
		input.data = inputs.data[i];
		fwd(net, input);
		printcol(net->layers[net->numlayers - 1]);
		printf(",");
		i++;
	}
	printf("]).T[0]\n");
}

void printerrors(mlp* net, mat inputs, mat outputs) {

	uint64 size = inputs.rows;
	vec input;
	input.size = inputs.cols;

	float TotalError = 0.;
	float MaxError = 0.;

	uint64 i = 0;
	while (i < size)
	{
		input.data = inputs.data[i];
		fwd(net, input);
		float error = fabs(net->layers[net->numlayers - 1].data[0] - outputs.data[i][0]);
		TotalError += error;
		if (i == 0 || error < MaxError) MaxError = error;
		i++;
	}
	printf("Total Error = %f\n", TotalError);
	printf("Max Error = %f\n", MaxError);
}

void printmlp(mlp* net)
{
	printf("\nPrinting multilayer perceptron state:\n");

	uint64 L = net->numlayers;
	uint64 N = net->numepochs;
	float  R = net->learningrate;

	uint64 numweights = net->numlayers - 1;

	printf("Number of layers: %lld\n", L);
	printf("Number of epochs: %lld\n", N);
	printf("Rate of learning: %f\n", R);

	uint64 l = 0;
	while (l < L)
	{
		printf("Layer %lld:\n", l);
		printvec(net->layers[l]);

		printf("Layer errors %lld:\n", l);
		printvec(net->layersErr[l]);
		l++;
	}

	uint64 w = 0;
	while (w < numweights)
	{
		printf("Weights %lld:\n", w);
		printmat(net->weights[w]);
		printf("Weight errors %lld:\n", w);
		printmat(net->weightsErr[w]);

		printf("Biases %lld:\n", w);
		printvec(net->biases[w]);

		printf("Bias errors %lld:\n", w);
		printvec(net->biasesErr[w]);
		w++;
	}
	printf("\n");
}


