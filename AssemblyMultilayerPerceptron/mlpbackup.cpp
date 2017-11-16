//void test(int input)
//{
//	if (!input) {
//		printf("Error!");
//	}
//}
//
//mlp* makemlp(uint64 numlayers, uint64* layersizes, float learningrate, uint64 numepochs)
//{
//	mlp* result = malloc(sizeof(mlp));
//
//	result->numlayers = numlayers;
//	result->learningrate = learningrate;
//	result->numepochs = numepochs;
//
//	uint64 numweights = numlayers - 1;
//
//	result->layers    = (vec*)malloc(numlayers * sizeof(vec));
//	result->layersErr = (vec*)malloc(numlayers * sizeof(vec));
//
//	result->weights    = (mat*)malloc(numweights * sizeof(mat));
//	result->weightsErr = (mat*)malloc(numweights * sizeof(mat));
//
//	result->biases    = (vec*)malloc(numweights * sizeof(vec));
//	result->biasesErr = (vec*)malloc(numweights * sizeof(vec));
//
//	uint64 l = 0;
//	while (l < numlayers)
//	{
//		result->layers[l] = veczeros(layersizes[l]);
//		result->layersErr[l] = veczeros(layersizes[l]);
//		l++;
//	}
//
//	uint64 w = 0;
//	while (w < numweights)
//	{
//		mat m = matrand(layersizes[w + 1], layersizes[w]);
//		mat m1 = matones(m.rows, m.cols);
//		scalemat(.5f, m1, m1);
//		submat(m, m1, m);
//		result->weights[w] = m;
//		result->weightsErr[w] = matzeros(layersizes[w + 1], layersizes[w]);
//		matfree(m1);
//
//		vec v = vecrand(layersizes[w + 1]);
//		vec v1 = vecones(v.size);
//		scalevec(.5f, v1, v1);
//		subvec(v, v1, v);
//		result->biases[w] = v;
//		result->biasesErr[w] = veczeros(layersizes[w + 1]);
//		vecfree(v1);
//
//		w++;
//	}
//
//	return result;
//}
//void freemlp(mlp* net) 
//{
//	uint64 numweights = net->numlayers - 1;
//
//	uint64 l = 0;
//	while (l < net->numlayers)
//	{
//		vecfree(net->layers[l]);
//		vecfree(net->layersErr[l]);
//		l++;
//	}
//
//	uint64 w = 0;
//	while (w < numweights)
//	{
//		matfree(net->weights[w]);
//		matfree(net->weightsErr[w]);
//		vecfree(net->biases[w]);
//		vecfree(net->biasesErr[w]);
//		w++;
//	}
//
//	free(net->layers);
//	free(net->layersErr);
//	free(net->weights);
//	free(net->weightsErr);
//	free(net->biases);
//	free(net->biasesErr);
//}
//
//void sigmoid_activation(vec v, vec result)
//{
//	uint64 index = 0;
//	uint64 size = min(v.size, result.size);
//	while (index < size)
//	{
//		result.data[index] = 1.f / (1.f + expf( -v.data[index]));
//		index++;
//	}
//}
//
//void transfer(mlp* net, uint64 l, int bActivate)
//{
//	test(net->weights[l - 1].cols == net->layers[l - 1].size);
//	test(net->weights[l - 1].rows == net->layers[l].size);
//	matdotvec(net->weights[l - 1], net->layers[l - 1], net->layers[l]);
//
//	test(net->layers[l].size == net->biases[l - 1].size);
//	addvec(net->layers[l], net->biases[l - 1], net->layers[l]); // Add bias term
//	
//	//  - Apply activation function
//	if (bActivate) { sigmoid_activation(net->layers[l], net->layers[l]); }
//}
//
//void fwd(mlp* net, vec input)
//{
//	// Insert input into first layer
//	veccopy(input, net->layers[0]);
//
//	uint64 L = net->numlayers - 1ULL;
//	// For each layer
//	uint64 l = 1;
//	while (l < net->numlayers)
//	{
//		// Skip activation function for last layer
//		int bActivate = l < L;
//
//		//  - transfer to next layer ( from l-1 to l )
//		transfer(net, l, bActivate);
//		l++;
//	}
//}
//
//void derivative_sigmoid_activation(vec v, vec result)
//{
//	uint64 index = 0;
//	uint64 size = min(v.size, result.size);
//	while (index < size)
//	{
//		float sigma = 1.f / (1.f + expf(-v.data[index]));
//		result.data[index] = sigma * (1.f - sigma);
//		index++;
//	}
//}
//
//void bwd(mlp* net, vec output)
//{
//	uint64 L = net->numlayers - 1ULL;
//	float rate = net->learningrate;
//
//	// calculate output error gradient
//	subvec(net->layers[L], output, net->layersErr[L]); // dE/dx_L = eps_L = x_L - y
//
//	uint64 l = L;
//	while (l > 0) // For all weights L-1, ..., 0
//	{
//		// Go to next backwards layer
//		l--;
//
//		// Propagate backwards
//		// - calc propagation factor 
//		// rho_l = eps_{l+1} * f'(x^{l+1})
//		vec rho = vecalloc(net->layers[l + 1].size);
//		if (l == L - 1) {
//			veccopy(net->layersErr[l + 1], rho);
//		}
//		else {
//			//derivative_sigmoid_activation(net->layers[l + 1], rho); // incorrect, activation applied twice
//			veccopy(net->layers[l], rho);
//			matdotvec(net->weights[l], rho, rho);
//			addvec(net->biases[l], rho, rho);
//			derivative_sigmoid_activation(rho, rho);
//			mulvec(rho, net->layersErr[l + 1], rho);
//		}
//
//		// - calc weight error
//		// dE/dW_l = rho_l x_l^T
//		outervec(rho, net->layers[l], net->weightsErr[l]);
//		scalemat(rate, net->weightsErr[l], net->weightsErr[l]); // Use Err as temp storage
//
//		// - calc bias error
//		// dE/db_l = rho_l
//		veccopy(rho, net->biasesErr[l]);
//		scalevec(rate, net->biasesErr[l], net->biasesErr[l]); // Use Err as temp storage
//
//		// Skip input error calc
//		// - calc next backwards layer error
//		if (l > 0)
//		{
//			// eps_l = W_l^T rho_l
//			matTdotvec(net->weights[l], rho, net->layersErr[l]);
//		}
//
//		// update weights by learning rate
//		subvec(net->biases[l], net->biasesErr[l], net->biases[l]);
//		submat(net->weights[l], net->weightsErr[l], net->weights[l]);
//		// free allocated vectors
//		vecfree(rho);
//	}
//}
//
//void learn(mlp* net, vec input, vec output)
//{
//	fwd(net, input);
//	bwd(net, output);
//	// system("cls");
//	// printmlp(net);
//}
//
//void train(mlp* net, mat inputs, mat outputs)
//{
//	uint64 size = min(inputs.rows, outputs.rows);
//	uint64 e = 0;
//	vec input;
//	vec output;
//	input.size = inputs.cols;
//	output.size = outputs.cols;
//	while (e < net->numepochs)
//	{
//		//printf("Epoch: %I64d\n", e);
//		uint64 i = 0;
//		while (i < size)
//		{
//			int idx = rand() % size;
//			input.data = inputs.data[idx];
//			output.data = outputs.data[idx];
//			learn(net, input, output);
//			//printvec(net->layers[net->numlayers - 1]);
//			//printf(", ");
//			i++;
//		}
//		e++;
//	}
//}
//
//void printoutputs(mlp* net, mat inputs)
//{
//	uint64 size = inputs.rows;
//	vec input;
//	input.size = inputs.cols;
//	printf("outputs = ");
//	printf("np.array([");
//
//	uint64 i = 0;
//	while (i < size)
//	{
//		input.data = inputs.data[i];
//		fwd(net, input);
//		printvec(net->layers[net->numlayers - 1]);
//		printf(",");
//		i++;
//	}
//	printf("]).T[0]");
//}
//
//void printmlp(mlp* net)
//{
//	printf("\nPrinting multilayer perceptron state:\n");
//
//	uint64 L = net->numlayers;
//	uint64 N = net->numepochs;
//	float  R = net->learningrate;
//
//	uint64 numweights =net->numlayers - 1;
//
//	printf("Number of layers: %lld\n", L);
//	printf("Number of epochs: %lld\n", N);
//	printf("Rate of learning: %f\n", R);
//
//	uint64 l = 0;
//	while (l < L)
//	{
//		printf("Layer %lld:\n", l);
//		printvec(net->layers[l]);
//
//		printf("Layer errors %lld:\n", l);
//		printvec(net->layersErr[l]);
//		l++;
//	}
//
//	uint64 w = 0;
//	while (w < numweights)
//	{
//		printf("Weights %lld:\n", w);
//		printmat(net->weights[w]);
//		printf("Weight errors %lld:\n", w);
//		printmat(net->weightsErr[w]);
//
//		printf("Biases %lld:\n", w);
//		printvec(net->biases[w]);
//
//		printf("Bias errors %lld:\n", w);
//		printvec(net->biasesErr[w]);
//		w++;
//	}
//	printf("\n");
//}