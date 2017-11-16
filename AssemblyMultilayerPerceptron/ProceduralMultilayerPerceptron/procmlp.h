#pragma once

#include "../ProceduralMathematicsLibrary/procmth.h"
#include "../ProceduralMathematicsLibrary/procmthutils.h"

#include "../AssemblyMathematicsLibrary/asmmth.h"

struct mlp {
	uint64 numlayers;
	vec* layers;
	vec* biases;
	mat* weights;
	
	vec* layersErr;
	vec* biasesErr;
	mat* weightsErr;

	float learningrate;
	uint64 numepochs;
};
typedef struct mlp mlp;

mlp* makemlp(uint64 numlayers, uint64* layersizes, float learningrate, uint64 numepochs);
void freemlp(mlp* net);
void sigmoid_activation(vec v, vec result);
void derivative_sigmoid_activation(vec v, vec result);
void fwd(mlp* net, vec input);
void bwd(mlp* net, vec output);
void learn(mlp* net, vec input, vec output);
void train(mlp* net, mat inputs, mat outputs);
void printoutputs(mlp* net, mat inputs);
void printerrors(mlp* net, mat inputs, mat outputs);
void printmlp(mlp* net);