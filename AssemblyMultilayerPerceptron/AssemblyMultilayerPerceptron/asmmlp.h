#pragma once

#include "../AssemblyMathematicsLibrary/asmmth.h"

/*
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
*/
extern mlp* makemlp_asm(uint64 numlayers, uint64* layersizes, uint64 numepochs, float learningrate);
extern void freemlp_asm(mlp* net);
extern vec* sigmoid_activation_asm(vec* v, vec* result);
extern void derivative_sigmoid_activation_asm(vec* v, vec* result);
extern void fwd_asm(mlp* net, vec* input);
extern void bwd_asm(mlp* net, vec* output);
extern void learn_asm(mlp* net, vec* input, vec* output);
extern void train_asm(mlp* net, mat* inputs, mat* outputs);
// void printoutputs(mlp* net, mat inputs);
// void printerrors(mlp* net, mat inputs, mat outputs);
// void printmlp(mlp* net);
