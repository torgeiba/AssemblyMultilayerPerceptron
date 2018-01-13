;
; x64 Assembly multilayer perceptron
;
; Include in C projects using asmmlp.h
; Requires asmmth math library amd asmclib bindings for C library helper functions
;

; asmclib
extern exp_asm:proc    ; exponential function from C math.h implementation
extern malloc_asm:proc ; dynamic memory allocation on the heap from the C standard library implementation
extern calloc_asm:proc ; dynamic zero-initialized memory allocation on the heap from the C standard library implementation
extern free_asm:proc   ; freeing dynamically allocated heap memory, as implemented in the C standard library
extern rand_asm:proc   ; random integer generator function from C stdlib

; asmmath lib
extern veccopy_asm:proc
extern vecones_asm:proc
extern scalevec_asm:proc
extern expvec_asm:proc
extern addvec_asm:proc
extern subvec_asm:proc
extern mulvec_asm:proc
extern divvec_asm:proc
extern matdotvec_asm:proc
extern vecfree_asm:proc
extern matfree_asm:proc


.data
.code

; ---- C data types ---- ;
; struct mlp {			 ;
; 	uint64 numlayers;	 ;
; 	vec* layers;		 ;
; 	vec* biases;		 ;
; 	mat* weights;		 ;
; 						 ;
; 	vec* layersErr;		 ;
; 	vec* biasesErr;		 ;
; 	mat* weightsErr;	 ;
; 						 ;
; 	float learningrate;	 ;
; 	uint64 numepochs;	 ;
; };					 ;
;						 ;
; struct vec {			 ;
;	float* data;		 ;
;	uint64 size;		 ;
; };					 ;
; 						 ;
; struct mat {			 ;
;	float** data;		 ;
;	uint64 rows;		 ;
;	uint64 cols;		 ;
; };					 ;
; ---------------------- ;

; Calling convention: Microsoft x64.
;
; Summmary:
;
; The first four integer or pointer parameters are passed in the rcx, rdx, r8, and r9 registers.
; The first four floating-point parameters are passed in the first four SSE registers, xmm0-xmm3.
; The caller reserves shadow space on the stack for arguments passed in registers. The called function can use this space to spill the contents of registers to the stack.
; Any additional arguments are passed on the stack.
; An integer or pointer return value is returned in the rax register, while a floating-point return value is returned in xmm0.
; rax, rcx, rdx, r8-r11 are volatile.
; rbx, rbp, rdi, rsi, r12-r15 are nonvolatile.
;
; For more info, go to:
; Calling convention: https://en.wikipedia.org/wiki/X86_calling_conventions#Microsoft_x64_calling_convention
;                     https://msdn.microsoft.com/en-us/library/9b372w95.aspx
; x64 architecture:   https://docs.microsoft.com/en-us/windows-hardware/drivers/debugger/x64-architecture
; Register usage:     https://msdn.microsoft.com/en-us/library/9z1stfyw.aspx

; extern mlp* makemlp_asm(uint64 numlayers, uint64* layersizes, uint64 numepochs, float learningrate);
; input rcx, rdx, r8, r9
; output rax
makemlp_asm proc
	; mlp* result = malloc(sizeof(mlp));
	; 
	; result->numlayers = numlayers;
	; result->learningrate = learningrate;
	; result->numepochs = numepochs;
	; 
	; uint64 numweights = numlayers - 1;
	; 
	; result->layers = (vec*)malloc(numlayers * sizeof(vec));
	; result->layersErr = (vec*)malloc(numlayers * sizeof(vec));
	; 
	; result->weights = (mat*)malloc(numweights * sizeof(mat));
	; result->weightsErr = (mat*)malloc(numweights * sizeof(mat));
	; 
	; result->biases = (vec*)malloc(numweights * sizeof(vec));
	; result->biasesErr = (vec*)malloc(numweights * sizeof(vec));
	; 
	; uint64 l = 0;
	; while (l < numlayers) {
	; 	result->layers[l] = veczeros_asm(layersizes[l]);
	; 	result->layersErr[l] = veczeros_asm(layersizes[l]);
	; 	l++;
	; }
	; 
	; uint64 w = 0;
	; while (w < numweights) {
	; 	mat m = matrand_asm(layersizes[w + 1], layersizes[w]);
	; 	mat m1 = matones_asm(m.rows, m.cols);
	; 	scalemat_asm(.5f, &m1, &m1);
	; 	submat_asm(&m, &m1, &m);
	; 	result->weights[w] = m;
	; 	result->weightsErr[w] = matzeros_asm(layersizes[w + 1], layersizes[w]);
	; 	matfree_asm(&m1);
	; 
	; 	vec v = vecrand_asm(layersizes[w + 1]);
	; 	vec v1 = vecones_asm(v.size);
	; 	scalevec_asm(.5f, &v1, &v1);
	; 	subvec_asm(&v, &v1, &v);
	; 	result->biases[w] = v;
	; 	result->biasesErr[w] = veczeros_asm(layersizes[w + 1]);
	; 	vecfree_asm(&v1);
	; 
	; 	w++;
	; }
	; return result;
	ret
makemlp_asm endp


; extern void freemlp_asm(mlp* net);
; rcx
; struct mlp {
;  	uint64 numlayers;	+ 0
;  	vec* layers;		+ 8
;  	vec* biases;		+ 16
;  	mat* weights;		+ 24
;  	
;  	vec* layersErr;		+ 32
;  	vec* biasesErr;		+ 40
;  	mat* weightsErr;	+ 48
;
;  	float learningrate; + 56 ; note!
;  	uint64 numepochs;   + 60 ; prev was 4 byte float
;  };
freemlp_asm proc

	push r12
	push r13
	push r14

	mov r12, rcx; r12 = mlp* net

	mov r13, [r12 + 0] 		; uint64 l = net->numlayers;
	dec r13
	layer_loop:								; while (l > 0) {
		mov r14, qword ptr [r12 + 8]		;	r14 = &mlp->layers[0]
		mov rdx, r13
		shl rdx, 4							; rdx = r13 * 16
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->layers[l] ; sizeof(vec) = 16
		call vecfree_asm					;	vecfree_asm(&net->layers[l]);

		mov r14, qword ptr [r12 + 32]		;	r14 = mlp->layersErr
		mov rdx, r13
		shl rdx, 4							; rdx = r13 * 16
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->layersErr[l] ; sizeof(vec) = 16
		call vecfree_asm					;	vecfree_asm(&net->layersErr[l]);

	dec r13									;	l--;
	jnz layer_loop							; }
	
	;uint64 w = numweights = net->numlayers - 1;					
	mov r13, [r12 + 0]
	dec r13
	dec r13									;	w--;										
	weight_loop:							; while (w > 0) {	
		mov r14, qword ptr [r12 + 48]		;	r14 = mlp->weightsErr
		mov rdx, r13
		shl rdx, 2
		sub rdx, r13
		shl rdx, 3							; rdx = r13 * 3 * 8 = rdx * 24
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->weightsErr[w]		 ; sizeof(vec) = 24
		call matfree_asm					;	matfree_asm(&net->weightsErr[w]);	
		
		mov r14, qword ptr [r12 + 40]		;	r14 = mlp->biasErr
		mov rdx, r13
		shl rdx, 4							; rdx = r13 * 16
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->biasesErr[w]		; sizeof(vec) = 16
		call vecfree_asm					;	vecfree_asm(&net->biasesErr[w]);	
		
		mov r14, qword ptr [r12 + 24]		;	r14 = mlp->weights
		mov rdx, r13
		shl rdx, 2
		sub rdx, r13
		shl rdx, 3							; rdx = r13 * 3 * 8 = rdx * 24
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->weights[w]			; sizeof(vec) = 24
		call matfree_asm					;	matfree_asm(&net->weights[w]);
		
		mov r14, qword ptr [r12 + 16]		;	r14 = mlp->bias
		mov rdx, r13
		shl rdx, 4							; rdx = r13 * 16
		lea rcx, qword ptr [r14 + rdx]		;	rcx = &net->biases[w]			; sizeof(vec) = 16
		call vecfree_asm					;	vecfree_asm(&net->biases[w]);
		
	dec r13									;	w--;								
	jnz weight_loop							; }										
			
	mov rcx, qword ptr [r12 + 8]	; vec* layers;		+ 8													
	call free_asm					; free(net->layers);
	
	mov rcx, qword ptr [r12 + 32]	; vec* layersErr;		+ 32							
	call free_asm					; free(net->layersErr);
	
	mov rcx, qword ptr [r12 + 24]	; mat* weights;		+ 24								
	call free_asm					; free(net->weights);
	
	mov rcx, qword ptr [r12 + 48]	; mat* weightsErr;	+ 48				
	call free_asm					; free(net->weightsErr);
	
	mov rcx, qword ptr [r12 + 16]	; vec* biases;		+ 16;					
	call free_asm					; free(net->biases);

	mov rcx, qword ptr [r12 + 40]	; vec* biasesErr;		+ 40
	call free_asm					; free(net->biasesErr);					
											
	pop r14									
	pop r13									
	pop r12									
											
	ret
freemlp_asm endp

; extern vec* sigmoid_activation_asm(vec* v, vec* result);
; input rcx, rdx
sigmoid_activation_asm proc
	push r12
	push r13
	mov r12, rdx ; store pointer to result vec in r12

	;veccopy_asm(&v, &result);
	; uses args to this procedure directly
	; no need to save pointer to v in rcx, as it is not used anymore
	call veccopy_asm
	
	;uint64 size = min(v.size, result.size);
	;vec one = vecones_asm((hidden return parameter), size);
	sub rsp, 16 ;vec one
	mov rcx, rsp
	mov rdx, [r12 + 8] ; set arg to vecones_asm to size of result
	call vecones_asm
	mov r13, rax ; store pointer to one vec in r13

	;scalevec_asm(-1.f, &result, &result);
	mov rcx, 1
	neg rcx
	cvtsi2ss xmm0, rcx; set first arg to scalevec to -1.f
	mov rdx, r12		; set second arg to scalevec to pointer to result
	mov r8, r12			; set third arg to scalevec to pointer to result
	call scalevec_asm

	;expvec_asm(&result, &result);
	mov rcx, r12 ; set first arg to expvec_asm to pointer to result
	mov rdx, r12 ; set second arg to expvec_asm to pointer to result
	call expvec_asm

	;addvec_asm(&result, &one, &result);
	mov rcx, r12 ; set first arg to result vec
	mov rdx, r13  ; set second arg to one vec
	mov r8, r12  ; set third arg to result vec
	call addvec_asm

	; divvec_asm(&one, &result, &result);
	mov rcx, r13 ; set first arg to one vec
	mov rdx, r12  ; set second arg to result vec
	mov r8, r12  ; set third arg to result vec
	call divvec_asm

	; vecfree_asm(&one);
	mov rcx, r13 ; set first arg to pointer to one vec
	call vecfree_asm

	mov rax, r12 ; set pointer to result as return value

	add rsp, 16
	pop r13		; restore non-volatile reg r13
	pop r12		; restore non-volatile reg r12

	ret
sigmoid_activation_asm endp

; extern void derivative_sigmoid_activation_asm(vec* v, vec* result);
; input rcx, rdx
derivative_sigmoid_activation_asm proc

	push r12
	push r13
	mov r12, rdx ; store pointer to result vec in r12
	
	; veccopy_asm(&v, &result);
	; uses args to this procedure directly
	; no need to save pointer to v in rcx, as it is not used anymore
	call veccopy_asm

	; uint64 size = min(v.size, result.size);
	; vec one = vecones_asm((hidden return parameter), size);
	sub rsp, 16 ;vec one
	mov rcx, rsp
	mov rdx, [r12 + 8] ; set arg to vecones_asm to size of result
	call vecones_asm
	mov r13, rax ; store pointer to one vec in r13

	; scalevec_asm(-1.f, &result, &result);
	mov rcx, 1
	neg rcx
	cvtsi2ss xmm0, rcx; set first arg to scalevec to -1.f
	mov rdx, r12		; set second arg to scalevec to pointer to result
	mov r8, r12			; set third arg to scalevec to pointer to result
	call scalevec_asm

	; expvec_asm(&result, &result);
	mov rcx, r12 ; set first arg to expvec_asm to pointer to result
	mov rdx, r12 ; set second arg to expvec_asm to pointer to result
	call expvec_asm

	; addvec_asm(&result, &one, &result);
	mov rcx, r12 ; set first arg to result vec
	mov rdx, r13  ; set second arg to one vec
	mov r8, r12  ; set third arg to result vec
	call addvec_asm

	; divvec_asm(&one, &result, &result);
	mov rcx, r13 ; set first arg to one vec
	mov rdx, r12  ; set second arg to result vec
	mov r8, r12  ; set third arg to result vec
	call divvec_asm

	; subvec_asm(&one, &result, &one);
	mov rcx, r13  ; set first arg to one vec
	mov rdx, r12  ; set second arg to result vec
	mov r8, r13   ; set third arg to result vec
	call subvec_asm

	; mulvec_asm(&one, &result, &result);
	mov rcx, r13  ; set first arg to one vec
	mov rdx, r12  ; set second arg to result vec
	mov r8, r12   ; set third arg to result vec
	call mulvec_asm


	;vecfree_asm(&one);
	mov rax, r12 ; set pointer to result as return value

	add rsp, 16
	pop r13		; restore non-volatile reg r13
	pop r12		; restore non-volatile reg r12

	ret
derivative_sigmoid_activation_asm endp


; internal void transfer(mlp* net, uint64 l, int bActivate)
; input rcx, rdx, r8
transfer_asm proc

	push r12
	push r13
	push r14
	push r15

	mov r12, rcx ; r12 = mlp* net
	mov r13, rdx ; r13 = uint64 l
	mov r14, r8  ; r14 = bool bActivate

	; r15 = l-1
	mov r15, r13;
	dec r15


	; struct mlp {
	;  	uint64 numlayers;	+ 0
	;  	vec* layers;		+ 8
	;  	vec* biases;		+ 16
	;  	mat* weights;		+ 24
	;  	
	;  	vec* layersErr;		+ 32
	;  	vec* biasesErr;		+ 40
	;  	mat* weightsErr;	+ 48
	;
	;  	float learningrate; + 56 ; note!
	;  	uint64 numepochs;   + 60 ; prev was 4 byte float
	;  };


	mov r9,  qword ptr [r12 + 8]  ; &net-layers[0]
	mov r10, qword ptr [r12 + 24] ; &net-weights[0]


	; matdotvec
	mov r9, [r12 + 24]	; r9 = &net-weights[0] (mlp+24)
	mov r10, r15		; r10 = l-1
	shl r10, 2				
	sub r10, r15		; r10 *= 3
	shl r10, 3			; r10 *= 8, r10 = 24
	; r10 = (l-1) * 24 = (l-1) * sizeof(mat) 		  
	lea rcx, [r9 + r10] ; second param &net-weights[l - 1]

	mov r9, [r12 + 8]	; r9 = &net-layers[0] (mlp+8)
	mov r10, r15		; r10 = l-1
	shl r10, 4			; r10 = (l-1) * 16 = (l-1) * sizeof(vec) 		  
	lea rdx, [r9 + r10] ; second param &net-layers[l - 1]


	mov r9, [r12 + 8]	; r9 = &net-layers[0]
	mov r10, r13		; r10 = l
	shl r10, 4			; r10 = l * 16 = l * sizeof(vec) 		  
	lea r8, [r9 + r10] ; first param &net-layers[l]

	call matdotvec_asm		; matdotvec_asm(&net->weights[l - 1], &net->layers[l - 1], &net->layers[l]);
	

	; addvec
	mov r9, [r12 + 8]	; r9 = &net-layers[0]
	mov r10, r13		; r10 = l
	shl r10, 4			; r10 = l * 16 = l * sizeof(vec) 		  
	lea rcx, [r9 + r10] ; first param &net-layers[l]

	mov r9, [r12 + 16]	; r9 = &net-biases[0] (mlp+16)
	mov r10, r15		; r10 = l-1
	shl r10, 4			; r10 = (l-1) * 16 = (l-1) * sizeof(vec) 		  
	lea rdx, [r9 + r10] ; second param &net-biases[l - 1]

	mov r8, rcx			; third param &net-layers[l]
	call addvec_asm		; addvec_asm(&net->layers[l], &net->biases[l - 1], &net->layers[l]); // Add bias term
	
	
	; //  - Apply activation function
	; if (bActivate) { sigmoid_activation_asm(&net->layers[l], &net->layers[l]); }
	mov r8, 0
	cmp r14, r8
	jle skip_activation
		mov r9, [r12 + 8]	; &net-layers[0]
		mov r10, r13		; r10 = l
		shl r10, 4			; r10 = l * 16 = l * sizeof(vec*) 		  
		lea rcx, [r9 + r10] ; first param &net-layers[l]
		mov rdx, rcx		; second param &net-layers[l]
		call sigmoid_activation_asm
	skip_activation:

	pop r15
	pop r14
	pop r13
	pop r12

	ret
transfer_asm endp

; extern void fwd_asm(mlp* net, vec* input);
; input rcx, rdx
fwd_asm proc

	; // Insert input into first layer
	; veccopy_asm(&input, &net->layers[0]);
	; 
	; uint64 L = net->numlayers - 1ULL;
	; // For each layer
	; uint64 l = 1;
	; while (l < net->numlayers)
	; {
	; 	// Skip activation function for last layer
	; 	int bActivate = l < L;
	; 
	; 	//  - transfer to next layer ( from l-1 to l )
	; 	transfer(net, l, bActivate);
	; 	l++;
	; }

	ret
fwd_asm endp

; extern void bwd_asm(mlp* net, vec* output);
; input rcx, rdx
bwd_asm proc
	
	;uint64 L = net->numlayers - 1ULL;
	;float rate = net->learningrate;
	;
	;// calculate output error gradient
	;subvec_asm(&net->layers[L], &output, &net->layersErr[L]); // dE/dx_L = eps_L = x_L - y
	;
	;uint64 l = L;
	;while (l > 0) { // For all weights L-1, ..., 0
	;	l--; // Go to next backwards layer
	;
	;	// Propagate backwards
	;	// - calc propagation factor 
	;	// rho_l = eps_{l+1} * f'(x^{l+1})
	;	vec rho = vecalloc(net->layers[l + 1].size);
	;	if (l == L - 1) {
	;		veccopy_asm(&net->layersErr[l + 1], &rho);
	;	}
	;	else {
	;		matdotvec_asm(&net->weights[l], &net->layers[l], &rho);
	;		addvec_asm(&net->biases[l], &rho, &rho);
	;		derivative_sigmoid_activation(rho, rho);
	;		mulvec_asm(&rho, &net->layersErr[l + 1], &rho);
	;	}
	;
	;	// - calc weight error
	;	// dE/dW_l = rho_l x_l^T
	;	outervec_asm(&rho, &net->layers[l], &net->weightsErr[l]);
	;	scalemat_asm(rate, &net->weightsErr[l], &net->weightsErr[l]); // Use Err as temp storage
	;
	;	// - calc bias error
	;	// dE/db_l = rho_l
	;	veccopy_asm(&rho, &net->biasesErr[l]);
	;	scalevec_asm(rate, &net->biasesErr[l], &net->biasesErr[l]); // Use Err as temp storage
	;
	;	// Skip input error calc
	;	// - calc next backwards layer error
	;	if (l > 0)
	;	{
	;		// eps_l = W_l^T rho_l
	;		matTdotvec_asm(&net->weights[l], &rho, &net->layersErr[l]);
	;	}
	;
	;	// update weights by learning rate
	;	subvec_asm(&net->biases[l], &net->biasesErr[l], &net->biases[l]);
	;	submat_asm(&net->weights[l], &net->weightsErr[l], &net->weights[l]);
	;	// free allocated vectors
	;	vecfree_asm(&rho);
	;}

	ret
bwd_asm endp

; extern void learn_asm(mlp* net, vec* input, vec* output);
; input rcx, rdx, r8
learn_asm proc

	; call fwd(net, input);

	; call bwd(net, output);

	ret
learn_asm endp

; extern void train_asm(mlp* net, mat* inputs, mat* outputs);
; input rcx, rdx, r8
train_asm proc

	; uint64 size = min(inputs.rows, outputs.rows);
	; uint64 e = 0;
	; vec input;
	; vec output;
	; input.size = inputs.cols;
	; output.size = outputs.cols;
	; while (e < net->numepochs)
	; {
	; 	uint64 i = 0;
	; 	while (i < size)
	; 	{
	; 		int idx = rand() % size;
	; 		input.data = inputs.data[idx];
	; 		output.data = outputs.data[idx];
	; 		learn(net, input, output);
	; 		i++;
	; 	}
	; 	e++;
	; }

	ret
train_asm endp

end