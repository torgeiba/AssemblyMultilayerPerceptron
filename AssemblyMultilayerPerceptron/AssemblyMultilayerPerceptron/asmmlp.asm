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
extern vecalloc_asm:proc
extern veccopy_asm:proc
extern vecones_asm:proc
extern scalevec_asm:proc
extern expvec_asm:proc
extern addvec_asm:proc
extern subvec_asm:proc
extern mulvec_asm:proc
extern divvec_asm:proc
extern matdotvec_asm:proc
extern matTdotvec_asm:proc
extern outervec_asm:proc
extern vecfree_asm:proc
extern matfree_asm:proc
extern submat_asm:proc
extern scalemat_asm:proc



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

	push r12 ; L
	push r13 ; l
	push r14 ; mlp
	push r15 ; numlayers

	mov r14, rcx ; r14 = mlp* net

	; // Insert input into first layer

	; veccopy_asm(&input, &net->layers[0]);
	mov rcx, rdx		; first param &input
	mov rdx, [r14 + 8]	; second param &net->layers[0];  	vec* layers; (mlp + 8)
	call veccopy_asm

	; uint64 L = net->numlayers - 1ULL; = &mlp + 0
	mov r15, [r14 + 0] ; net->numlayers
	mov r12, r15 ; L = numlayers - 1
	dec r12

	; // For each layer
	; uint64 l = 1;
	mov r13, 1
	layer_loop: ; while (l < net->numlayers) {
		; 	//  - transfer to next layer ( from l-1 to l )
		mov rcx, r14 ; net
		mov rdx, r13 ; l

		; // Skip activation function for last layer
		xor r8, r8
		mov r9, 1
		cmp r13, r12 ; int bActivate = l < L;
		cmovl r8, r9  ; bActivate
		call transfer_asm ; transfer_asm(&net, l, bActivate);
		
	inc r13; 	l++;
	cmp r13, r15
	jl layer_loop
	; }

	pop r15
	pop r14
	pop r13
	pop r12

	ret
fwd_asm endp

; extern void bwd_asm(mlp* net, vec* output);
; input:
; rcx = &net
; rdx = &output
bwd_asm proc											; struct mlp {
														;  	uint64 numlayers;	+ 0
	push r12 ; L										;  	vec* layers;		+ 8
	push r13 ; l										;  	vec* biases;		+ 16
	push r14 ; &rho										;  	mat* weights;		+ 24
	push r15 ; &net										;  	vec* layersErr;		+ 32
	push rbx ; rate										;
														;  	vec* biasesErr;		+ 40
	mov r15, rcx ; r15 = &net							;  	mat* weightsErr;	+ 48
	;uint64 L = net->numlayers - 1ULL;(+0)				;
	mov r12, [r15 + 0]									;  	float learningrate; + 56 ; note!
	dec r12												;  	uint64 numepochs;   + 60 ; prev was 4 byte float
														;  };
	mov rbx, [r15 + 56] ; float rate = net->learningrate;	(+56)
	;// calculate output error gradient			
												
	mov rcx, ; TODO
	mov rdx, rdx; &output
	mov r8 , ; TODO
	call subvec_asm ; subvec_asm(&net->layers[L], &output, &net->layersErr[L]); // dE/dx_L = eps_L = x_L - y

	; vec rho;
	mov r14, rsp ; r14 = &rho
	sub rsp, 16; move sizeof(vec) = 16 bytes
	
	mov r13, r12 ;uint64 l = L;

	bwd_layer_loop: ;while (l > 0) { // For all weights L-1, ..., 0
		dec r13 ;	l--; // Go to next backwards layer
	;
	;	// Propagate backwards
	;	// - calc propagation factor 
	;	// rho_l = eps_{l+1} * f'(x^{l+1})
		mov rcx, r14; &rho return pointer
		mov rdx, [r15 + 8] ; &net->layers[0]

		mov rdx, [r15 + 8]; &net->layers[0];  	vec* layers;		+ 8
		mov r9, r13; r9 = l
		inc r9     ; r9 = l + 1
		shl r9, 4  ; r9 = l * sizeof(vec); sizeof(vec) = 16 = 1000b
		add rdx, r9; &net->layers[l + 1]
		mov rdx, [rdx + 8]; &net->layers[0]->size
		call vecalloc_asm ;	vec rho = vecalloc_asm(/*hidden return addr  &rho */ net->layers[l + 1].size);
	
		mov r9, r12 ; r9 = L
		dec r9      ; r9 = L - 1
		cmp r13, r9 ; l == L - 1
		jne bwd_not_last_layer ;	if (l == L - 1)
		;{
			mov rcx, [r15 + 32]; &net->layersErr[0];  	vec* layersErr;		+ 32
			mov r9, r13; r9 = l
			inc r9     ; r9 = l + 1
			shl r9, 4  ; r9 = (l + 1) * sizeof(vec); sizeof(vec) = 16 = 1000b
			add rcx, r9; &net->layersErr[l + 1] 

			mov rdx, r14; &rho
			call veccopy_asm ;		veccopy_asm(&net->layersErr[l + 1], &rho);
		;}
		jmp bwd_end_last_layer
		bwd_not_last_layer: ;	else {

			mov rcx, [r15 + 24]; &net->weights[0];  	mat* weights;		+ 24
			mov r9, r13; r9 = l
			shl r9, 3  ; r9 = l* sizeof(mat); sizeof(mat) = 24 = 1100b = 1000b + 100b
			add rcx, r9; &net->weights[l]
			shl r9, 1  ; r9 = l* sizeof(mat); sizeof(mat) = 24 = 1100b = 1000b + 100b
			add rcx, r9; &net->weights[l]

			mov rdx, [r15 + 8]; &net->layers[0];  	vec* layers;		+ 8
			mov r9, r13; r9 = l
			shl r9, 4  ; r9 = l * sizeof(vec); sizeof(vec) = 16 = 1000b
			add rdx, r9; &net->layers[l];

			mov r8 , r14; &rho
			call addvec_asm ;		matdotvec_asm(&net->weights[l], &net->layers[l], &rho);


			mov rcx, [r15 + 16]; &net->biases[0];  	vec* biases;		+ 16
			mov r9, r13; r9 = l
			shl r9, 4  ; r9 = l * sizeof(vec); sizeof(vec) = 16 = 1000b
			add rcx, r9; &net->biases[l] 

			mov rdx, r14; &rho
			mov r8 , r14; &rho
			call addvec_asm ;		addvec_asm(&net->biases[l], &rho, &rho);

			mov rcx, r14; &rho
			mov rdx, r14; &rho
			call derivative_sigmoid_activation_asm;		derivative_sigmoid_activation_asm(&rho, &rho);

			mov rcx, r14; &rho
			mov rdx, ; &net->layersErr[l + 1]
			mov r8 , r14; &rho
			call mulvec_asm ; mulvec_asm(&rho, &net->layersErr[l + 1], &rho);
		bwd_end_last_layer:
	;	}
	;
	;	// - calc weight error
	;	// dE/dW_l = rho_l x_l^T
		mov rcx, r14; &rho
		mov rdx, ; TODO
		mov r8 , ; TODO
		call outervec_asm ;	outervec_asm(&rho, &net->layers[l], &net->weightsErr[l]);
	;	

		mov rcx, rbx; rate
		mov rdx, ; TODO
		mov r8 , ; TODO
		call scalemat_asm ; scalemat_asm(rate, &net->weightsErr[l], &net->weightsErr[l]); // Use Err as temp storage
	;
	;	// - calc bias error
	;	// dE/db_l = rho_l
		mov rcx, r14; &rho
		mov rdx, ; TODO
		call veccopy_asm ;	veccopy_asm(&rho, &net->biasesErr[l]);


		mov rcx, rbx; rate
		mov rdx, ; TODO
		mov r8 , ; TODO
		call scalevec_asm ;	scalevec_asm(rate, &net->biasesErr[l], &net->biasesErr[l]); // Use Err as temp storage
	;
	;	// Skip input error calc
	;	// - calc next backwards layer error
		xor r8, r8 
		cmp r13, r8;	if (l > 0) {
		jng skip_input_error_calc
	;		// eps_l = W_l^T rho_l
			mov rcx, ; TODO
			mov rdx, r14 ; &rho
			mov r8 , ; TODO
			call matTdotvec_asm ;		matTdotvec_asm(&net->weights[l], &rho, &net->layersErr[l]);
		skip_input_error_calc:
	;	}
	;
	;	// update weights by learning rate
		mov rcx, ; TODO &net->biases[l]
		mov rdx, ; TODO &net->biasesErr[l]
		mov r8 , rcx; &net->biases[l]
		call subvec_asm ;	subvec_asm(&net->biases[l], &net->biasesErr[l], &net->biases[l]);

		mov rcx, ; TODO &net->weights[l]
		mov rdx, ; TODO &net->weightsErr[l]
		mov r8 , rcx; &net->weights[l]
		call submat_asm ;	submat_asm(&net->weights[l], &net->weightsErr[l], &net->weights[l]);
	
		;	// free allocated vectors
		mov rcx, r14 ; &rho
		call vecfree_asm ;	vecfree_asm(&rho);

	xor r8, r8
	cmp r13, r8 ;(l > 0)
	jg bwd_layer_loop
	;}

	add rsp, 16; move sizeof(vec) = 16 bytes
	pop r15
	pop r14
	pop r13
	pop r12

	ret
bwd_asm endp

; extern void learn_asm(mlp* net, vec* input, vec* output);
; input rcx, rdx, r8
learn_asm proc

	push r12
	push r13

	mov r12, rcx ; r12 = mlp* net
	mov r13, r8  ; r13 = vec* output

	; call fwd(net, input);
	; params already in place
	call fwd_asm

	; call bwd(net, output);
	mov rcx, r12 ; net
	mov rdx, r13 ; output
	call bwd_asm

	pop r13
	pop r12

	ret
learn_asm endp

; extern void train_asm(mlp* net, mat* inputs, mat* outputs);
; input rcx, rdx, r8
train_asm proc

	push r12 ; &input
	push r13 ; &output
	push r14 ; epoch counter
	push r15 ; traningset counter

	; uint64 size = min(inputs.rows, outputs.rows);
	; uint64 e = 0;
	
	sub rsp, 32 ; alloc input and output vecs
	; vec input;  (16 bytes)
	; vec output; (16 bytes)

	mov r12, rsp ; r9 = &input
	
	mov r13, rsp
	add r13, 16  ; r10 = &output

	; input.size = inputs.cols;
	mov r9, [rdx + 16]
	mov [r12 + 8], r9

	; output.size = outputs.cols;
	mov r9, [r8 + 16]
	mov [r13 + 8], r9

	; r14 = e = net->numepochs
	mov r14, [rcx + 60] ;  	uint64 numepochs;   + 60
	epoch_iterations: ; while (e > 0) {
		mov r15, [rdx + 8]; 	uint64 i = size; inputs.rows
		trainingset_iterations: ; 	while (i > 0)
		
		; TODO: optimize, use rsp offsets for input/output vecs
		; to free up non-volatile regs for input pointers to this function
		; so volatile regs do not have to be pushed and popped

		; prep to call rand_asm
		push rcx
		push rdx
		push r8						; 	{
	
		call rand_asm				; 		int idx = rand() % size;
	
		mov r10, [rcx + 8]
		div r10						;		rdx	= idx
		mov r9, rdx					;		r9  = idx

		pop r8
		pop rdx
		pop rcx
		
		; prep to call learn_asm
		push rcx
		push rdx
		push r8

		mov r10, [rdx]			;		&inputs.data[0]
		mov r11, [r8]			;		&output.data[0]
		shl r9, 3				;		multiply index by pointer size(8) to get offsets

		add r10, r9
		add r11, r9
		
		mov [r12], r10			; 		input.data  = inputs.data[idx];
		mov [r12], r10			; 		output.data = outputs.data[idx];

								; first param mlp* net already in rcx
		mov rdx, r12			; second param input
		mov r8 , r13			; third param output

		call learn_asm			; 		learn_asm(&net, &input, &output);
		
		pop r8
		pop rdx
		pop rcx
		
		dec r15					; 		i--;
		jnz trainingset_iterations
	; 	}
	; 	e--;
	dec r14
	jnz epoch_iterations
	; }

	add rsp, 32 ; dealloc input and output vecs
	pop r15
	pop r14
	pop r13
	pop r12


	ret
train_asm endp

end