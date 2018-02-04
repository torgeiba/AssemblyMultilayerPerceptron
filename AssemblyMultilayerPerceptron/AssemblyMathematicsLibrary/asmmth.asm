;
; x64 Assembly math routines
;
; Include in C projects using asmmth.h
;

; ---- C data types ---- ;
; struct vec {			 ;
;	float* data;		 ;
;	uint64 size;		 ;
; };					 ;
; typedef struct vec vec;;
; 						 ;
; struct mat {			 ;
;	float** data;		 ;
;	uint64 rows;		 ;
;	uint64 cols;		 ;
; };					 ;
; typedef struct mat mat;;
; -----------------------;

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

; Internal calling conventions:
; The datatypes vec and mat are passed by memory address.
; This means that they will be passed using registers rcx, rdx, and r8.
; No function should more than 4 parameters, so that they are all passed by register.
; External programs call wrapper functions which in turn calls helper procedures and the core procedure.
; Helper procedures contain parameter reorganization logic which is common for multiple core procedures.
; Core procedures can take arguments from registers RAX, RCX, RDX, R8, R9, R10, and R11 for addresses and integers,
; where RCX, RDX, R8 and R9, are used for the first 4 parameters, to match the Microsoft x64 calling convention, then R10, R11, and RAX
; and XMM0, XMM1, XMM2, XMM3 for floatingpoint arguments.
; The internal parameter passing registers are also all considered volatile when calling the core procedures.
; In procedures with multiple vec or mat parameters, some of the dimensional sizes are redundant.
; The redundant parameters are removed in the helper procedures.
; Return values are handled as per the Microsoft x64 convention.

extern exp_asm:proc    ; exponential function from C math.h implementation
extern malloc_asm:proc ; dynamic memory allocation on the heap from the C standard library implementation
extern calloc_asm:proc ; dynamic zero-initialized memory allocation on the heap from the C standard library implementation
extern free_asm:proc   ; freeing dynamically allocated heap memory, as implemented in the C standard library
;extern rand_asm:proc   ; random integer generator function from C stdlib

.data
asm_rnd_state qword 7777

.code

; C standard library function replacements:

; void srand_asm(uint64 seed); seeds rand_asm() random number generator by setting global variable
srand_asm proc
	mov qword ptr [asm_rnd_state], rcx
	ret
srand_asm endp

; int rand_asm(); produces 15 bit random int, corresponding to C-stdlib rand()
rand_asm proc

	mov rcx, qword ptr [asm_rnd_state]
	mov rdx, rcx
	shl rcx, 13
	xor rdx, rcx
	shr rdx, 17
	xor rcx, rdx
	shl rcx, 5
	xor rdx, rcx

	mov qword ptr [asm_rnd_state], rdx

	mov rax, 7FFFh
	and rax, rdx

	ret
rand_asm endp


; End C standard library function replacements.





; Helper procedures

; vec parameters

; For when the procedure takes a single vec* argument 'u'.
; input:
; rcx = address of u
pass_vec_parameter proc
	mov rdx, qword ptr [rcx + 8] ; Load size of u into rdx
	mov rcx, qword ptr [rcx]     ; load data pointer of u into rcx
	ret
pass_vec_parameter endp

; For when the procedure takes two vec* arguments 'u' and 'v' with the same size.
; input:
; rcx = address of u
; rdx = address of v
pass_two_vec_parameters proc
	mov r8 , qword ptr [rcx + 8] ; Load size of u into r8

	mov rcx, qword ptr [rcx]     ; load data pointer of u into rcx
	mov rdx, qword ptr [rdx]     ; load data pointer of v into rdx
	ret
pass_two_vec_parameters endp

; For when the procedure takes three vec* arguments 'u', 'v', 'w' with the same size.
; input:
; rcx = address of u
; rdx = address of v
; r8  = address of w
pass_three_vec_parameters proc
	mov r9 , qword ptr [rcx + 8] ; Load size of u into R9

	mov rcx, qword ptr [rcx]     ; load data pointer of u into rcx
	mov rdx, qword ptr [rdx]     ; load data pointer of v into rdx
	mov r8 , qword ptr [r8 ]     ; load data pointer of w into r8
	ret
pass_three_vec_parameters endp


; mat parameters

; For when the procedure takes a single mat* argument 'm'.
; input:
; rcx = address of m
pass_mat_parameter proc
	mov r8 , qword ptr [rcx + 10h] ; Load col count of m into RCX
	mov rdx, qword ptr [rcx +  8]  ; Load row count of m into RCX
	mov rcx, qword ptr [rcx]       ; load data pointer of m into rdx
	ret
pass_mat_parameter endp

; For when the procedure takes two mat* arguments 'm' and 'n' with the same size.
; input:
; rcx = address of m
; rdx  = address of n
pass_two_mat_parameters proc
	mov r9, qword ptr [rcx + 10h] ; load col count of m into r9
	mov r8, qword ptr [rcx +  8]  ; Load row count of m into r8

	mov rcx, qword ptr [rcx]     ; load data pointer of m into rcx
	mov rdx, qword ptr [rdx]     ; load data pointer of n into rdx
	ret
pass_two_mat_parameters endp

; For when the procedure takes three mat* arguments 'm', 'n', 'o' with the same size.
; input:
; rcx = address of m
; rdx = address of n
; r8  = address of o
pass_three_mat_parameters proc
	mov r10, qword ptr [rcx + 10h]  ; Load col count of m into r10
	mov r9 , qword ptr [rcx +  8]   ; Load row count of m into r9

	mov rcx, qword ptr [rcx]       ; load data pointer of m into rcx
	mov rdx, qword ptr [rdx]       ; load data pointer of n into rdx
	mov r8 , qword ptr [r8 ]       ; load data pointer of o into r8
	ret
pass_three_mat_parameters endp

; For when the procedure takes three mat* arguments 'm', 'n', 'o' with matrix product compatible sizes.
; m.rows = o.rows = i
; m.cols = n.rows = j
; n.cols = o.cols = k
; (i, j) * (j, k) = (i, k)
; input:
; rcx = address of m 
; rdx = address of n
; r8  = address of o
pass_to_matrix_product proc
	mov r11, qword ptr [rcx + 10h]   ; Load j (col count of m) into R11
	mov r10, qword ptr [rdx + 10h]   ; Load k (col count of n) into r10
	mov r9 , qword ptr [r8  +  8]    ; Load i (row count of o) into r9

	mov rcx, qword ptr [rcx]       ; load data pointer of m into rcx
	mov rdx, qword ptr [rdx]       ; load data pointer of n into rdx
	mov r8 , qword ptr [r8 ]       ; load data pointer of o into r8
	ret
pass_to_matrix_product endp

; Converts from arguments passed by pass_two_mat_parameters
; to arguments passed by pass_two_vector_parameters
; For use with elementwise operations on matrices to reuse vector code.
; Assumes mat rows 
convert_two_mat_to_vec_parameters proc
	mov r10, rdx ; save rdx, as it will be overwritten
	mov rax, r8  ; load rows as an operand for mul
	mul r9	     ; execute mul on rows and cols
	mov rdx, r10 ; restore rdx
	mov rcx, [rcx] ; dereference double (mat) data pointer to single (vec) data pointer
	mov rdx, [rdx] ; dereference double (mat) data pointer to single (vec) data pointer
	mov r8,  rax ; use result of row * col mul as size parameter of veccpy
	ret
convert_two_mat_to_vec_parameters endp


; Converts from arguments passed by pass_three_mat_parameters
; to arguments passed by pass_three_vector_parameters
; For use with elementwise operations on matrices to reuse vector code.
; Assumes mat rows 
; Input : r10 = m.cols, r9  = m.rows, rcx = m.data**, rdx = n.data**, r8  = o.data**
; Output: rcx = u.data*, rdx = v.data* , r8  = w.data* , r9  = u.size
convert_three_mat_to_vec_parameters proc
	mov r11, rdx ; save rdx, as it will be overwritten
	mov rax, r9  ; load rows as an operand for mul
	mul r10	     ; execute mul on rows and cols
	mov rdx, r11 ; restore rdx
	mov rcx, [rcx] ; dereference double (mat) data pointer to single (vec) data pointer
	mov rdx, [rdx] ; dereference double (mat) data pointer to single (vec) data pointer
	mov r8 , [r8]  ; dereference double (mat) data pointer to single (vec) data pointer
	mov r9,  rax   ; use result of row * col mul as vec size parameter
	ret
convert_three_mat_to_vec_parameters endp

; Direct memory to memory copy for contiguous data 
; uses move string instruction (movsb, movsw , movsd)
; uses rep prefix (repeat by counting down RCX)
; Source Index address in RSI register
; Desitantion Index address in RDI register
; RSI and EDI are callee save
; Parameters:
; Count: ECX
; Source: RDC
; Destination: R8
datacpy proc
	; Allocate on stack before pushing regs to avoid poping other stack data instead of saved register values
	; allocation is left to the caller.
	push rsi
	push rdi

	mov rsi, rdx
	mov rdi, r8
	rep movsd

	pop rdi
	pop rsi
	ret 
datacpy endp


; Wrapper and Core procedures

; vec vecalloc_asm(/*return pointer (hidden),*/ uint64 size);
; (hidden parameter) return address: rcx = vec
; rdx = size
vecalloc_asm proc
	mov qword ptr [rcx + 8], rdx
	
	push rcx

	imul rcx, rdx, 4 ; move malloc parameter to first parameter register (4 bytes per float)

	sub rsp, 20h
	call malloc_asm
	add rsp, 20h

	pop rcx

	mov qword ptr [rcx], rax ; store vector data pointer using return value from malloc
	mov rax, rcx ; store pointer to returned struct in
	ret
vecalloc_asm endp

; void vecfree_asm(vec* u);
; rcx = vec* u
; Note: this does not set data pointer to null, which could be desirable
vecfree_asm proc
	mov qword ptr [rcx + 8], 0 ; set size to zero
	mov rcx, qword ptr[rcx] ; get pointer u.data*, whose memory we will free
	sub rsp, 20h
	call free_asm
	add rsp, 20h
	ret
vecfree_asm endp


; vec veczeros_asm(/*return pointer (hidden),*/ uint64 size);
; rcx = size
veczeros_asm proc
	mov qword ptr [rcx + 8], rdx
	
	push rcx

	mov rcx, rdx ; load move count / vector size as first argument
	mov rdx, 4   ; move bytesize / elementsize to 4 bytes per float as second argument

	sub rsp, 20h
	call calloc_asm
	add rsp, 20h

	pop rcx

	mov qword ptr [rcx], rax ; store vector data pointer using return value from malloc
	mov rax, rcx ; store pointer to returned struct in rax
	ret
veczeros_asm endp

; vec vecones_asm(/*return pointer (hidden),*/ uint64 size);
vecones_asm proc
	mov qword ptr [rcx + 8], rdx
	
	push rcx

	mov rcx, rdx ; load move count / vector size as first argument
	mov rdx, 4   ; move bytesize / elementsize to 4 bytes per float as second argument

	sub rsp, 20h
	call calloc_asm
	add rsp, 20h

	pop rcx
	mov qword ptr [rcx], rax ; store vector data pointer using return value from malloc
	mov rax, rcx ; store pointer to returned struct in

	mov r8, qword ptr [rcx+8] ; size in r8
	mov r9, qword ptr [rcx]   ; store data ptr in r9
	
	onesloop:
		mov dword ptr [r9], 3f800000h ; Should be 32-bit single presicion floatingpoint 1.0f
		add r9, 4	
		dec r8
	jnz onesloop 

	ret
vecones_asm endp

; vec* vecrand_asm(/*return pointer (hidden),*/ uint64 size);
; rdx = size
; use stdlib's rand(), returning an int, and convert to float in range [0.f, 1.f)
; can use conversion instructions cvtsi2sd / cvtsi2ss and cvtsd2ss to perform this conversion
; Needs RAND_MAX to scale with. RAND_MAX is a simple define: #define RAND_MAX 0x7fff
vecrand_asm proc

	mov qword ptr [rcx + 8], rdx
	
	push rcx

	mov rcx, rdx ; load move count / vector size as first argument
	mov rdx, 4   ; move bytesize / elementsize to 4 bytes per float as second argument

	sub rsp, 20h
	call calloc_asm
	add rsp, 20h

	pop rcx
	mov qword ptr [rcx], rax ; store vector data pointer using return value from malloc
	
	push rcx ; store pointer to returned struct on stack, to later pop to rax

	
	push r12 ; store in non-volatile regs to prevent interference from rand()
	push r13
	push r14

	mov r12, qword ptr [rcx+8] ; size in r12
	mov r13, qword ptr [rcx]   ; store data ptr in r13
	mov eax, 7fffh
	cvtsi2ss xmm0,  eax; (float) RAND_MAX
	movd r14, xmm0 ;store in non-volatile register
	randloop:
		call rand_asm ; load random int in eax
		cvtsi2ss xmm0, eax ; convert random int to random float
		movd xmm1, r14    ; retrieve (float)RAND_MAX 
		divss xmm0, xmm1   ; scale to interval [0.f, 1.f)
		movss dword ptr [r13], xmm0 ;
		add r13, 4	
		dec r12
	jnz randloop 

	pop r14
	pop r13
	pop r12

	pop rax ; get stored pointer to vec from stack
	
	ret
vecrand_asm endp

; vec* veccopy_asm(vec* fromvec, *vec tovec);
veccopy_asm proc
	mov rax, rdx ; store pointer to result as return value
	call pass_two_vec_parameters ; (fromvec.data* rcx, tovec.data* rdx, size r8)
	call veccpy
	ret
veccopy_asm endp

; (fromvec.data* rcx, tovec.data* rdx, size r8)
veccpy proc
	cpyloop:
		movss xmm0, dword ptr [rcx]
		movss dword ptr [rdx], xmm0

		add rcx, 4 ; step to next dword
		add rdx, 4
		dec r8
		jnz cpyloop
	ret
veccpy endp

; mat* matcopy_asm(mat* frommat, mat* tomat);
matcopy_asm proc
	mov r11, rdx ; store pointer to dest for return
	call pass_two_mat_parameters ; r9 = frommat.col, r8 = frommat.row , rcx  = frommat.data**, rdx = tomat.data**
	call convert_two_mat_to_vec_parameters ; The mat version of cpy assumes that matrix data is stored contiguously
	call veccpy
	mov rax, r11
	ret
matcopy_asm endp

; vec linspace_asm(/*return pointer (hidden),*/ float start, float end, uint64 steps);
linspace_asm proc
	mov qword ptr [rcx + 8], r9
	
	cmp r9, 0
	jg nonzero_size
		mov qword ptr [rcx], 0
		ret
	nonzero_size:

	push rcx

	mov rcx, r9  ; load move count / vector size as first argument
	mov rdx, 4   ; move bytesize / elementsize to 4 bytes per float as second argument

	push r13 ; 
	push r14
	movd r13, xmm1 ; save volatile value (start)
	movd r14, xmm2 ; save volatile value (end)

	sub rsp, 20h
	call calloc_asm
	add rsp, 20h

	
	movd xmm1, r13 ; retrieve volatile value ( start )
	movd xmm2, r14 ; retrieve volatile value ( end )
	pop r14
	pop r13 ; 

	pop rcx

	subss xmm2, xmm1
	mov r11, qword ptr [rcx + 8]
	dec r11
	cvtsi2ss xmm0, r11
	divss xmm2, xmm0

	mov qword ptr [rcx], rax ; store vector data pointer using return value from calloc
	mov rax, rcx ; store pointer to returned struct in

	mov r8, qword ptr [rcx+8] ; size in r8
	mov r9, qword ptr [rcx]   ; store data ptr in r9

	linspaceloop:
		movss dword ptr [r9], xmm1 ; Should be 32-bit single presicion floatingpoint 1.0f
		add r9, 4
		addss xmm1, xmm2
		dec r8
	jnz linspaceloop 
	ret
linspace_asm endp

; vec* addvec_asm(vec* u, vec* v, vec* result);
addvec_asm proc
	mov rax, r8 ; store pointer to result for return
	call pass_three_vec_parameters ; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
	call addvec_impl
	ret
addvec_asm endp

; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
addvec_impl proc
	addloop:
		movss xmm0, dword ptr [rcx] ; load u.data[i] to reg
		movss xmm1, dword ptr [rdx] ; load v.data[i] to reg
		addss xmm0, xmm1			; add u.data[i] and v.data[i]
		movss dword ptr [r8], xmm0  ; store result to w.data[i]

		add rcx, 4 ; step to next dword
		add rdx, 4
		add r8 , 4
		dec r9
		jnz addloop
	ret
addvec_impl endp

; vec* subvec_asm(vec* u, vec* v, vec* result);
subvec_asm proc
	mov rax, r8 ; store pointer to result in return value
	call pass_three_vec_parameters ; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
	call subvec_impl
	ret
subvec_asm endp

; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
subvec_impl proc
	subloop:
		movss xmm0, dword ptr [rcx] ; load u.data[i] to reg
		movss xmm1, dword ptr [rdx] ; load v.data[i] to reg
		subss xmm0, xmm1			; sub u.data[i] and v.data[i]
		movss dword ptr [r8], xmm0  ; store result to w.data[i]

		add rcx, 4 ; step to next dword
		add rdx, 4
		add r8 , 4
		dec r9
		jnz subloop
	ret
subvec_impl endp

; vec* mulvec_asm(vec* u, vec* v, vec* result);
mulvec_asm proc
	mov rax, r8 ; store pointer to result in return value
	call pass_three_vec_parameters ; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
	call mulvec_impl
	ret
mulvec_asm endp

; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
mulvec_impl proc
	mulloop:
		movss xmm0, dword ptr [rcx] ; load u.data[i] to reg
		movss xmm1, dword ptr [rdx] ; load v.data[i] to reg
		mulss xmm0, xmm1			; mul u.data[i] and v.data[i]
		movss dword ptr [r8], xmm0  ; store result to w.data[i]

		add rcx, 4 ; step to next dword
		add rdx, 4
		add r8 , 4
		dec r9
		jnz mulloop
	ret
mulvec_impl endp

; vec* divvec_asm(vec* u, vec* v, vec* result);
divvec_asm proc
	mov rax, r8 ; store pointer to result in return value
	call pass_three_vec_parameters ; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
	call divvec_impl
	ret
divvec_asm endp

; rcx = u.data*, rdx = v.data*, r8  = w.data*, r9 = u.size
divvec_impl proc
	divloop:
		movss xmm0, dword ptr [rcx] ; load u.data[i] to reg
		movss xmm1, dword ptr [rdx] ; load v.data[i] to reg
		divss xmm0, xmm1			; div u.data[i] and v.data[i]
		movss dword ptr [r8], xmm0  ; store result to w.data[i]

		add rcx, 4 ; step to next dword
		add rdx, 4
		add r8 , 4
		dec r9
		jnz divloop
	ret
divvec_impl endp

; vec* scalevec_asm(float a, vec* v, vec* result);
scalevec_asm proc
	; float a in xmm0
	mov rax, r8 ; store pointer to result in return value
	mov r9 , qword ptr [rdx + 8] ; Load size of u into R9
	call sclvec
	ret		
scalevec_asm endp

sclvec proc
	mov rdx, qword ptr [rdx]     ; load data pointer of u into rdx
	mov r8 , qword ptr [r8 ]     ; load data pointer of v into r8
	scalevecloop:
		movss xmm1, dword ptr [rdx] ; load u.data[i] to reg
		mulss xmm1, xmm0			; div u.data[i] and v.data[i]
		movss dword ptr [r8], xmm1  ; store result to w.data[i]

		add rdx, 4
		add r8 , 4
		dec r9
		jnz scalevecloop
	ret	
sclvec endp



; vec* expvec_asm(vec* u, vec* result)
; rcx = vec* u, rdx = vec* result
expvec_asm proc
	push rdx
	
	; Use callee save registers to avoid having to push and pop
	; registers before and after each call to exp
	push r12
	push r13
	push r14

	mov r14, qword ptr [rcx+8]  ; r14  = u.size
	mov r13, qword ptr [rcx]	; r13 = u.data*
	mov r12, qword ptr [rdx]	; r12 = result.data*

	exploop:
		movss xmm0, dword ptr [r13] ; load u[i] as float argument (xmm0) to exp
		sub rsp, 20h ; reserve shadowspace
		call exp_asm
		add rsp, 20h ; deallocate shadowspace

		movss dword ptr [r12], xmm0 ; store result

		add r12, 4 ; result.data* += 1
		add r13, 4 ; u.data* += 1
		dec r14
	jnz exploop

	pop r14
	pop r13
	pop r12
	pop rax ; get result pointer as return value from first rdx push
	ret
expvec_asm endp


; float vecdotvec_asm(vec* u, vec* v);
; rcx = vec* u, rdx = vec* v
vecdotvec_asm proc
	call pass_two_vec_parameters
	call vdotv
	ret ; returns with result in xmm0
vecdotvec_asm endp

; rcx = u.data*, rdx = v.data*, r8 = u.size
vdotv proc
	pxor xmm0, xmm0 ; Initialize result variable to 0
	vecdotloop:
		movss xmm1, dword ptr [rcx] ; load u.data[i]
		movss xmm2, dword ptr [rdx] ; load v.data[i]
		mulss xmm1, xmm2  ; store product in xmm1
		addss xmm0, xmm1  ; add result to result

		add rcx, 4
		add rdx, 4
		dec r8 ; update loop counter
	jnz vecdotloop
	ret
vdotv endp


; vec* matdotvec_asm(mat* m, vec* v, vec* result);
; rcx = mat* m, rdx = vec* v, r8 = vec* result
; m.rows = result.size
; m.cols = v.size
matdotvec_asm proc
	push r8 ; Store pointer to result to stack
	push rbx ; save rbx to free up for use
	push rbp     ; store stack frame (base) pointer of calling procedure
	mov rbp, rsp ; set stack frame base


	; Store info in regs not touched by vdotv ( rax, rbx, r9, r10, r11 )
	mov rbx, qword ptr [rcx + 16]	; load v.size = m.cols into rbx
	mov r11, qword ptr [rcx + 8]	; load result.size = m.rows into r11 to use as row counter
	mov r10, qword ptr [r8]			; load result.data* into r10
	mov r9 , qword ptr [rdx]		; load v.data* into r9
	mov rax, qword ptr [rcx]		; load m.data**

	; 'allocate' stack data for result.size dwords ( 4 * result.size bytes )
	imul rcx, rbx, 4
	sub rsp, rcx

	; data copy parameters
	mov rcx, rbx		; load size argument for datacpy
	mov rdx, r9			; load source address for datacpy
	mov r8 , rsp		; load destination address for datacpy
	mov r9 , r8			; replace source data for vdotv with stack stored data to enable direct writing to result data (destination)

	CLD ; Clear direction flag to make sure we copy from lowest address to highest address. See https://en.wikipedia.org/wiki/Direction_flag

	call datacpy ; copy data of v to stack

	rowloop:
		mov rcx, qword ptr [rax]; load data of next row of m as u.data*
		mov rdx, r9			; store pointer to v.data* in rdx
		mov r8 , rbx		; store size of 'u' in r8

		call vdotv
		movss dword ptr [r10], xmm0 ; Store computed element of result

		add r10, 4 ; move to next result index
		add rax, 8 ; move to next row. (8-byte / 64 bit) pointer to float*
		dec r11    ; update row counter
	jnz rowloop

	mov rsp, rbp
	pop rbp ; restore stack base pointer
	pop rbx ; restore rbx
	pop rax ; get pointer to result from stack
	ret
matdotvec_asm endp


; vec* matdotvec_asm(mat* m, vec* v, vec* result);
; rcx = mat* m, rdx = vec* v, r8 = vec* result
; m.rows = result.size
; m.cols = v.size
matdotvec_unsafe_asm proc
	push r8 ; Store pointer to result to stack
	push rbx ; save rbx to free up for use

	; Store info in regs not touched by vdotv ( rax, rbx, r9, r10, r11 )
	mov rbx, [rcx + 16] ; load u.size = m.cols into rbx
	mov r11, [rcx + 8]	; load result.size = m.rows into r11 to use as row counter
	mov r10, [r8]		; load result.data* into r10
	mov r9 , [rdx]		; load v.data* into r9
	mov rax, [rcx]		; load m.data**

	rowloop:
		mov rcx, [rax]		; load data of next row of m as u.data*
		mov rdx, r9			; store pointer to v.data* in rdx
		mov r8 , rbx		; store size of 'u' in r8

		call vdotv
		movss dword ptr [r10], xmm0 ; Store computed element of result

		add r10, 4 ; move to next result index
		add rax, 8 ; move to next row. (8-byte / 64 bit) pointer to float*
		dec r11    ; update row counter
	jnz rowloop

	pop rbx ; restore rbx
	pop rax ; get pointer to result from stack
	ret
matdotvec_unsafe_asm endp

; vec* matTdotvec_asm(mat* m, vec* v, vec* result);
matTdotvec_asm proc
	push r8  ; Store pointer to result to stack
	push rbx ; save rbx to free up for use
	push r12 ; save rbx to free up for use

	push rbp     ; store stack frame (base) pointer of calling procedure
	mov rbp, rsp ; set stack frame base

	; Const (r10, r11, r12)
	mov r12, [rcx + 16] ; r12 <- u.size = m.cols      ; const
	mov r11, [r8]		; r11 <- result.data*		  ; const
	mov r10, [rdx]		; r10 <- v.data*			  ; const

	mov rax, [rcx]		; rax <- m.data**  ; Variable (rax, rbx, rcx, rdx, r8)
										   ; rax = m row pointer 
										   ; rbx = m col pointer
										   ; rcx = v element pointer
										   ; rdx = result element pointer
										   ; r8  = i 
										   ; r9  = j

	;'allocate' stack data for result.size dwords ( 4 * result.size bytes )
	
	imul r8, [rdx + 8], 4 ; [rdx + 8] = v.size
	sub rsp, r8

	mov r9, rcx ; store rcx for after datacpy

	; data copy parameters
	mov rcx, [rdx + 8]	; load size argument for datacpy
	mov rdx, [rdx]		; load source address for datacpy
	mov r8 , rsp		; load destination address for datacpy

	CLD ; Clear direction flag to make sure we copy from lowest address to highest address. See https://en.wikipedia.org/wiki/Direction_flag

	call datacpy ; copy data of v to stack

	mov rcx, r9 ; restore rcx from before datacpy

	mov r8, [rcx + 8] ; initialize i = m.rows

	; Set result data to zero
	mov rcx, [rcx + 16] ; load result.size
	mov r9, r11      ; load result.data*
	zeroloop:
		mov dword ptr [r9], 0
		add r9, 4
		dec rcx
	jnz zeroloop

	; loop over rows
	mov rcx, rsp     ; initialize v element pointer
	iloop:
		mov rdx, r11     ; reset result element pointer
		movss xmm0, dword ptr [rcx] ; load v.data[i] into xmm1
		mov rbx, [rax] ; load pointer to u.data* = m.data[i] into rbx
		mov r9, r12    ; load u.size as counter
		jloop:
			movss xmm1, dword ptr [rbx]	; load m.data[i][j] into xmm1
			movss xmm2, dword ptr [rdx] ; load result[i] into xmm2
			mulss xmm1, xmm0			; compute m.data[i][j] * v.data[i]
			addss xmm1, xmm2			; add product to result.data[j]
			movss dword ptr [rdx], xmm1 ; store in result[i]

			add rdx, 4  ; go to next element of result
			add rbx, 4  ; go to next element of u
			dec r9		; decrement col counter
		jnz jloop

		add rcx, 4 ; go to next element of v
		add rax, 8 ; go to next row in m
		dec r8; i--
	jnz iloop

	mov rsp, rbp ; clean up temp stack allocation
	pop rbp ; restore stack frame base
	pop r12
	pop rbx ; restore rbx
	pop rax ; get pointer to result from stack
	ret
matTdotvec_asm endp

; vec* matTdotvec_asm(mat* m, vec* v, vec* result);
; Can not have v and result data overlapping
matTdotvec_unsafe_asm proc
	push r8  ; Store pointer to result to stack
	push rbx ; save rbx to free up for use
	push r12 ; save rbx to free up for use

	; Const (r10, r11, r12)
	mov r12, [rcx + 16] ; r12 <- u.size = m.cols      ; const
	mov r11, [r8]		; r11 <- result.data*		  ; const
	mov r10, [rdx]		; r10 <- v.data*			  ; const

	mov rax, [rcx]		; rax <- m.data**  ; Variable (rax, rbx, rcx, rdx, r8)
										   ; rax = m row pointer 
										   ; rbx = m col pointer
										   ; rcx = v element pointer
										   ; rdx = result element pointer
										   ; r8  = i 
										   ; r9  = j
	mov r8, [rcx + 8] ; initialize i = m.rows

	; Set result data to zero
	mov rcx, [rcx + 16] ; load result.size
	mov r9, r11      ; load result.data*
	zeroloop:
		mov dword ptr [r9], 0
		add r9, 4
		dec rcx
	jnz zeroloop
	
	; loop over rows
	mov rcx, r10     ; initialize v element pointer
	iloop:
		mov rdx, r11     ; reset result element pointer
		movss xmm0, dword ptr [rcx] ; load v.data[i] into xmm1
		mov rbx, [rax] ; load pointer to u.data* = m.data[i] into rbx
		mov r9, r12    ; load u.size as counter
		jloop:
			movss xmm1, dword ptr [rbx]	; load m.data[i][j] into xmm1
			movss xmm2, dword ptr [rdx] ; load result[i] into xmm2
			mulss xmm1, xmm0			; compute m.data[i][j] * v.data[i]
			addss xmm1, xmm2			; add product to result.data[j]
			movss dword ptr [rdx], xmm1 ; store in result[i]

			add rdx, 4  ; go to next element of result
			add rbx, 4  ; go to next element of u
			dec r9		; decrement col counter
		jnz jloop

		add rcx, 4 ; go to next element of v
		add rax, 8 ; go to next row in m
		dec r8; i--
	jnz iloop

	pop r12
	pop rbx ; restore rbx
	pop rax ; get pointer to result from stack
	ret
matTdotvec_unsafe_asm endp

; mat* outervec_asm(vec* u, vec* v, mat* result);
; TODO: return result pointer in rax
outervec_asm proc
	push r8 ; push pointer to result for return value
	mov r10, [rcx + 8] ; load u.size into r10 to use as counter
	mov r11, [rdx + 8] ; load v.size into r11 to reset v.size counter (r9)

	mov rcx, [rcx]     ; load u.data* into rcx to use as index
	mov rdx, [rdx]     ; load v.data* into rdx to reset v.data index (r8)

	mov rax, [r8]	   ; load result.data** into rax
	mov rax, [rax]	   ; load result.data* as vec into rax
	; Loop over rows ( copies ov v scaled by elements of u )
	rowloop:
		; movss xmm0, [rcx] ; load u.data[i] as float argument a 
		movss xmm0, dword ptr [rcx]
		mov r9, r11 ; load v.size into r10 to use as counter
		mov r8, rdx ; load v.data* into r8 to use as index
		colloop:
			movss xmm1, dword ptr [r8]; load next element of v
			
			mulss xmm1, xmm0 ; compute u.data[i] * v.data[j] and store in xmm1
			movss dword ptr [rax], xmm1 ; store in result.data[i][j]

			add rax, 4 ; step result index
			add r8 , 4 ; step v index
			dec r9     ; decrement loop counter over columns
		jnz colloop
		
		add rcx, 4 ; step u index
		dec r10    ; decrement loop counter over rows
	jnz rowloop
	pop rax ; get pointer to result for return from stack
	ret
outervec_asm endp

; uint64 equalsvec_asm(vec* u, vec* v);
; rcx = u.data*, rdx = v.data*,  r9 = u.size
equalsvec_asm proc
	call pass_two_vec_parameters
	call equalsvec_impl
	ret
equalsvec_asm endp

equalsvec_impl proc
	eqloop:
		movss xmm0, dword ptr [rcx] ; load u.data[i] to reg
		movss xmm1, dword ptr [rdx] ; load v.data[i] to reg
		comiss xmm0, xmm1			; sub u.data[i] and v.data[i]
		jne noteq
		add rcx, 4 ; step to next dword
		add rdx, 4
		dec r8
		jnz eqloop
	mov rax, 1
	ret
	 
	noteq:
	mov rax, 0
	ret

equalsvec_impl endp

; mat matalloc_asm(/*return pointer (hidden),*/ uint64 rows, uint64 cols);
; (hidden parameter) return address: rcx = vec
; rdx = rowsize
; r8 = colsize
matalloc_asm proc
	mov qword ptr [rcx + 8], rdx ; store rowsize in return mat
	mov qword ptr [rcx + 16], r8 ; store colsize in return mat
	
	push r12
	push r13
	push r14
	push rbx

	mov r12, rcx ; store mat pointer in nonvolatile register
	mov r13, rdx ; store rowsize in nonvolatile register to use as counter later
	mov r14, r8 ; store colsize in noncolatile resister to use as step size to initialize row pointers

	imul rcx, rdx, 8 ; move malloc parameter to first parameter register (8 bytes per address for x64)

	sub rsp, 20h
	call malloc_asm ; allocate row pointers
	add rsp, 20h

	mov rbx, rax ; store pointer to row pointers to nonvolatile register to use as index
	mov qword ptr [r12], rax ; store pointer to row pointers using return value from malloc
	
	; compute total size of matrix
	mov rax, 4   ; load sizeof(float)
	mul r14		 ; multiply by colsize 
	mul r13		 ; multiply by rowsize 
	mov rcx, rax ; pass as parameter to malloc

	sub rsp, 20h
	call malloc_asm ; allocate matrix elements
	add rsp, 20h

	xor rcx, rcx ; zero rcx to use as index

	rowallocloop:
		lea r8, qword ptr [rax + rcx * 4] ; get address of base of matrix elements (rax) offset by element index (rcx) times sizeof(float)
		mov qword ptr [rbx], r8 ; store address of row
		add rcx, r14; increase index (rcx) by colsize (r14)
		add rbx, 8 ; go to next row
		dec r13
	jnz rowallocloop
	
	mov rax, r12 ; return mat pointer in rax
	pop rbx
	pop r14
	pop r13
	pop r12
	
	ret
matalloc_asm endp

; mat matzeros_asm(/*return pointer (hidden),*/ uint64 rows, uint64 cols);
; (hidden parameter) return address: rcx = vec
; rdx = rowsize
; r8 = colsize
matzeros_asm proc
	mov qword ptr [rcx + 8], rdx ; store rowsize in return mat
	mov qword ptr [rcx + 16], r8 ; store colsize in return mat
	
	push r12
	push r13
	push r14
	push rbx

	mov r12, rcx ; store mat pointer in nonvolatile register
	mov r13, rdx ; store rowsize in nonvolatile register to use as counter later
	mov r14, r8 ; store colsize in noncolatile resister to use as step size to initialize row pointers

	imul rcx, rdx, 8 ; move malloc parameter to first parameter register (8 bytes per address for x64)

	sub rsp, 20h
	call malloc_asm ; allocate row pointers
	add rsp, 20h

	mov rbx, rax ; store pointer to row pointers to nonvolatile register to use as index
	mov qword ptr [r12], rax ; store pointer to row pointers using return value from malloc
	
	; compute total size of matrix
	mov rax, r14   ; load colsize 
	mul r13		 ; multiply by rowsize 
	mov rcx, rax ; pass as first parameter to calloc ( count / number of elements) 
	mov rdx, 4   ; pass as second parameter to calloc ( size in bytes per element)
	sub rsp, 20h
	call calloc_asm ; allocate zeroed matrix elements
	add rsp, 20h

	xor rcx, rcx ; zero rcx to use as index

	rowcallocloop:
		lea r8, qword ptr [rax + rcx * 4] ; get address of base of matrix elements (rax) offset by element index (rcx) times sizeof(float)
		mov qword ptr [rbx], r8 ; store address of row
		add rcx, r14; increase index (rcx) by colsize (r14)
		add rbx, 8 ; go to next row
		dec r13
	jnz rowcallocloop
	
	mov rax, r12 ; return mat pointer in rax
	pop rbx
	pop r14
	pop r13
	pop r12
	
	ret
matzeros_asm endp


; mat matones_asm(/*return pointer (hidden),*/ uint64 rows, uint64 cols)
matones_asm proc
	call matalloc_asm
	
	mov r11, rax ; save return value of matalloc_asm in r11

	mov rax, qword ptr [r11+8] ; rowsize in rax
	mul qword ptr [r11+10h] ; colsize 
	mov r8, rax
	mov r9, qword ptr [r11]   ; store rows ptr in r9
	mov r9, qword ptr [r9]   ; get ptr to first element of mat
	
	matonesloop:
		mov dword ptr [r9], 3f800000h ; Should be 32-bit single presicion floatingpoint 1.0f
		add r9, 4	
		dec r8
	jnz matonesloop 

	mov rax, r11; restore return value of matalloc
	ret
matones_asm endp

; mat matrand_asm(/*return pointer (hidden),*/ uint64 rows, uint64 cols);
matrand_asm proc
	call matalloc_asm

	push rax ; store pointer to returned struct on stack, to later pop to rax
	push r12 ; store in non-volatile regs to prevent interference from rand()
	push r13
	push r14
	
	mov r11, rax	; save return value of matalloc_asm in r11

	mov rax, qword ptr [r11 + 8]	; load rowsize
	mul qword ptr [r11 + 10h]		; rax = rowsize *  colsize
	mov r12, rax					; matrix size in r12
	mov r11, qword ptr [r11]	; get row pointers
	mov r13, qword ptr [r11]	; store data pointer in r13
	mov eax, 7fffh		; load rand max
	cvtsi2ss xmm0,  eax	; (float) RAND_MAX
	movd r14, xmm0		; store in non-volatile register
	matrandloop:
		call rand_asm		; load random int in eax
		cvtsi2ss xmm0, eax  ; convert random int to random float
		movd xmm1, r14      ; retrieve (float)RAND_MAX 
		divss xmm0, xmm1    ; scale to interval [0.f, 1.f)
		movss dword ptr [r13], xmm0 ;
		add r13, 4	
		dec r12
	jnz matrandloop 

	pop r14
	pop r13
	pop r12
	pop rax ; get stored pointer to vec from stack

	ret
matrand_asm endp

; mat matidentity_asm(/*return pointer (hidden),*/ uint64 rows, uint64 cols);
matidentity_asm proc
	call matzeros_asm
	
	mov r11, rax ; save return value of matzeros_asm in r11

	mov rax, qword ptr [r11 + 8] ; rowsize in rax
	mov rcx, qword ptr [r11 + 10h] ; colsize 
	mov r8, rax
	cmp rcx, rax
	cmovl r8, rcx 
	mov r9, qword ptr [r11]   ; store rows ptr in r9
	mov r9, qword ptr [r9]   ; get ptr to first element ofmat
	xor rdx, rdx
	matidloop:
		lea r10, [r9 + rdx * 4] ; 
		mov dword ptr [r10], 3f800000h ; Should be 32-bit single presicion floatingpoint 1.0f
		add r9, 4 ; move to next col not needed ?
		add rdx, rcx ; move to next row
		dec r8
	jnz matidloop 

	mov rax, r11; restore return value of matalloc
	ret
matidentity_asm endp

; void matfree_asm(mat* m);
matfree_asm proc
	mov qword ptr [rcx +  8], 0 ; set rowsize to zero
	mov qword ptr [rcx + 16], 0 ; set colsize to zero
	mov rcx, qword ptr[rcx] ; get pointer u.data*, whose memory we will free

	push rcx ; save pointer to row pointers to free later
	
	mov rcx, qword ptr[rcx] ; get pointer u.data, whose memory we will free
	sub rsp, 20h
	call free_asm ; free matrix element data
	add rsp, 20h
	
	pop rcx ; restore pointer to row pointers

	sub rsp, 20h
	call free_asm ; free row pointer data
	add rsp, 20h

	ret
matfree_asm endp

; mat* matdotmat_asm(mat* m, mat* n, mat* result);
; rcx = mat* m, rdx = mat* n, r8 = mat* result
; compute matrix product sum_j [i, j]x[j, k] = [i, k]
matdotmat_asm proc

	; prolog ( push callee save regs to stack )
	push r8 ; push result pointer to use as return value
	push r12
	push r13
	push r14
	push r15

	push rbp     ; store stack frame (base) pointer of calling procedure
	mov rbp, rsp ; set stack frame base

	; Initialize constants
	mov r11, [rcx +  8]; get left outer size  ( of dimension i ) Si = m.rowsize = result.rowsize
	mov r12, [rdx +  8]; get inner size	      ( of dimension j ) Sj = m.colsize = n.rowsize
	mov r13, [r8  + 16]; get right outer size ( of dimension k ) Sk = m.colsize = result.colsize

	mov r10, [r8]   ; Result.data
	mov r8 , [rcx]  ; M.data
	mov r9 , [rdx]  ; N.data

	; 'allocate' stack data for result.size dwords ( 4 * result.size bytes )
	mov rax, r11 ; mul operand 1 = Si
	mul r13      ; mul operand 2 = Sk
	mov r15, rax ; save number of dwords (floats) in result
	mov rdx, 4   ; 4 bytes per float
	mul rdx

	sub rsp, rax ; move stack pointer to make room for temp result
	
	mov r14, rsp ; write to stack via r14

	push r10 ; push pointer to result, to temporarily write to stack
	push r15 ; push allocated memory size for datacopy ( in number of dwords! )

	mov r10, r14 ; save pointer to stack data (temp result[0])

	; rax =  i, rcx =  j, rdx =  k
	; r11 = Si, r12 = Sj, r13 = Sk

	; fill result data with zeros
	; for i..Si
	xor rax, rax ; zero i
	zero_i_loop:
	;   for k..Sk
		xor rdx, rdx 
		zero_k_loop:
			; result[i][k] = 0
			
			;mov r14, qword ptr [r14]
			;lea r14, qword ptr [r14 + rdx * 4] ; get temp result[i][k]
			mov dword ptr [r14], 0

			add r14, 4
			inc rdx
			cmp rdx, r13 ; ++k <= Sk
		jl zero_k_loop
	; end for k
		inc rax
		cmp rax, r11 ; ++i <= Si
	jl zero_i_loop
	; end for i

	; rax, rcx, rdx, r8, r9, r10, r11

	; for i = 0..Si (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;   for k = 0..Sk (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;     for j = 0..Sj (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;       result[i][k] +=  m[i][j] * n[j][k]

	
	mov r14, r10; reset r14 to point to start of temp result.data on stack

		; for i..Si
	xor rax, rax ; zero i
	i_loop:
	;   for k..Sk
		xor rdx, rdx ; zero k
		k_loop:
			
			;lea r14, qword ptr [r10 + rax * 8] ; get temp &result[i]
			;mov r14, qword ptr [r14]
			;lea r14, qword ptr [r14 + rdx * 4] ; get temp &result[i][k]

			; for j = 0..Sj
			xor rcx, rcx ; zero j
			j_loop:
				; load m[i][j] into xmm0
				lea r15, qword ptr [r8  + rax * 8] ; get &m[i]
				mov r15, qword ptr [r15]
				lea r15, qword ptr [r15 + rcx * 4] ; get &m[i][j]
				movss xmm0, dword ptr [r15]

				; load n[j][k] into xmm1
				lea r15, qword ptr [r9  + rcx * 8] ; get &n[j]
				mov r15, qword ptr [r15]
				lea r15, qword ptr [r15 + rdx * 4] ; get &n[j][k]
				movss xmm1, dword ptr [r15]

				; compute m[i][j] * n[j][k]
				mulss xmm0, xmm1

				; temp result[i][k] += m[i][j] * n[j][k]
				movss xmm1, dword ptr [r14] ; load partial result
				addss xmm0, xmm1  ; add to running sum
				movss dword ptr [r14], xmm0 ; store updated partial result

				inc rcx
				cmp rcx, r12 ; ++j < Sj
			jl j_loop
			; end for j

			add r14, 4 ; Go to next element of result

			inc rdx
			cmp rdx, r13 ; ++k < Sk
		jl k_loop
	; end for k
		inc rax
		cmp rax, r11 ; ++i < Si
	jl i_loop
	; end for i

	pop rcx ; get size of stack allocated matrix to perform datacopy
	pop rax ; restore result.data pointer

	; perform data copy to result
	; data copy parameters
	mov rdx, r10				; load source address for datacpy ( data from stack )
	mov r8, qword ptr [rax]		; load destination address for datacpy ( result.data )

	CLD ; Clear direction flag to make sure we copy from lowest address to highest address. See https://en.wikipedia.org/wiki/Direction_flag

	call datacpy ; copy data from stack to result

	mov rsp, rbp ; reset stack frame base
	pop rbp      ; restore stack frame (base) pointer of calling procedure

	; epilog ( pop callee save regs from stack )
	pop r12
	pop r13
	pop r14
	pop r15

	; return result pointer in rax
	pop rax ; pop result return value

	ret
matdotmat_asm endp


; mat* matdotmat_asm(mat* m, mat* n, mat* result);
; rcx = mat* m, rdx = mat* n, r8 = mat* result
; compute matrix product sum_j [i, j]x[j, k] = [i, k]
; unsafe: output and input can not be the same matrices
matdotmat_unsafe_asm proc

	; prolog ( push callee save regs to stack )
	push r8 ; push result pointer to use as return value
	push r12
	push r13
	push r14
	push r15

	; Initialize constants
	mov r11, [rcx +  8]; get left outer size  ( of dimension i ) Si = m.rowsize = result.rowsize
	mov r12, [rdx +  8]; get inner size	      ( of dimension j ) Sj = m.colsize = n.rowsize
	mov r13, [r8  + 16]; get right outer size ( of dimension k ) Sk = m.colsize = result.colsize

	mov r10, [r8]   ; Result.data
	mov r8 , [rcx]  ; M.data
	mov r9 , [rdx]  ; N.data

	; rax =  i, rcx =  j, rdx =  k
	; r11 = Si, r12 = Sj, r13 = Sk

	; fill result data with zeros
	; for i..Si
	xor rax, rax ; zero i
	zero_i_loop:
	;   for k..Sk
		xor rdx, rdx 
		zero_k_loop:
			; result[i][k] = 0
			lea r14, qword ptr [r10  + rax * 8] ; get result[i]
			mov r14, qword ptr [r14]
			lea r14, qword ptr [r14 + rdx * 4] ; get result[i][k]
			mov dword ptr [r14], 0

			inc rdx
			cmp rdx, r13 ; ++k <= Sk
		jl zero_k_loop
		; end for k
		inc rax
		cmp rax, r11 ; ++i <= Si
	jl zero_i_loop
	; end for i

	; rax, rcx, rdx, r8, r9, r10, r11

	; for i = 0..Si (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;   for k = 0..Sk (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;     for j = 0..Sj (initialize / zero variable, add label, loop body, increment index, conditional jump)
	;       result[i][k] +=  m[i][j] * n[j][k]

		; for i..Si
	xor rax, rax ; zero i
	i_loop:
	;   for k..Sk
		xor rdx, rdx ; zero k
		k_loop:
			
			lea r14, qword ptr [r10 + rax * 8] ; get &result[i]
			mov r14, qword ptr [r14]
			lea r14, qword ptr [r14 + rdx * 4] ; get &result[i][k]

			; for j = 0..Sj
			xor rcx, rcx ; zero j
			j_loop:
				; load m[i][j] into xmm0
				lea r15, qword ptr [r8  + rax * 8] ; get &m[i]
				mov r15, qword ptr [r15]
				lea r15, qword ptr [r15 + rcx * 4] ; get &m[i][j]
				movss xmm0, dword ptr [r15]

				; load n[j][k] into xmm1
				lea r15, qword ptr [r9  + rcx * 8] ; get &n[j]
				mov r15, qword ptr [r15]
				lea r15, qword ptr [r15 + rdx * 4] ; get &n[j][k]
				movss xmm1, dword ptr [r15]

				; compute m[i][j] * n[j][k]
				mulss xmm0, xmm1

				; result[i][k] += m[i][j] * n[j][k]
				movss xmm1, dword ptr [r14] ; load partial result
				addss xmm0, xmm1  ; add to running sum
				movss dword ptr [r14], xmm0 ; store updated partial result

				inc rcx
				cmp rcx, r12 ; ++j < Sj
			jl j_loop
			; end for j

			inc rdx
			cmp rdx, r13 ; ++k < Sk
		jl k_loop
	; end for k
		inc rax
		cmp rax, r11 ; ++i < Si
	jl i_loop
	; end for i


	; epilog ( pop callee save regs from stack )
	pop r12
	pop r13
	pop r14
	pop r15

	; return result pointer in rax
	pop rax ; pop result return value

	ret
matdotmat_unsafe_asm endp

; mat* addmat_asm(mat* m, mat* n, mat* result);
addmat_asm proc
	push r8 ; push pointer to result for return value
	call pass_three_mat_parameters
	call convert_three_mat_to_vec_parameters ; The mat version of add assumes that matrix data is stored contiguously
	call addvec_impl
	pop rax ; get pointer to result for return from stack
	ret
addmat_asm endp

; mat* submat_asm(mat* m, mat* n, mat* result);
submat_asm proc
	push r8 ; push pointer to result for return value
	call pass_three_mat_parameters
	call convert_three_mat_to_vec_parameters ; The mat version of sub assumes that matrix data is stored contiguously
	call subvec_impl
	pop rax ; get pointer to result for return from stack
	ret
submat_asm endp

; mat* mulmat_asm(mat* m, mat* n, mat* result);
mulmat_asm proc
	push r8 ; push pointer to result for return value
	call pass_three_mat_parameters
	call convert_three_mat_to_vec_parameters ; The mat version of mul assumes that matrix data is stored contiguously
	call mulvec_impl
	pop rax ; get pointer to result for return from stack
	ret
mulmat_asm endp

; mat* divmat_asm(mat* m, mat* n, mat* result);
divmat_asm proc
	push r8 ; push pointer to result for return value
	call pass_three_mat_parameters
	call convert_three_mat_to_vec_parameters ; The mat version of div assumes that matrix data is stored contiguously
	call divvec_impl
	pop rax ; get pointer to result for return from stack
	ret
divmat_asm endp

; mat* scalemat_asm(float a, mat* m, mat* result);
scalemat_asm proc
	; float a in xmm0
	push r8 ; push pointer to result for return value
	mov r10, rdx					; store rdx, as it will get overwritten by mul
	mov rax, qword ptr [rdx + 10h]  ; Load col count of m into r10
	mov r9 , qword ptr [rdx +  8]   ; Load row count of m into r9
	mul r9
	mov r9 , rax					; load size of m.rows * m.cols
	mov rdx, qword ptr [r10]        ; (restore rdx from earlier and) load data pointer of m into rdx
	mov rax, r8						; store result pointer as return value
	mov r8 , qword ptr [r8 ]        ; load data** pointer of result into r8
	call sclvec					; r9 = u.size, rdx = u.data*, r8 = v.data*
	pop rax ; get pointer to result for return from stack
	ret
scalemat_asm endp

expmat_asm proc

	push rdx

	call pass_two_mat_parameters ; r9 = frommat.col, r8 = frommat.row , rcx  = frommat.data**, rdx = tomat.data**
	call convert_two_mat_to_vec_parameters ; The mat version of cpy assumes that matrix data is stored contiguously
	
	; Use callee save registers to avoid having to push and pop
	; registers before and after each call to exp
	push r12
	push r13
	push r14

	mov r14, qword ptr r8  ; r14  = u.size
	mov r13, qword ptr rcx	; r13 = u.data*
	mov r12, qword ptr rdx	; r12 = result.data*

	exploop:
		movss xmm0, dword ptr [r13] ; load u[i] as float argument (xmm0) to exp
		sub rsp, 20h ; reserve shadowspace
		call exp_asm
		add rsp, 20h ; deallocate shadowspace

		movss dword ptr [r12], xmm0 ; store result

		add r12, 4 ; result.data* += 1
		add r13, 4 ; u.data* += 1
		dec r14
	jnz exploop

	pop r14
	pop r13
	pop r12
	pop rax ; get result pointer as return value from first rdx push
	ret

	ret
expmat_asm endp

; mat* transposemat(mat* m, mat* result);
transposemat_asm proc
	push rdx ; store result return value to stack

	push r12 ;
	push r13 ;
	push r14 ;

	; Get dimensions
	mov r8, [rcx +  8] ; Load Si = m.rows = result.cols
	mov r9, [rcx + 16] ; Load Sk = m.cols = result.rows

	mov r14, r8    ; Store Si in r14
	cmp r9, r8     ; Compare Si and Sk
	cmovl r14, r9  ; Set r14 to min(Si, Sk);

	mov r10, qword ptr [rcx] ; Load m.data**
	mov r11, qword ptr [rdx] ; Load result.data**

	; rax =  i, rdx =  k
	; r8 = Si, r9 = Sk

	; for i..Si
	xor rax, rax ; zero i
	transpose_i_loop:
	;   for k..Sk
		mov rdx, rax ; only run for k = i, ..., Sk
		transpose_k_loop:
			; Get elements
			lea r12, qword ptr [r10  + rax * 8] ; get m.data[i]
			mov r12, qword ptr [r12]
			lea r12, qword ptr [r12 + rdx * 4] ; get m.data[i][k]
			movss xmm0, dword ptr [r12]

			lea r12, qword ptr [r10  + rdx * 8] ; get m.data[k]
			mov r12, qword ptr [r12]
			lea r12, qword ptr [r12 + rax * 4] ; get m.data[k][i]
			movss xmm1, dword ptr [r12]

			; transpose and store elements
			
			lea r13, qword ptr [r11 + rax * 8] ; get result.data[k]
			mov r13, qword ptr [r13]
			lea r13, qword ptr [r13 + rdx * 4] ; get result.data[i][k]
			movss dword ptr [r13], xmm1

			lea r13, qword ptr [r11 + rdx * 8] ; get result.data[k]
			mov r13, qword ptr [r13]
			lea r13, qword ptr [r13 + rax * 4] ; get result.data[k][i]
			movss dword ptr [r13], xmm0

			inc rdx
			cmp rdx, r14 ; ++k <= min(Si, Sk)
		jl transpose_k_loop
		; end for k
		inc rax
		cmp rax, r14 ; ++i <= min(Si, Sk)
	jl transpose_i_loop
	; end for i

	; for i..Si
	xor rax, rax ; zero i
	xor r14, r14
	cmp r8, r9
	cmovg rax, r9 ; i = (Si > Sj) ? Sj : 0
	
	cmp r8, r9
	cmovl r14, r8 ; k_start = (Si > Sj) ? Sj : 0
	transpose_i_loop2:
	;   for k..Sk
		mov rdx, r14 ; reset k = k_start
		transpose_k_loop2:
			; Get elements
			lea r12, qword ptr [r10  + rax * 8] ; get m.data[i]
			mov r12, qword ptr [r12]
			lea r12, qword ptr [r12 + rdx * 4] ; get m.data[i][k]
			movss xmm0, dword ptr [r12]

			; transpose and store elements
			lea r13, qword ptr [r11 + rdx * 8] ; get result.data[k]
			mov r13, qword ptr [r13]
			lea r13, qword ptr [r13 + rax * 4] ; get result.data[k][i]
			movss dword ptr [r13], xmm0

			inc rdx
			cmp rdx, r9 ; ++k <= Sk
		jl transpose_k_loop2
		; end for k
		inc rax
		cmp rax, r8 ; ++i <= Si
	jl transpose_i_loop2
	; end for i

	pop r14 ;
	pop r13 ;
	pop r12 ;

	pop rax ; load result return value from stack

	ret
transposemat_asm endp

; uint64 equalsmat_asm(mat* m, mat* n);
equalsmat_asm proc
	call pass_two_mat_parameters
	call convert_two_mat_to_vec_parameters ; The mat version of sub assumes that matrix data is stored contiguously
	call equalsvec_impl
	ret
equalsmat_asm endp


end