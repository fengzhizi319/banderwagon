
.macro reduce_macro ra0, ra1, ra2, ra3, rb0, rb1, rb2, rb3
    movq    \ra0, \rb0
    subq    q(%rip), \ra0
    movq    \ra1, \rb1
    sbbq    q+8(%rip), \ra1
    movq    \ra2, \rb2
    sbbq    q+16(%rip), \ra2
    movq    \ra3, \rb3
    sbbq    q+24(%rip), \ra3
    cmovc   \rb0, \ra0
    cmovc   \rb1, \ra1
    cmovc   \rb2, \ra2
    cmovc   \rb3, \ra3
.endm

.section .rodata
.global q
q:
    .quad 0xffffffff00000001
    .quad 0x53bda402fffe5bfe
    .quad 0x3339d80809a1d805
    .quad 0x73eda753299d7d48

.section .rodata
.global q_inv
q_inv:
    .quad 0xfffffffeffffffff

.section .text
.global get_q

get_q:
    movq q(%rip), %rax
    ret

// parameters: The first six integer parameters are passed through the 
// rdi, rsi, rdx, rcx, r8, and r9 registers in order.
.section .text
.global fr_add

fr_add:
	push %rcx
	push %rbx
	push %rax
	push %r12
	push %r8
	push %r9
	push %r10
	push %r11
	
    movq 0(%rsi), %rcx
    movq 8(%rsi), %rbx
    movq 16(%rsi), %rax
    movq 24(%rsi), %r12

    addq 0(%rdx), %rcx
    adcq 8(%rdx), %rbx
    adcq 16(%rdx), %rax
    adcq 24(%rdx), %r12

    reduce_macro %rcx, %rbx, %rax, %r12, %r8, %r9, %r10, %r11

    movq %rcx, 0(%rdi)
    movq %rbx, 8(%rdi)
    movq %rax, 16(%rdi)
    movq %r12, 24(%rdi)

	pop %r11
	pop %r10
	pop %r9
	pop %r8
	pop %r12
	pop %rax
	pop %rbx
	pop %rcx

    ret

.section .text
.global mul_by_5

mul_by_5:
	push %rdx
	push %rcx
	push %rbx
	push %rsi
	push %r8
	push %r9
	push %r10
	push %rax
	
    movq 0(%rdi), %rdx
    movq 8(%rdi), %rcx
    movq 16(%rdi), %rbx
    movq 24(%rdi), %rsi

    addq %rdx, %rdx
    adcq %rcx, %rcx
    adcq %rbx, %rbx
    adcq %rsi, %rsi

	reduce_macro %rdx, %rcx, %rbx, %rsi, %rax, %r8, %r9, %r10

    addq %rdx, %rdx
    adcq %rcx, %rcx
    adcq %rbx, %rbx
    adcq %rsi, %rsi

	reduce_macro %rdx, %rcx, %rbx, %rsi, %rax, %r8, %r9, %r10

    addq 0(%rdi), %rdx
    adcq 8(%rdi), %rcx
    adcq 16(%rdi), %rbx
    adcq 24(%rdi), %rsi

	reduce_macro %rdx, %rcx, %rbx, %rsi, %rax, %r8, %r9, %r10

    movq %rdx, 0(%rdi)
    movq %rcx, 8(%rdi)
    movq %rbx, 16(%rdi)
    movq %rsi, 24(%rdi)

	pop %rax
	pop %r10
	pop %r9
	pop %r8
	pop %rsi
	pop %rbx
	pop %rcx
	pop %rdx
	ret

// rdi, rsi, rdx, rcx, r8, and r9 registers in order.
.section .text
.global fr_mul
fr_mul:
    push %rbp
	push %rcx
	push %rbx
	push %r8
	push %r9
	push %r10
	push %rax
	push %r11
	push %r12
	push %r13
	push %r14
	push %r15
	
    movq 0(%rsi), %r12
    movq 8(%rsi), %r8
    movq 16(%rsi), %r9
    movq 24(%rsi), %r10
    
    movq %rdx, %r11

    # clear the flags
    xorq %rax, %rax
    movq 0(%r11), %rdx
	

	// (A,t[0])  := x[0]*y[0] + A
	mulxq %r12, %r14, %r13

	// (A,t[1])  := x[1]*y[0] + A
	mulxq %r8, %rax, %rcx
	adoxq %rax, %r13

	// (A,t[2])  := x[2]*y[0] + A
	mulxq %r9, %rax, %rbx
	adoxq %rax, %rcx

	// (A,t[3])  := x[3]*y[0] + A
	mulxq %r10, %rax, %rbp
	adoxq %rax, %rbx

	// A += carries from ADCXQ and ADOXQ
	movq $0, %rax
	adoxq %rax, %rbp

	// m := t[0]*q'[0] mod W
	movq q_inv(%rip), %rdx
	imulq %r14, %rdx

	// clear the flags
	xorq %rax, %rax

	// C,_ := t[0] + m*q[0]
	mulxq q(%rip), %rax, %rsi
	adcxq %r14, %rax
	movq %rsi, %r14

	// (C,t[0]) := t[1] + m*q[1] + C
	adcxq %r13, %r14
	mulxq q+8(%rip), %rax, %r13
	adoxq %rax, %r14

	// (C,t[1]) := t[2] + m*q[2] + C
	adcxq %rcx, %r13
	mulxq q+16(%rip), %rax, %rcx
	adoxq %rax, %r13

	// (C,t[2]) := t[3] + m*q[3] + C
	adcxq %rbx, %rcx
	mulxq q+24(%rip), %rax, %rbx
	adoxq %rax, %rcx

	// t[3] = C + A
	movq $0, %rax
	adcxq %rax, %rbx
	adoxq %rbp, %rbx

	// clear the flags
	xorq %rax, %rax
	movq 8(%r11), %rdx

	// (A,t[0])  := t[0] + x[0]*y[1] + A
	mulxq %r12, %rax, %rbp
	adoxq %rax, %r14

	// (A,t[1])  := t[1] + x[1]*y[1] + A
	adcxq %rbp, %r13
    mulxq %r8, %rax, %rbp
    adoxq %rax, %r13

	// (A,t[2])  := t[2] + x[2]*y[1] + A
	adcxq %rbp, %rcx
	mulxq %r9, %rax, %rbp
	adoxq %rax, %rcx

	// (A,t[3])  := t[3] + x[3]*y[1] + A
	adcxq %rbp, %rbx
	mulxq %r10, %rax, %rbp
	adoxq %rax, %rbx

	// A += carries from ADCXQ and ADOXQ
	movq $0, %rax
	adcxq %rax, %rbp
	adoxq %rax, %rbp

	// m := t[0]*q'[0] mod W
	movq q_inv(%rip), %rdx
	imulq %r14, %rdx

	// clear the flags
	xorq %rax, %rax

	// C,_ := t[0] + m*q[0]
	mulxq q(%rip), %rax, %rsi
	adcxq %r14, %rax
	movq %rsi, %r14

	// (C,t[0]) := t[1] + m*q[1] + C
	adcxq %r13, %r14
	mulxq q+8(%rip), %rax, %r13
	adoxq %rax, %r14

	// (C,t[1]) := t[2] + m*q[2] + C
	adcxq %rcx, %r13
	mulxq q+16(%rip), %rax, %rcx
	adoxq %rax, %r13

	// (C,t[2]) := t[3] + m*q[3] + C
	adcxq %rbx, %rcx
	mulxq q+24(%rip), %rax, %rbx
	adoxq %rax, %rcx

	// t[3] = C + A
	movq $0, %rax
	adcxq %rax, %rbx
	adoxq %rbp, %rbx

	// clear the flags
	xorq %rax, %rax
	movq 16(%r11), %rdx

	// (A,t[0])  := t[0] + x[0]*y[2] + A
	mulxq %r12, %rax, %rbp
	adoxq %rax, %r14

	// (A,t[1])  := t[1] + x[1]*y[2] + A
	adcxq %rbp, %r13
	mulxq %r8, %rax, %rbp
	adoxq %rax, %r13

	// (A,t[2])  := t[2] + x[2]*y[2] + A
	adcxq %rbp, %rcx
	mulxq %r9, %rax, %rbp
	adoxq %rax, %rcx

	// (A,t[3])  := t[3] + x[3]*y[2] + A
	adcxq %rbp, %rbx
	mulxq %r10, %rax, %rbp
	adoxq %rax, %rbx

	// A += carries from ADCXQ and ADOXQ
	movq  $0, %rax
	adcxq %rax, %rbp
	adoxq %rax, %rbp

	// m := t[0]*q'[0] mod W
	movq  q_inv(%rip), %rdx
	imulq %r14, %rdx

	// clear the flags
	xorq %rax, %rax

	// C,_ := t[0] + m*q[0]
	mulxq q(%rip), %rax, %rsi
	adcxq %r14, %rax
	movq  %rsi, %r14

	// (C,t[0]) := t[1] + m*q[1] + C
	adcxq %r13, %r14
	mulxq q+8(%rip), %rax, %r13
	adoxq %rax, %r14

	// (C,t[1]) := t[2] + m*q[2] + C
	adcxq %rcx, %r13
	mulxq q+16(%rip), %rax, %rcx
	adoxq %rax, %r13

	// (C,t[2]) := t[3] + m*q[3] + C
	adcxq %rbx, %rcx
	mulxq q+24(%rip), %rax, %rbx
	adoxq %rax, %rcx

	// t[3] = C + A
	movq  $0, %rax
	adcxq %rax, %rbx
	adoxq %rbp, %rbx

	// clear the flags
	xorq %rax, %rax
	movq 24(%r11), %rdx

	// (A,t[0])  := t[0] + x[0]*y[3] + A
	mulxq %r12, %rax, %rbp
	adoxq %rax, %r14

	// (A,t[1])  := t[1] + x[1]*y[3] + A
	adcxq %rbp, %r13
	mulxq %r8, %rax, %rbp
	adoxq %rax, %r13

	// (A,t[2])  := t[2] + x[2]*y[3] + A
	adcxq %rbp, %rcx
	mulxq %r9, %rax, %rbp
	adoxq %rax, %rcx

	// (A,t[3])  := t[3] + x[3]*y[3] + A
	adcxq %rbp, %rbx
	mulxq %r10, %rax, %rbp
	adoxq %rax, %rbx

	// A += carries from ADCXQ and ADOXQ
	movq  $0, %rax
	adcxq %rax, %rbp
	adoxq %rax, %rbp

	// m := t[0]*q'[0] mod W
	movq  q_inv(%rip), %rdx
	imulq %r14, %rdx

	// clear the flags
	xorq %rax, %rax

	// C,_ := t[0] + m*q[0]
	mulxq q(%rip), %rax, %rsi
	adcxq %r14, %rax
	movq  %rsi, %r14

	// (C,t[0]) := t[1] + m*q[1] + C
	adcxq %r13, %r14
	mulxq q+8(%rip), %rax, %r13
	adoxq %rax, %r14

	// (C,t[1]) := t[2] + m*q[2] + C
	adcxq %rcx, %r13
	mulxq q+16(%rip), %rax, %rcx
	adoxq %rax, %r13

	// (C,t[2]) := t[3] + m*q[3] + C
	adcxq %rbx, %rcx
	mulxq q+24(%rip), %rax, %rbx
	adoxq %rax, %rcx

	// t[3] = C + A
	movq  $0, %rax
	adcxq %rax, %rbx
	adoxq %rbp, %rbx

	// reduce element(R14,R13,CX,BX) using temp registers (SI,R12,R11,AX)
	reduce_macro %r14, %r13, %rcx, %rbx, %rsi, %r12, %r11, %rax

	movq %r14, 0(%rdi)
	movq %r13, 8(%rdi)
	movq %rcx, 16(%rdi)
	movq %rbx, 24(%rdi)
	
	pop %r15
	pop %r14
	pop %r13
	pop %r12
	pop %r11
	pop %rax
	pop %r10
	pop %r9
	pop %r8
	pop %rbx
	pop %rcx
	pop %rbp
	ret

