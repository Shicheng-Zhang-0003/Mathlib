	.file	"fma_audit.c"
	.text
	.p2align 4
	.globl	audit_exp
	.type	audit_exp, @function
audit_exp:
.LFB24:
	.cfi_startproc
	endbr64
	vxorpd	%xmm1, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	jp	.L12
	je	.L10
.L12:
	vcomisd	.LC3(%rip), %xmm0
	ja	.L20
	vmovsd	.LC4(%rip), %xmm2
	vxorpd	%xmm1, %xmm1, %xmm1
	vcomisd	%xmm0, %xmm2
	ja	.L1
	vdivsd	.LC5(%rip), %xmm0, %xmm2
	vmovq	%xmm2, %rdx
	movabsq	$9218868437227405312, %rax
	notq	%rdx
	testq	%rax, %rdx
	jne	.L21
.L6:
	vmulsd	.LC7(%rip), %xmm2, %xmm1
	leaq	144+inv_fact.0(%rip), %rax
	leaq	-144(%rax), %rdx
	vsubsd	%xmm1, %xmm0, %xmm0
	vmulsd	.LC8(%rip), %xmm2, %xmm1
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	.LC0(%rip), %xmm1
	.p2align 4,,10
	.p2align 3
.L9:
	vfmadd213sd	(%rax), %xmm0, %xmm1
	subq	$16, %rax
	vfmadd213sd	8(%rax), %xmm0, %xmm1
	cmpq	%rax, %rdx
	jne	.L9
	vfmadd213sd	.LC1(%rip), %xmm1, %xmm0
	vcvttsd2sil	%xmm2, %eax
	salq	$52, %rax
	vmovq	%xmm0, %rcx
	addq	%rcx, %rax
	vmovq	%rax, %xmm1
	vmovsd	%xmm1, %xmm1, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L20:
	vmovsd	.LC1(%rip), %xmm1
	vxorpd	%xmm0, %xmm0, %xmm0
	vdivsd	%xmm0, %xmm1, %xmm1
.L1:
	vmovsd	%xmm1, %xmm1, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L10:
	vmovsd	.LC1(%rip), %xmm1
	vmovsd	%xmm1, %xmm1, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L21:
	vcomisd	%xmm1, %xmm2
	vxorps	%xmm3, %xmm3, %xmm3
	jb	.L18
	vaddsd	.LC6(%rip), %xmm2, %xmm2
	vcvttsd2siq	%xmm2, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm2
	jmp	.L6
.L18:
	vsubsd	.LC6(%rip), %xmm2, %xmm2
	vcvttsd2siq	%xmm2, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm2
	jmp	.L6
	.cfi_endproc
.LFE24:
	.size	audit_exp, .-audit_exp
	.p2align 4
	.globl	audit_sin
	.type	audit_sin, @function
audit_sin:
.LFB25:
	.cfi_startproc
	endbr64
	movabsq	$9218868437227405312, %rax
	vmovq	%xmm0, %rdx
	andq	%rax, %rdx
	cmpq	%rax, %rdx
	jne	.L23
	movabsq	$4503599627370495, %rcx
	vxorpd	%xmm0, %xmm0, %xmm0
	vdivsd	%xmm0, %xmm0, %xmm0
	vmovq	%xmm0, %rax
	testq	%rcx, %rax
	je	.L48
	notq	%rax
	testq	%rdx, %rax
	jne	.L48
	ret
	.p2align 4,,10
	.p2align 3
.L23:
	vmulsd	.LC10(%rip), %xmm0, %xmm1
	vmovq	%xmm1, %rdx
	notq	%rdx
	testq	%rax, %rdx
	jne	.L49
.L27:
	vmulsd	.LC13(%rip), %xmm1, %xmm2
	vcvttsd2sil	%xmm1, %eax
	movabsq	$9218868437227405312, %rcx
	andl	$3, %eax
	vsubsd	%xmm2, %xmm0, %xmm0
	vmulsd	.LC33(%rip), %xmm1, %xmm2
	vaddsd	%xmm2, %xmm0, %xmm0
	vmovq	%xmm0, %rdx
	movq	%rdx, %rsi
	notq	%rsi
	testq	%rcx, %rsi
	jne	.L30
	movabsq	$4503599627370495, %rcx
	testq	%rcx, %rdx
	jne	.L50
.L30:
	vmulsd	%xmm0, %xmm0, %xmm1
	cmpl	$2, %eax
	je	.L32
	cmpl	$3, %eax
	je	.L33
	cmpl	$1, %eax
	jne	.L26
	vmovsd	.LC16(%rip), %xmm0
	vfmadd213sd	.LC17(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC18(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC19(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC20(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC21(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC22(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC23(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC24(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC1(%rip), %xmm1, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L33:
	vmovsd	.LC16(%rip), %xmm0
	vfmadd213sd	.LC17(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC18(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC19(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC20(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC21(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC22(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC23(%rip), %xmm1, %xmm0
	vfmadd213sd	.LC24(%rip), %xmm1, %xmm0
	vfnmsub213sd	.LC1(%rip), %xmm1, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L48:
	vmulsd	%xmm0, %xmm0, %xmm1
.L26:
	leaq	64+minimax_sin_coeffs(%rip), %rax
	vmovsd	.LC9(%rip), %xmm2
	leaq	-64(%rax), %rcx
	.p2align 4,,10
	.p2align 3
.L35:
	movq	%rax, %rdx
	vfmadd213sd	(%rax), %xmm1, %xmm2
	subq	$8, %rax
	cmpq	%rcx, %rdx
	jne	.L35
	vmulsd	%xmm2, %xmm0, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L49:
	vcomisd	.LC11(%rip), %xmm1
	ja	.L27
	vmovsd	.LC12(%rip), %xmm2
	vcomisd	%xmm1, %xmm2
	ja	.L27
	vxorpd	%xmm2, %xmm2, %xmm2
	vxorps	%xmm3, %xmm3, %xmm3
	vcomisd	%xmm2, %xmm1
	jb	.L47
	vaddsd	.LC6(%rip), %xmm1, %xmm1
	vcvttsd2siq	%xmm1, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm1
	jmp	.L27
	.p2align 4,,10
	.p2align 3
.L32:
	vmovsd	.LC9(%rip), %xmm2
	vfmadd213sd	.LC25(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC26(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC27(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC28(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC29(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC30(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC31(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC32(%rip), %xmm1, %xmm2
	vfmadd213sd	.LC1(%rip), %xmm2, %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vxorpd	.LC14(%rip), %xmm0, %xmm0
	ret
	.p2align 4,,10
	.p2align 3
.L47:
	vsubsd	.LC6(%rip), %xmm1, %xmm1
	vcvttsd2siq	%xmm1, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm1
	jmp	.L27
.L50:
	vxorpd	%xmm0, %xmm0, %xmm0
	vdivsd	%xmm0, %xmm0, %xmm0
	ret
	.cfi_endproc
.LFE25:
	.size	audit_sin, .-audit_sin
	.section	.rodata
	.align 32
	.type	inv_fact.0, @object
	.size	inv_fact.0, 160
inv_fact.0:
	.long	0
	.long	1072693248
	.long	0
	.long	1072693248
	.long	0
	.long	1071644672
	.long	1431655765
	.long	1069897045
	.long	1431655765
	.long	1067799893
	.long	286331153
	.long	1065423121
	.long	381774871
	.long	1062650220
	.long	436314138
	.long	1059717536
	.long	436314138
	.long	1056571808
	.long	-1521039564
	.long	1053236707
	.long	-1216831652
	.long	1049787983
	.long	1744127204
	.long	1046144581
	.long	-268904296
	.long	1042411224
	.long	329805065
	.long	1038488134
	.long	-1463780195
	.long	1034500468
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1026222067
	.long	1882238282
	.long	1021924039
	.long	1673100695
	.long	1017545336
	.long	1182875991
	.long	1013118107
	.align 32
	.type	minimax_sin_coeffs, @object
	.size	minimax_sin_coeffs, 80
minimax_sin_coeffs:
	.long	0
	.long	1072693248
	.long	1431655765
	.long	-1077586603
	.long	286331153
	.long	1065423121
	.long	436314138
	.long	-1087766112
	.long	-1521039564
	.long	1053236707
	.long	1744127204
	.long	-1101339067
	.long	329805065
	.long	1038488134
	.long	-416040929
	.long	-1117067277
	.long	1882238283
	.long	1021924039
	.long	1552785989
	.long	-1134365541
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	1182875991
	.long	1013118107
	.align 8
.LC1:
	.long	0
	.long	1072693248
	.align 8
.LC3:
	.long	1889785610
	.long	1082535485
	.align 8
.LC4:
	.long	1030792151
	.long	-1064875766
	.align 8
.LC5:
	.long	-17155601
	.long	1072049730
	.align 8
.LC6:
	.long	0
	.long	1071644672
	.align 8
.LC7:
	.long	-18874368
	.long	1072049730
	.align 8
.LC8:
	.long	821600374
	.long	1038760431
	.align 8
.LC9:
	.long	1552785989
	.long	-1134365541
	.align 8
.LC10:
	.long	1841940611
	.long	1071931184
	.align 8
.LC11:
	.long	1235175040
	.long	1138752769
	.align 8
.LC12:
	.long	1235175040
	.long	-1008730879
	.align 8
.LC13:
	.long	1413754136
	.long	1073291771
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC14:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.section	.rodata.cst8
	.align 8
.LC16:
	.long	1672545861
	.long	-1129938312
	.align 8
.LC17:
	.long	-416040929
	.long	1026222067
	.align 8
.LC18:
	.long	-1463780195
	.long	-1112983180
	.align 8
.LC19:
	.long	-268904296
	.long	1042411224
	.align 8
.LC20:
	.long	-1216831652
	.long	-1097695665
	.align 8
.LC21:
	.long	436314138
	.long	1056571808
	.align 8
.LC22:
	.long	381774871
	.long	-1084833428
	.align 8
.LC23:
	.long	1431655765
	.long	1067799893
	.align 8
.LC24:
	.long	0
	.long	-1075838976
	.align 8
.LC25:
	.long	1882238283
	.long	1021924039
	.align 8
.LC26:
	.long	-416040929
	.long	-1117067277
	.align 8
.LC27:
	.long	329805065
	.long	1038488134
	.align 8
.LC28:
	.long	1744127204
	.long	-1101339067
	.align 8
.LC29:
	.long	-1521039564
	.long	1053236707
	.align 8
.LC30:
	.long	436314138
	.long	-1087766112
	.align 8
.LC31:
	.long	286331153
	.long	1065423121
	.align 8
.LC32:
	.long	1431655765
	.long	-1077586603
	.align 8
.LC33:
	.long	856972295
	.long	-1131305434
	.ident	"GCC: (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
