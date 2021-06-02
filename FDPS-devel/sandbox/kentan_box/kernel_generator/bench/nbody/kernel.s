	.text
	.file	"kernel.cpp"
	.globl	_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force // -- Begin function _Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force
	.p2align	3
	.type	_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force,@function
_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force: // @_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force
.Lfunc_begin0:
	.file	1 "/home/users/go19/go0011/kernel_generator/bench/nbody" "kernel.cpp"
	.loc	1 5 0                   // kernel.cpp:5:0
	.cfi_startproc
// %bb.0:
	addvl	sp, sp, #-10
	stp	d15, d14, [sp, #-80]!   // 16-byte Folded Spill
	stp	d13, d12, [sp, #16]     // 16-byte Folded Spill
	stp	d11, d10, [sp, #32]     // 16-byte Folded Spill
	stp	d9, d8, [sp, #48]       // 16-byte Folded Spill
	str	x28, [sp, #64]          // 8-byte Folded Spill
	.cfi_escape 0x0f, 0x0f, 0x92, 0x1f, 0x00, 0x11, 0xd0, 0x00, 0x22, 0x11, 0xd0, 0x00, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfa = sp + 80 + 80 * VG
	.cfi_escape 0x10, 0x1c, 0x0b, 0x11, 0x70, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(X28) = cfa + -16 + -80 * VG
	.cfi_escape 0x10, 0x48, 0x0b, 0x11, 0x68, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D8) = cfa + -24 + -80 * VG
	.cfi_escape 0x10, 0x49, 0x0b, 0x11, 0x60, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D9) = cfa + -32 + -80 * VG
	.cfi_escape 0x10, 0x4a, 0x0b, 0x11, 0x58, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D10) = cfa + -40 + -80 * VG
	.cfi_escape 0x10, 0x4b, 0x0b, 0x11, 0x50, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D11) = cfa + -48 + -80 * VG
	.cfi_escape 0x10, 0x4c, 0x0b, 0x11, 0x48, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D12) = cfa + -56 + -80 * VG
	.cfi_escape 0x10, 0x4d, 0x0b, 0x11, 0x40, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D13) = cfa + -64 + -80 * VG
	.cfi_escape 0x10, 0x4e, 0x0c, 0x11, 0xb8, 0x7f, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D14) = cfa + -72 + -80 * VG
	.cfi_escape 0x10, 0x4f, 0x0c, 0x11, 0xb0, 0x7f, 0x22, 0x11, 0xb0, 0x7f, 0x92, 0x2e, 0x00, 0x1e, 0x22 // cfi(D15) = cfa + -80 + -80 * VG
.Ltmp0:
	.loc	1 12 1 prologue_end     // kernel.cpp:12:1
	cmp	w1, #1                  // =1
	b.lt	.LBB0_18
// %bb.1:
	.loc	1 12 19 is_stmt 0       // kernel.cpp:12:19
	add	w8, w1, #15             // =15
	.loc	1 12 23                 // kernel.cpp:12:23
	add	w9, w1, #30             // =30
	cmp	w8, #0                  // =0
	csel	w8, w9, w8, lt
	and	w8, w8, #0xfffffff0
	.loc	1 13 7 is_stmt 1        // kernel.cpp:13:7
	sxtw	x8, w8
	.loc	1 12 1                  // kernel.cpp:12:1
	cmp	w3, #0                  // =0
	b.le	.LBB0_11
// %bb.2:
	.loc	1 0 1 is_stmt 0         // kernel.cpp:0:1
	mov	z26.s, #0               // =0x0
	fmov	z16.s, #1.00000000
	mov	w12, w3
	ptrue	p0.s
	and	x11, x12, #0x3
	fmov	z17.s, #0.37500000
	mov	x9, xzr
	.loc	1 13 7 is_stmt 1        // kernel.cpp:13:7
	sub	x10, x12, #1            // =1
	sub	x12, x12, x11
	fmov	z27.s, #0.50000000
	add	x13, x2, #64            // =64
	add	x14, x2, #16            // =16
	neg	x15, x11
	.p2align	2
.LBB0_3:                                // =>This Loop Header: Depth=1
                                        //     Child Loop BB0_6 Depth 2
                                        //     Child Loop BB0_9 Depth 2
	whilelt	p1.s, w9, w1
	.loc	1 14 31                 // kernel.cpp:14:31
	lsl	x17, x9, #3
	.loc	1 20 33                 // kernel.cpp:20:33
	add	x16, x4, x9, lsl #4
	.loc	1 14 31                 // kernel.cpp:14:31
	ld4w	{ z22.s, z23.s, z24.s, z25.s }, p1/z, [x0, x17, lsl #2]
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z4.s, z5.s, z6.s, z7.s }, p1/z, [x16]
	mov	z21.d, z5.d
	mov	z20.d, z6.d
	mov	z18.d, z7.d
	sel	z28.s, p1, z22.s, z26.s
	sel	z29.s, p1, z23.s, z26.s
	sel	z30.s, p1, z24.s, z26.s
	sel	z25.s, p1, z25.s, z26.s
	.loc	1 28 1                  // kernel.cpp:28:1
	cmp	x10, #3                 // =3
	b.hs	.LBB0_5
// %bb.4:                               //   in Loop: Header=BB0_3 Depth=1
	.loc	1 0 1 is_stmt 0         // kernel.cpp:0:1
	mov	x17, xzr
	.loc	1 28 1                  // kernel.cpp:28:1
	cbnz	x11, .LBB0_8
	b	.LBB0_10
	.p2align	2
.LBB0_5:                                //   in Loop: Header=BB0_3 Depth=1
	.loc	1 0 1                   // kernel.cpp:0:1
	add	x28, sp, #80            // =80
	str	z28, [x28, #2, mul vl]  // 16-byte Folded Spill
	add	x28, sp, #80            // =80
	str	z29, [x28, #3, mul vl]  // 16-byte Folded Spill
	add	x28, sp, #80            // =80
	mov	x17, xzr
	mov	x18, x13
	str	z30, [x28, #4, mul vl]  // 16-byte Folded Spill
	add	x28, sp, #80            // =80
	str	z25, [x28, #5, mul vl]  // 16-byte Folded Spill
	.p2align	2
.LBB0_6:                                //   Parent Loop BB0_3 Depth=1
                                        // =>  This Inner Loop Header: Depth=2
	.loc	1 30 30 is_stmt 1       // kernel.cpp:30:30
	add	x28, sp, #80            // =80
	str	z20, [x28, #1, mul vl]  // 16-byte Folded Spill
	add	x28, sp, #80            // =80
	str	z4, [x28, #6, mul vl]
	str	z5, [x28, #7, mul vl]
	str	z6, [x28, #8, mul vl]
	str	z7, [x28, #9, mul vl]
	add	x28, sp, #80            // =80
	str	z21, [x28]              // 16-byte Folded Spill
	ldp	d23, d24, [x18, #-64]
	.loc	1 30 19 is_stmt 0       // kernel.cpp:30:19
	fcvt	s23, d23
	.loc	1 31 19 is_stmt 1       // kernel.cpp:31:19
	fcvt	s26, d24
	.loc	1 30 9                  // kernel.cpp:30:9
	mov	z23.s, s23
	.loc	1 32 30                 // kernel.cpp:32:30
	ldur	d24, [x18, #-48]
	.loc	1 32 19 is_stmt 0       // kernel.cpp:32:19
	fcvt	s28, d24
	add	x28, sp, #80            // =80
	.loc	1 34 23 is_stmt 1       // kernel.cpp:34:23
	ldp	s24, s25, [x18, #-40]
	.loc	1 36 9                  // kernel.cpp:36:9
	mov	z25.s, s25
	.loc	1 31 9                  // kernel.cpp:31:9
	mov	z26.s, s26
	fmov	z17.s, #1.00000000
	.loc	1 30 30                 // kernel.cpp:30:30
	ldp	d27, d29, [x18, #-32]
	.loc	1 30 19 is_stmt 0       // kernel.cpp:30:19
	fcvt	s27, d27
	.loc	1 30 9                  // kernel.cpp:30:9
	mov	z27.s, s27
	.loc	1 31 19 is_stmt 1       // kernel.cpp:31:19
	fcvt	s31, d29
	.loc	1 32 9                  // kernel.cpp:32:9
	mov	z29.s, s28
	.loc	1 34 23                 // kernel.cpp:34:23
	ldp	s30, s28, [x18, #-8]
	.loc	1 36 9                  // kernel.cpp:36:9
	mov	z28.s, s28
	fmov	z16.s, #0.50000000
	.loc	1 30 30                 // kernel.cpp:30:30
	ldp	d8, d9, [x18]
	.loc	1 31 9                  // kernel.cpp:31:9
	mov	z10.s, s31
	.loc	1 30 19                 // kernel.cpp:30:19
	fcvt	s31, d8
	.loc	1 30 9 is_stmt 0        // kernel.cpp:30:9
	mov	z8.s, s31
	.loc	1 32 30 is_stmt 1       // kernel.cpp:32:30
	ldur	d31, [x18, #-16]
	.loc	1 32 19 is_stmt 0       // kernel.cpp:32:19
	fcvt	s31, d31
	.loc	1 31 19 is_stmt 1       // kernel.cpp:31:19
	fcvt	s9, d9
	.loc	1 32 9                  // kernel.cpp:32:9
	mov	z11.s, s31
	.loc	1 31 9                  // kernel.cpp:31:9
	mov	z9.s, s9
	.loc	1 32 30                 // kernel.cpp:32:30
	ldr	d31, [x18, #16]
	.loc	1 32 19 is_stmt 0       // kernel.cpp:32:19
	fcvt	s12, d31
	.loc	1 34 6 is_stmt 1        // kernel.cpp:34:6
	ld1rw	{ z31.s }, p0/z, [x18, #24]
	.loc	1 36 9                  // kernel.cpp:36:9
	ld1rw	{ z13.s }, p0/z, [x18, #28]
	.loc	1 32 9                  // kernel.cpp:32:9
	mov	z12.s, s12
	fmov	z20.s, #0.37500000
	.loc	1 28 19                 // kernel.cpp:28:19
	add	x17, x17, #4            // =4
	.loc	1 30 30                 // kernel.cpp:30:30
	ldp	d14, d15, [x18, #32]
	.loc	1 30 19 is_stmt 0       // kernel.cpp:30:19
	fcvt	s14, d14
	.loc	1 30 9                  // kernel.cpp:30:9
	mov	z14.s, s14
	ldr	z19, [x28, #2, mul vl]  // 16-byte Folded Reload
	add	x28, sp, #80            // =80
	.loc	1 37 10 is_stmt 1       // kernel.cpp:37:10
	mov	z22.d, z19.d
	fsub	z22.s, p1/m, z22.s, z23.s
	.loc	1 34 6                  // kernel.cpp:34:6
	mov	z23.s, s24
	ldr	z24, [x28, #3, mul vl]  // 16-byte Folded Reload
	add	x28, sp, #80            // =80
	.loc	1 38 10                 // kernel.cpp:38:10
	mov	z21.d, z24.d
	fsub	z21.s, p1/m, z21.s, z26.s
	ldr	z26, [x28, #4, mul vl]  // 16-byte Folded Reload
	add	x28, sp, #80            // =80
	ldr	z4, [x28, #5, mul vl]   // 16-byte Folded Reload
	.loc	1 51 8                  // kernel.cpp:51:8
	add	x28, sp, #80            // =80
	.loc	1 39 10                 // kernel.cpp:39:10
	mov	z6.d, z26.d
	fsub	z6.s, p1/m, z6.s, z29.s
	.loc	1 40 96                 // kernel.cpp:40:96
	mov	z29.d, z4.d
	fadd	z29.s, p1/m, z29.s, z25.s
	.loc	1 40 66 is_stmt 0       // kernel.cpp:40:66
	movprfx	z25.s, p1/z, z29.s
	fmla	z25.s, p1/m, z22.s, z22.s
	.loc	1 40 36                 // kernel.cpp:40:36
	movprfx	z29.s, p1/z, z25.s
	fmla	z29.s, p1/m, z21.s, z21.s
	.loc	1 40 6                  // kernel.cpp:40:6
	movprfx	z25.s, p1/z, z29.s
	fmla	z25.s, p1/m, z6.s, z6.s
	.loc	1 34 6 is_stmt 1        // kernel.cpp:34:6
	mov	z30.s, s30
	.loc	1 42 8                  // kernel.cpp:42:8
	frsqrte	z29.s, z25.s
	.loc	1 43 17                 // kernel.cpp:43:17
	movprfx	z25.s, p1/z, z25.s
	fmul	z25.s, p1/m, z25.s, z29.s
	.loc	1 44 5                  // kernel.cpp:44:5
	movprfx	z25.s, p1/z, z25.s
	fmsb	z25.s, p1/m, z29.s, z17.s
	.loc	1 45 20                 // kernel.cpp:45:20
	movprfx	z1.s, p1/z, z16.s
	fmla	z1.s, p1/m, z25.s, z20.s
	.loc	1 46 8                  // kernel.cpp:46:8
	movprfx	z1.s, p1/z, z1.s
	fmul	z1.s, p1/m, z1.s, z25.s
	.loc	1 47 8                  // kernel.cpp:47:8
	movprfx	z25.s, p1/z, z29.s
	fmla	z25.s, p1/m, z29.s, z1.s
	.loc	1 37 10                 // kernel.cpp:37:10
	mov	z1.d, z19.d
	fsub	z1.s, p1/m, z1.s, z27.s
	.loc	1 38 10                 // kernel.cpp:38:10
	mov	z27.d, z24.d
	fsub	z27.s, p1/m, z27.s, z10.s
	.loc	1 40 96                 // kernel.cpp:40:96
	mov	z10.d, z4.d
	fadd	z10.s, p1/m, z10.s, z28.s
	.loc	1 39 10                 // kernel.cpp:39:10
	mov	z29.d, z26.d
	.loc	1 49 9                  // kernel.cpp:49:9
	movprfx	z23.s, p1/z, z23.s
	fmul	z23.s, p1/m, z23.s, z25.s
	.loc	1 40 66                 // kernel.cpp:40:66
	movprfx	z28.s, p1/z, z10.s
	fmla	z28.s, p1/m, z1.s, z1.s
	.loc	1 39 10                 // kernel.cpp:39:10
	fsub	z29.s, p1/m, z29.s, z11.s
	.loc	1 40 36                 // kernel.cpp:40:36
	movprfx	z10.s, p1/z, z28.s
	fmla	z10.s, p1/m, z27.s, z27.s
	.loc	1 40 6 is_stmt 0        // kernel.cpp:40:6
	movprfx	z28.s, p1/z, z10.s
	fmla	z28.s, p1/m, z29.s, z29.s
	.loc	1 42 8 is_stmt 1        // kernel.cpp:42:8
	frsqrte	z10.s, z28.s
	.loc	1 43 17                 // kernel.cpp:43:17
	movprfx	z28.s, p1/z, z28.s
	fmul	z28.s, p1/m, z28.s, z10.s
	.loc	1 44 5                  // kernel.cpp:44:5
	movprfx	z28.s, p1/z, z28.s
	fmsb	z28.s, p1/m, z10.s, z17.s
	.loc	1 45 20                 // kernel.cpp:45:20
	movprfx	z11.s, p1/z, z16.s
	fmla	z11.s, p1/m, z28.s, z20.s
	.loc	1 46 8                  // kernel.cpp:46:8
	movprfx	z28.s, p1/z, z28.s
	fmul	z28.s, p1/m, z28.s, z11.s
	.loc	1 47 8                  // kernel.cpp:47:8
	movprfx	z11.s, p1/z, z10.s
	fmla	z11.s, p1/m, z10.s, z28.s
	.loc	1 40 96                 // kernel.cpp:40:96
	mov	z10.d, z4.d
	.loc	1 49 9                  // kernel.cpp:49:9
	movprfx	z28.s, p1/z, z30.s
	fmul	z28.s, p1/m, z28.s, z11.s
	.loc	1 37 10                 // kernel.cpp:37:10
	mov	z30.d, z19.d
	.loc	1 40 96                 // kernel.cpp:40:96
	fadd	z10.s, p1/m, z10.s, z13.s
	.loc	1 37 10                 // kernel.cpp:37:10
	fsub	z30.s, p1/m, z30.s, z8.s
	.loc	1 38 10                 // kernel.cpp:38:10
	mov	z8.d, z24.d
	fsub	z8.s, p1/m, z8.s, z9.s
	.loc	1 39 10                 // kernel.cpp:39:10
	mov	z9.d, z26.d
	fsub	z9.s, p1/m, z9.s, z12.s
	.loc	1 40 66                 // kernel.cpp:40:66
	movprfx	z12.s, p1/z, z10.s
	fmla	z12.s, p1/m, z30.s, z30.s
	.loc	1 40 36 is_stmt 0       // kernel.cpp:40:36
	movprfx	z10.s, p1/z, z12.s
	fmla	z10.s, p1/m, z8.s, z8.s
	.loc	1 40 6                  // kernel.cpp:40:6
	movprfx	z12.s, p1/z, z10.s
	fmla	z12.s, p1/m, z9.s, z9.s
	.loc	1 50 30 is_stmt 1       // kernel.cpp:50:30
	movprfx	z10.s, p1/z, z23.s
	fmul	z10.s, p1/m, z10.s, z25.s
	.loc	1 42 8                  // kernel.cpp:42:8
	frsqrte	z13.s, z12.s
	.loc	1 43 17                 // kernel.cpp:43:17
	movprfx	z12.s, p1/z, z12.s
	fmul	z12.s, p1/m, z12.s, z13.s
	.loc	1 44 5                  // kernel.cpp:44:5
	movprfx	z12.s, p1/z, z12.s
	fmsb	z12.s, p1/m, z13.s, z17.s
	.loc	1 45 20                 // kernel.cpp:45:20
	movprfx	z2.s, p1/z, z16.s
	fmla	z2.s, p1/m, z12.s, z20.s
	.loc	1 46 8                  // kernel.cpp:46:8
	movprfx	z2.s, p1/z, z2.s
	fmul	z2.s, p1/m, z2.s, z12.s
	.loc	1 47 8                  // kernel.cpp:47:8
	movprfx	z12.s, p1/z, z13.s
	fmla	z12.s, p1/m, z13.s, z2.s
	.loc	1 31 19                 // kernel.cpp:31:19
	fcvt	s2, d15
	.loc	1 31 9 is_stmt 0        // kernel.cpp:31:9
	mov	z2.s, s2
	.loc	1 50 30 is_stmt 1       // kernel.cpp:50:30
	movprfx	z13.s, p1/z, z28.s
	fmul	z13.s, p1/m, z13.s, z11.s
	.loc	1 32 30                 // kernel.cpp:32:30
	ldr	d15, [x18, #48]
	.loc	1 32 19 is_stmt 0       // kernel.cpp:32:19
	fcvt	s15, d15
	.loc	1 38 10 is_stmt 1       // kernel.cpp:38:10
	fsub	z24.s, p1/m, z24.s, z2.s
	.loc	1 32 9                  // kernel.cpp:32:9
	mov	z15.s, s15
	.loc	1 34 6                  // kernel.cpp:34:6
	ld1rw	{ z3.s }, p0/z, [x18, #56]
	.loc	1 36 9                  // kernel.cpp:36:9
	ld1rw	{ z0.s }, p0/z, [x18, #60]
	.loc	1 50 14                 // kernel.cpp:50:14
	movprfx	z25.s, p1/z, z25.s
	fmul	z25.s, p1/m, z25.s, z10.s
	.loc	1 37 10                 // kernel.cpp:37:10
	mov	z10.d, z19.d
	fsub	z10.s, p1/m, z10.s, z14.s
	.loc	1 28 1                  // kernel.cpp:28:1
	add	x18, x18, #128          // =128
	.loc	1 39 10                 // kernel.cpp:39:10
	fsub	z26.s, p1/m, z26.s, z15.s
	.loc	1 40 96                 // kernel.cpp:40:96
	mov	z15.d, z4.d
	fadd	z15.s, p1/m, z15.s, z0.s
	.loc	1 40 66 is_stmt 0       // kernel.cpp:40:66
	movprfx	z0.s, p1/z, z15.s
	fmla	z0.s, p1/m, z10.s, z10.s
	.loc	1 40 36                 // kernel.cpp:40:36
	movprfx	z15.s, p1/z, z0.s
	fmla	z15.s, p1/m, z24.s, z24.s
	.loc	1 40 6                  // kernel.cpp:40:6
	movprfx	z0.s, p1/z, z15.s
	fmla	z0.s, p1/m, z26.s, z26.s
	.loc	1 50 14 is_stmt 1       // kernel.cpp:50:14
	movprfx	z11.s, p1/z, z11.s
	fmul	z11.s, p1/m, z11.s, z13.s
	.loc	1 42 8                  // kernel.cpp:42:8
	frsqrte	z13.s, z0.s
	.loc	1 43 17                 // kernel.cpp:43:17
	movprfx	z0.s, p1/z, z0.s
	fmul	z0.s, p1/m, z0.s, z13.s
	.loc	1 44 5                  // kernel.cpp:44:5
	movprfx	z0.s, p1/z, z0.s
	fmsb	z0.s, p1/m, z13.s, z17.s
	.loc	1 45 20                 // kernel.cpp:45:20
	movprfx	z15.s, p1/z, z16.s
	fmla	z15.s, p1/m, z0.s, z20.s
	.loc	1 46 8                  // kernel.cpp:46:8
	movprfx	z0.s, p1/z, z0.s
	fmul	z0.s, p1/m, z0.s, z15.s
	.loc	1 49 9                  // kernel.cpp:49:9
	movprfx	z31.s, p1/z, z31.s
	fmul	z31.s, p1/m, z31.s, z12.s
	.loc	1 47 8                  // kernel.cpp:47:8
	movprfx	z15.s, p1/z, z13.s
	fmla	z15.s, p1/m, z13.s, z0.s
	.loc	1 50 30                 // kernel.cpp:50:30
	movprfx	z0.s, p1/z, z31.s
	fmul	z0.s, p1/m, z0.s, z12.s
	.loc	1 50 14 is_stmt 0       // kernel.cpp:50:14
	movprfx	z0.s, p1/z, z0.s
	fmul	z0.s, p1/m, z0.s, z12.s
	.loc	1 54 7 is_stmt 1        // kernel.cpp:54:7
	movprfx	z7.s, p1/z, z18.s
	fsub	z7.s, p1/m, z7.s, z23.s
	.loc	1 51 8                  // kernel.cpp:51:8
	ldr	z16, [x28, #6, mul vl]
	ldr	z17, [x28, #7, mul vl]
	ldr	z18, [x28, #8, mul vl]
	ldr	z19, [x28, #9, mul vl]
	.loc	1 52 8                  // kernel.cpp:52:8
	add	x28, sp, #80            // =80
	.loc	1 51 8                  // kernel.cpp:51:8
	movprfx	z4.s, p1/z, z16.s
	fmls	z4.s, p1/m, z25.s, z22.s
	.loc	1 52 8                  // kernel.cpp:52:8
	ldr	z5, [x28]               // 16-byte Folded Reload
	.loc	1 53 8                  // kernel.cpp:53:8
	add	x28, sp, #80            // =80
	.loc	1 52 8                  // kernel.cpp:52:8
	movprfx	z5.s, p1/z, z5.s
	fmls	z5.s, p1/m, z25.s, z21.s
	.loc	1 53 8                  // kernel.cpp:53:8
	ldr	z16, [x28, #1, mul vl]  // 16-byte Folded Reload
	movprfx	z6.s, p1/z, z6.s
	fmsb	z6.s, p1/m, z25.s, z16.s
	.loc	1 54 7                  // kernel.cpp:54:7
	movprfx	z7.s, p1/z, z7.s
	fsub	z7.s, p1/m, z7.s, z28.s
	.loc	1 51 8                  // kernel.cpp:51:8
	movprfx	z1.s, p1/z, z1.s
	fmsb	z1.s, p1/m, z11.s, z4.s
	.loc	1 52 8                  // kernel.cpp:52:8
	movprfx	z4.s, p1/z, z5.s
	fmls	z4.s, p1/m, z11.s, z27.s
	.loc	1 53 8                  // kernel.cpp:53:8
	movprfx	z5.s, p1/z, z6.s
	fmls	z5.s, p1/m, z11.s, z29.s
	.loc	1 54 7                  // kernel.cpp:54:7
	movprfx	z18.s, p1/z, z7.s
	fsub	z18.s, p1/m, z18.s, z31.s
	.loc	1 51 8                  // kernel.cpp:51:8
	movprfx	z1.s, p1/z, z1.s
	fmls	z1.s, p1/m, z0.s, z30.s
	.loc	1 52 8                  // kernel.cpp:52:8
	movprfx	z16.s, p1/z, z4.s
	fmls	z16.s, p1/m, z0.s, z8.s
	.loc	1 53 8                  // kernel.cpp:53:8
	movprfx	z0.s, p1/z, z0.s
	fmsb	z0.s, p1/m, z9.s, z5.s
	.loc	1 49 9                  // kernel.cpp:49:9
	movprfx	z3.s, p1/z, z3.s
	fmul	z3.s, p1/m, z3.s, z15.s
	.loc	1 50 30                 // kernel.cpp:50:30
	movprfx	z4.s, p1/z, z3.s
	fmul	z4.s, p1/m, z4.s, z15.s
	.loc	1 50 14 is_stmt 0       // kernel.cpp:50:14
	movprfx	z19.s, p1/z, z4.s
	fmul	z19.s, p1/m, z19.s, z15.s
	.loc	1 51 8 is_stmt 1        // kernel.cpp:51:8
	movprfx	z4.s, p1/z, z1.s
	fmls	z4.s, p1/m, z19.s, z10.s
	.loc	1 52 8                  // kernel.cpp:52:8
	movprfx	z21.s, p1/z, z16.s
	fmls	z21.s, p1/m, z19.s, z24.s
	.loc	1 53 8                  // kernel.cpp:53:8
	movprfx	z20.s, p1/z, z0.s
	fmls	z20.s, p1/m, z19.s, z26.s
	.loc	1 54 7                  // kernel.cpp:54:7
	movprfx	z18.s, p1/z, z18.s
	fsub	z18.s, p1/m, z18.s, z3.s
	.loc	1 28 1                  // kernel.cpp:28:1
	cmp	x12, x17
	b.ne	.LBB0_6
// %bb.7:                               //   in Loop: Header=BB0_3 Depth=1
	.loc	1 0 1 is_stmt 0         // kernel.cpp:0:1
	add	x18, sp, #80            // =80
	ldr	z28, [x18, #2, mul vl]  // 16-byte Folded Reload
	add	x18, sp, #80            // =80
	ldr	z29, [x18, #3, mul vl]  // 16-byte Folded Reload
	add	x18, sp, #80            // =80
	ldr	z30, [x18, #4, mul vl]  // 16-byte Folded Reload
	add	x18, sp, #80            // =80
	ldr	z25, [x18, #5, mul vl]  // 16-byte Folded Reload
	mov	z26.s, #0               // =0x0
	fmov	z16.s, #1.00000000
	fmov	z17.s, #0.37500000
	fmov	z27.s, #0.50000000
	.loc	1 28 1                  // kernel.cpp:28:1
	cbz	x11, .LBB0_10
.LBB0_8:                                //   in Loop: Header=BB0_3 Depth=1
	.loc	1 30 30 is_stmt 1       // kernel.cpp:30:30
	add	x17, x14, x17, lsl #5
	mov	x18, x15
	.p2align	2
.LBB0_9:                                //   Parent Loop BB0_3 Depth=1
                                        // =>  This Inner Loop Header: Depth=2
	ldp	d0, d1, [x17, #-16]
	.loc	1 30 19 is_stmt 0       // kernel.cpp:30:19
	fcvt	s0, d0
	.loc	1 31 19 is_stmt 1       // kernel.cpp:31:19
	fcvt	s1, d1
	.loc	1 30 9                  // kernel.cpp:30:9
	mov	z0.s, s0
	.loc	1 31 9                  // kernel.cpp:31:9
	mov	z1.s, s1
	.loc	1 32 30                 // kernel.cpp:32:30
	ldr	d2, [x17]
	.loc	1 34 6                  // kernel.cpp:34:6
	ld1rw	{ z3.s }, p0/z, [x17, #8]
	.loc	1 36 9                  // kernel.cpp:36:9
	ld1rw	{ z19.s }, p0/z, [x17, #12]
	.loc	1 28 1                  // kernel.cpp:28:1
	add	x17, x17, #32           // =32
	.loc	1 32 19                 // kernel.cpp:32:19
	fcvt	s2, d2
	.loc	1 37 10                 // kernel.cpp:37:10
	mov	z23.d, z28.d
	fsub	z23.s, p1/m, z23.s, z0.s
	.loc	1 32 9                  // kernel.cpp:32:9
	mov	z0.s, s2
	.loc	1 38 10                 // kernel.cpp:38:10
	mov	z2.d, z29.d
	fsub	z2.s, p1/m, z2.s, z1.s
	.loc	1 39 10                 // kernel.cpp:39:10
	mov	z1.d, z30.d
	fsub	z1.s, p1/m, z1.s, z0.s
	.loc	1 40 96                 // kernel.cpp:40:96
	mov	z0.d, z25.d
	fadd	z0.s, p1/m, z0.s, z19.s
	.loc	1 40 66 is_stmt 0       // kernel.cpp:40:66
	movprfx	z19.s, p1/z, z0.s
	fmla	z19.s, p1/m, z23.s, z23.s
	.loc	1 40 36                 // kernel.cpp:40:36
	movprfx	z0.s, p1/z, z19.s
	fmla	z0.s, p1/m, z2.s, z2.s
	.loc	1 40 6                  // kernel.cpp:40:6
	movprfx	z19.s, p1/z, z0.s
	fmla	z19.s, p1/m, z1.s, z1.s
	.loc	1 42 8 is_stmt 1        // kernel.cpp:42:8
	frsqrte	z0.s, z19.s
	.loc	1 43 17                 // kernel.cpp:43:17
	movprfx	z19.s, p1/z, z19.s
	fmul	z19.s, p1/m, z19.s, z0.s
	.loc	1 44 5                  // kernel.cpp:44:5
	movprfx	z19.s, p1/z, z19.s
	fmsb	z19.s, p1/m, z0.s, z16.s
	.loc	1 45 20                 // kernel.cpp:45:20
	movprfx	z24.s, p1/z, z27.s
	fmla	z24.s, p1/m, z19.s, z17.s
	.loc	1 46 8                  // kernel.cpp:46:8
	movprfx	z19.s, p1/z, z19.s
	fmul	z19.s, p1/m, z19.s, z24.s
	.loc	1 47 8                  // kernel.cpp:47:8
	movprfx	z24.s, p1/z, z0.s
	fmla	z24.s, p1/m, z0.s, z19.s
	.loc	1 49 9                  // kernel.cpp:49:9
	movprfx	z0.s, p1/z, z3.s
	fmul	z0.s, p1/m, z0.s, z24.s
	.loc	1 50 30                 // kernel.cpp:50:30
	movprfx	z3.s, p1/z, z0.s
	fmul	z3.s, p1/m, z3.s, z24.s
	.loc	1 50 14 is_stmt 0       // kernel.cpp:50:14
	movprfx	z3.s, p1/z, z3.s
	fmul	z3.s, p1/m, z3.s, z24.s
	.loc	1 51 8 is_stmt 1        // kernel.cpp:51:8
	movprfx	z4.s, p1/z, z4.s
	fmls	z4.s, p1/m, z3.s, z23.s
	.loc	1 52 8                  // kernel.cpp:52:8
	movprfx	z21.s, p1/z, z21.s
	fmls	z21.s, p1/m, z3.s, z2.s
	.loc	1 53 8                  // kernel.cpp:53:8
	movprfx	z20.s, p1/z, z20.s
	fmls	z20.s, p1/m, z3.s, z1.s
	.loc	1 54 7                  // kernel.cpp:54:7
	movprfx	z18.s, p1/z, z18.s
	fsub	z18.s, p1/m, z18.s, z0.s
	.loc	1 28 1                  // kernel.cpp:28:1
	adds	x18, x18, #1            // =1
	b.ne	.LBB0_9
.LBB0_10:                               //   in Loop: Header=BB0_3 Depth=1
	.loc	1 60 1                  // kernel.cpp:60:1
	mov	z5.d, z21.d
	mov	z6.d, z20.d
	mov	z7.d, z18.d
	.loc	1 12 32                 // kernel.cpp:12:32
	add	x9, x9, #16             // =16
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z4.s, z5.s, z6.s, z7.s }, p1, [x16]
	.loc	1 12 1                  // kernel.cpp:12:1
	cmp	x9, x8
	b.lt	.LBB0_3
	b	.LBB0_18
.LBB0_11:
	.loc	1 13 7                  // kernel.cpp:13:7
	cmp	x8, #16                 // =16
	orr	w9, wzr, #0x10
	csel	x8, x8, x9, gt
	sub	x8, x8, #1              // =1
	lsr	x11, x8, #4
	add	w9, w11, #1             // =1
	and	x9, x9, #0x7
	cmp	x8, #112                // =112
	b.hs	.LBB0_13
// %bb.12:
	.loc	1 0 7 is_stmt 0         // kernel.cpp:0:7
	mov	x8, xzr
	cbnz	x9, .LBB0_16
	b	.LBB0_18
.LBB0_13:
	.loc	1 13 7                  // kernel.cpp:13:7
	mvn	x13, x11
	mov	x8, xzr
	add	x10, x4, #1792          // =1792
	mov	x11, #-448
	mov	x12, #-384
	mov	x14, #-320
	add	x13, x9, x13
	orr	x15, xzr, #0xffffffffffffff00
	mov	x16, #-192
	mov	w17, #80
	orr	x18, xzr, #0xffffffffffffff80
	orr	x0, xzr, #0xffffffffffffffc0
	.p2align	2
.LBB0_14:                               // =>This Inner Loop Header: Depth=1
	whilelt	p0.s, w8, w1
	.loc	1 20 33 is_stmt 1       // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x11, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	orr	w2, w8, #0x10
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x11, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, #0x20
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x12, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x12, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, #0x30
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x14, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x14, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, #0x40
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x15, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x15, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, w17
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x16, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x16, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, #0x60
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x18, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x18, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	orr	w2, w8, #0x70
	.loc	1 12 1                  // kernel.cpp:12:1
	add	x8, x8, #128            // =128
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10, x0, lsl #2]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10, x0, lsl #2]
	.loc	1 13 7                  // kernel.cpp:13:7
	whilelt	p0.s, w2, w1
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10]
	.loc	1 12 1                  // kernel.cpp:12:1
	add	x10, x10, #2048         // =2048
	adds	x13, x13, #8            // =8
	b.ne	.LBB0_14
// %bb.15:
	.loc	1 0 1 is_stmt 0         // kernel.cpp:0:1
	cbz	x9, .LBB0_18
.LBB0_16:
	.loc	1 13 7 is_stmt 1        // kernel.cpp:13:7
	add	x10, x4, x8, lsl #4
	neg	x9, x9
	.p2align	2
.LBB0_17:                               // =>This Inner Loop Header: Depth=1
	whilelt	p0.s, w8, w1
	.loc	1 12 1                  // kernel.cpp:12:1
	add	w8, w8, #16             // =16
	.loc	1 20 33                 // kernel.cpp:20:33
	ld4w	{ z0.s, z1.s, z2.s, z3.s }, p0/z, [x10]
	.loc	1 60 1                  // kernel.cpp:60:1
	st4w	{ z0.s, z1.s, z2.s, z3.s }, p0, [x10]
	.loc	1 12 1                  // kernel.cpp:12:1
	add	x10, x10, #256          // =256
	adds	x9, x9, #1              // =1
	b.ne	.LBB0_17
.LBB0_18:
	.loc	1 62 1                  // kernel.cpp:62:1
	ldp	d9, d8, [sp, #48]       // 16-byte Folded Reload
	ldp	d11, d10, [sp, #32]     // 16-byte Folded Reload
	ldp	d13, d12, [sp, #16]     // 16-byte Folded Reload
	ldr	x28, [sp, #64]          // 8-byte Folded Reload
	ldp	d15, d14, [sp], #80     // 16-byte Folded Reload
	addvl	sp, sp, #10
	ret
.Ltmp1:
.Lfunc_end0:
	.size	_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force, .Lfunc_end0-_Z17CalcForceLongEpEpPK3EPIiPK3EPJiP5Force
	.cfi_endproc
                                        // -- End function
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4               // -- Begin function _GLOBAL__sub_I_kernel.cpp
.LCPI1_0:
	.word	3204448256              // float -0.5
	.word	1056964608              // float 0.5
	.word	3204448256              // float -0.5
	.word	1056964608              // float 0.5
.LCPI1_1:
	.word	3204448256              // float -0.5
	.word	3204448256              // float -0.5
	.word	1056964608              // float 0.5
	.word	1056964608              // float 0.5
.LCPI1_2:
	.word	1056964608              // float 0.5
	.word	3204448256              // float -0.5
	.word	3204448256              // float -0.5
	.word	1056964608              // float 0.5
.LCPI1_3:
	.word	3204448256              // float -0.5
	.word	1056964608              // float 0.5
	.word	1056964608              // float 0.5
	.word	1056964608              // float 0.5
	.section	.text.startup,"ax",@progbits
	.p2align	3
	.type	_GLOBAL__sub_I_kernel.cpp,@function
_GLOBAL__sub_I_kernel.cpp:              // @_GLOBAL__sub_I_kernel.cpp
.Lfunc_begin1:
	.loc	1 0 0                   // kernel.cpp:0:0
	.cfi_startproc
// %bb.0:
	stp	x19, x30, [sp, #-16]!   // 16-byte Folded Spill
	.cfi_def_cfa_offset 16
	.cfi_offset w30, -8
	.cfi_offset w19, -16
.Ltmp2:
	.file	2 "/home/users/go19/go0011/kernel_generator/bench/nbody" "/opt/FJSVxos/devkit/aarch64/rfs/lib/gcc/aarch64-redhat-linux/8/../../../../include/c++/8/iostream"
	.loc	2 74 25 prologue_end    // /opt/FJSVxos/devkit/aarch64/rfs/lib/gcc/aarch64-redhat-linux/8/../../../../include/c++/8/iostream:74:25
	adrp	x19, .L_MergedGlobals
	add	x19, x19, :lo12:.L_MergedGlobals
	mov	x0, x19
	bl	_ZNSt8ios_base4InitC1Ev
	adrp	x0, _ZNSt8ios_base4InitD1Ev
	adrp	x2, __dso_handle
	add	x0, x0, :lo12:_ZNSt8ios_base4InitD1Ev
	add	x2, x2, :lo12:__dso_handle
	mov	x1, x19
	bl	__cxa_atexit
	adrp	x8, .LCPI1_0
.Ltmp3:
	.file	3 "/home/users/go19/go0011/kernel_generator/bench/nbody" "/home/users/go19/go0011/FDPS/src/vector3.hpp"
	.loc	3 12 62                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:62
	ldr	q1, [x8, :lo12:.LCPI1_0]
.Ltmp4:
	.loc	3 12 55 is_stmt 0       // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:55
	movi	v0.4s, #191, lsl #24
	adrp	x9, .LCPI1_1
	adrp	x10, .LCPI1_2
	adrp	x8, .LCPI1_3
.Ltmp5:
	.loc	3 12 69                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:69
	ldr	q2, [x9, :lo12:.LCPI1_1]
.Ltmp6:
	.loc	3 12 55                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:55
	ldr	q3, [x10, :lo12:.LCPI1_2]
.Ltmp7:
	.loc	3 12 55                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:55
	stp	q0, q1, [x19, #16]
.Ltmp8:
	.loc	3 12 62                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:62
	ldr	q0, [x8, :lo12:.LCPI1_3]
.Ltmp9:
	.loc	3 12 69                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:69
	stp	q2, q3, [x19, #48]
.Ltmp10:
	.loc	3 12 62                 // /home/users/go19/go0011/FDPS/src/vector3.hpp:12:62
	stp	q0, q0, [x19, #80]
	ldp	x19, x30, [sp], #16     // 16-byte Folded Reload
	ret
.Ltmp11:
.Lfunc_end1:
	.size	_GLOBAL__sub_I_kernel.cpp, .Lfunc_end1-_GLOBAL__sub_I_kernel.cpp
	.cfi_endproc
	.file	4 "/home/users/go19/go0011/kernel_generator/bench/nbody" "/home/users/go19/go0011/FDPS/src/ps_defs.hpp"
                                        // -- End function
	.hidden	__dso_handle
	.section	.init_array,"aw",@init_array
	.p2align	3
	.xword	_GLOBAL__sub_I_kernel.cpp
	.type	.L_MergedGlobals,@object // @_MergedGlobals
	.local	.L_MergedGlobals
	.comm	.L_MergedGlobals,112,16
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	"clang: Fujitsu C/C++ Compiler 4.1.0 (Dec 23 2019 15:40:15) (based on LLVM 7.1.0)" // string offset=0
.Linfo_string1:
	.asciz	"kernel.cpp"            // string offset=81
.Linfo_string2:
	.asciz	"/home/users/go19/go0011/kernel_generator/bench/nbody" // string offset=92
.Linfo_string3:
	.asciz	"__cxx_global_var_init" // string offset=145
.Linfo_string4:
	.asciz	"Vector3"               // string offset=167
.Linfo_string5:
	.asciz	"__cxx_global_var_init.1" // string offset=175
.Linfo_string6:
	.asciz	"CalcForceLongEpEp"     // string offset=199
	.section	.debug_abbrev,"",@progbits
	.byte	1                       // Abbreviation Code
	.byte	17                      // DW_TAG_compile_unit
	.byte	1                       // DW_CHILDREN_yes
	.byte	37                      // DW_AT_producer
	.byte	14                      // DW_FORM_strp
	.byte	19                      // DW_AT_language
	.byte	5                       // DW_FORM_data2
	.byte	3                       // DW_AT_name
	.byte	14                      // DW_FORM_strp
	.byte	16                      // DW_AT_stmt_list
	.byte	23                      // DW_FORM_sec_offset
	.byte	27                      // DW_AT_comp_dir
	.byte	14                      // DW_FORM_strp
	.byte	17                      // DW_AT_low_pc
	.byte	1                       // DW_FORM_addr
	.byte	85                      // DW_AT_ranges
	.byte	23                      // DW_FORM_sec_offset
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	2                       // Abbreviation Code
	.byte	46                      // DW_TAG_subprogram
	.byte	1                       // DW_CHILDREN_yes
	.byte	17                      // DW_AT_low_pc
	.byte	1                       // DW_FORM_addr
	.byte	18                      // DW_AT_high_pc
	.byte	6                       // DW_FORM_data4
	.byte	3                       // DW_AT_name
	.byte	14                      // DW_FORM_strp
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	3                       // Abbreviation Code
	.ascii	"\200\340\003"          // DW_TAG_FJ_loop
	.byte	0                       // DW_CHILDREN_no
	.byte	58                      // DW_AT_decl_file
	.byte	11                      // DW_FORM_data1
	.ascii	"\200f"                 // DW_AT_FJ_loop_start_line
	.byte	11                      // DW_FORM_data1
	.ascii	"\201f"                 // DW_AT_FJ_loop_end_line
	.byte	11                      // DW_FORM_data1
	.ascii	"\202f"                 // DW_AT_FJ_loop_nest_level
	.byte	11                      // DW_FORM_data1
	.ascii	"\203f"                 // DW_AT_FJ_loop_type
	.byte	11                      // DW_FORM_data1
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	4                       // Abbreviation Code
	.byte	46                      // DW_TAG_subprogram
	.byte	1                       // DW_CHILDREN_yes
	.byte	17                      // DW_AT_low_pc
	.byte	1                       // DW_FORM_addr
	.byte	18                      // DW_AT_high_pc
	.byte	6                       // DW_FORM_data4
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	5                       // Abbreviation Code
	.byte	29                      // DW_TAG_inlined_subroutine
	.byte	0                       // DW_CHILDREN_no
	.byte	49                      // DW_AT_abstract_origin
	.byte	19                      // DW_FORM_ref4
	.byte	17                      // DW_AT_low_pc
	.byte	1                       // DW_FORM_addr
	.byte	18                      // DW_AT_high_pc
	.byte	6                       // DW_FORM_data4
	.byte	88                      // DW_AT_call_file
	.byte	11                      // DW_FORM_data1
	.byte	89                      // DW_AT_call_line
	.byte	11                      // DW_FORM_data1
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	6                       // Abbreviation Code
	.byte	29                      // DW_TAG_inlined_subroutine
	.byte	1                       // DW_CHILDREN_yes
	.byte	49                      // DW_AT_abstract_origin
	.byte	19                      // DW_FORM_ref4
	.byte	17                      // DW_AT_low_pc
	.byte	1                       // DW_FORM_addr
	.byte	18                      // DW_AT_high_pc
	.byte	6                       // DW_FORM_data4
	.byte	88                      // DW_AT_call_file
	.byte	11                      // DW_FORM_data1
	.byte	89                      // DW_AT_call_line
	.byte	11                      // DW_FORM_data1
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	7                       // Abbreviation Code
	.byte	29                      // DW_TAG_inlined_subroutine
	.byte	0                       // DW_CHILDREN_no
	.byte	49                      // DW_AT_abstract_origin
	.byte	19                      // DW_FORM_ref4
	.byte	85                      // DW_AT_ranges
	.byte	23                      // DW_FORM_sec_offset
	.byte	88                      // DW_AT_call_file
	.byte	11                      // DW_FORM_data1
	.byte	89                      // DW_AT_call_line
	.byte	11                      // DW_FORM_data1
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	8                       // Abbreviation Code
	.byte	46                      // DW_TAG_subprogram
	.byte	0                       // DW_CHILDREN_no
	.byte	3                       // DW_AT_name
	.byte	14                      // DW_FORM_strp
	.byte	0                       // EOM(1)
	.byte	0                       // EOM(2)
	.byte	0                       // EOM(3)
	.section	.debug_info,"",@progbits
.Lcu_begin0:
	.word	208                     // Length of Unit
	.hword	4                       // DWARF version number
	.word	.debug_abbrev           // Offset Into Abbrev. Section
	.byte	8                       // Address Size (in bytes)
	.byte	1                       // Abbrev [1] 0xb:0xc9 DW_TAG_compile_unit
	.word	.Linfo_string0          // DW_AT_producer
	.hword	4                       // DW_AT_language
	.word	.Linfo_string1          // DW_AT_name
	.word	.Lline_table_start0     // DW_AT_stmt_list
	.word	.Linfo_string2          // DW_AT_comp_dir
	.xword	0                       // DW_AT_low_pc
	.word	.Ldebug_ranges3         // DW_AT_ranges
	.byte	2                       // Abbrev [2] 0x2a:0x1e DW_TAG_subprogram
	.xword	.Lfunc_begin0           // DW_AT_low_pc
	.word	.Lfunc_end0-.Lfunc_begin0 // DW_AT_high_pc
	.word	.Linfo_string6          // DW_AT_name
	.byte	3                       // Abbrev [3] 0x3b:0x6 DW_TAG_FJ_loop
	.byte	1                       // DW_AT_decl_file
	.byte	12                      // DW_AT_FJ_loop_start_line
	.byte	61                      // DW_AT_FJ_loop_end_line
	.byte	1                       // DW_AT_FJ_loop_nest_level
	.byte	5                       // DW_AT_FJ_loop_type
	.byte	3                       // Abbrev [3] 0x41:0x6 DW_TAG_FJ_loop
	.byte	1                       // DW_AT_decl_file
	.byte	28                      // DW_AT_FJ_loop_start_line
	.byte	55                      // DW_AT_FJ_loop_end_line
	.byte	2                       // DW_AT_FJ_loop_nest_level
	.byte	5                       // DW_AT_FJ_loop_type
	.byte	0                       // End Of Children Mark
	.byte	4                       // Abbrev [4] 0x48:0x7c DW_TAG_subprogram
	.xword	.Lfunc_begin1           // DW_AT_low_pc
	.word	.Lfunc_end1-.Lfunc_begin1 // DW_AT_high_pc
	.byte	5                       // Abbrev [5] 0x55:0x13 DW_TAG_inlined_subroutine
	.word	196                     // DW_AT_abstract_origin
	.xword	.Ltmp2                  // DW_AT_low_pc
	.word	.Ltmp3-.Ltmp2           // DW_AT_high_pc
	.byte	1                       // DW_AT_call_file
	.byte	0                       // DW_AT_call_line
	.byte	6                       // Abbrev [6] 0x68:0x5b DW_TAG_inlined_subroutine
	.word	206                     // DW_AT_abstract_origin
	.xword	.Ltmp3                  // DW_AT_low_pc
	.word	.Ltmp11-.Ltmp3          // DW_AT_high_pc
	.byte	1                       // DW_AT_call_file
	.byte	0                       // DW_AT_call_line
	.byte	5                       // Abbrev [5] 0x7b:0x13 DW_TAG_inlined_subroutine
	.word	201                     // DW_AT_abstract_origin
	.xword	.Ltmp3                  // DW_AT_low_pc
	.word	.Ltmp4-.Ltmp3           // DW_AT_high_pc
	.byte	4                       // DW_AT_call_file
	.byte	153                     // DW_AT_call_line
	.byte	7                       // Abbrev [7] 0x8e:0xb DW_TAG_inlined_subroutine
	.word	201                     // DW_AT_abstract_origin
	.word	.Ldebug_ranges0         // DW_AT_ranges
	.byte	4                       // DW_AT_call_file
	.byte	153                     // DW_AT_call_line
	.byte	7                       // Abbrev [7] 0x99:0xb DW_TAG_inlined_subroutine
	.word	201                     // DW_AT_abstract_origin
	.word	.Ldebug_ranges1         // DW_AT_ranges
	.byte	4                       // DW_AT_call_file
	.byte	154                     // DW_AT_call_line
	.byte	5                       // Abbrev [5] 0xa4:0x13 DW_TAG_inlined_subroutine
	.word	201                     // DW_AT_abstract_origin
	.xword	.Ltmp6                  // DW_AT_low_pc
	.word	.Ltmp7-.Ltmp6           // DW_AT_high_pc
	.byte	4                       // DW_AT_call_file
	.byte	155                     // DW_AT_call_line
	.byte	7                       // Abbrev [7] 0xb7:0xb DW_TAG_inlined_subroutine
	.word	201                     // DW_AT_abstract_origin
	.word	.Ldebug_ranges2         // DW_AT_ranges
	.byte	4                       // DW_AT_call_file
	.byte	155                     // DW_AT_call_line
	.byte	0                       // End Of Children Mark
	.byte	0                       // End Of Children Mark
	.byte	8                       // Abbrev [8] 0xc4:0x5 DW_TAG_subprogram
	.word	.Linfo_string3          // DW_AT_name
	.byte	8                       // Abbrev [8] 0xc9:0x5 DW_TAG_subprogram
	.word	.Linfo_string4          // DW_AT_name
	.byte	8                       // Abbrev [8] 0xce:0x5 DW_TAG_subprogram
	.word	.Linfo_string5          // DW_AT_name
	.byte	0                       // End Of Children Mark
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.xword	.Ltmp4
	.xword	.Ltmp5
	.xword	.Ltmp7
	.xword	.Ltmp8
	.xword	0
	.xword	0
.Ldebug_ranges1:
	.xword	.Ltmp5
	.xword	.Ltmp6
	.xword	.Ltmp9
	.xword	.Ltmp10
	.xword	0
	.xword	0
.Ldebug_ranges2:
	.xword	.Ltmp8
	.xword	.Ltmp9
	.xword	.Ltmp10
	.xword	.Ltmp11
	.xword	0
	.xword	0
.Ldebug_ranges3:
	.xword	.Lfunc_begin0
	.xword	.Lfunc_end0
	.xword	.Lfunc_begin1
	.xword	.Lfunc_end1
	.xword	0
	.xword	0
	.section	.debug_macinfo,"",@progbits
	.byte	0                       // End Of Macro List Mark

.set _ZStL8__ioinit, .L_MergedGlobals
	.size	_ZStL8__ioinit, 1
.set _ZN17ParticleSimulatorL12SHIFT_CENTERE, .L_MergedGlobals+16
	.size	_ZN17ParticleSimulatorL12SHIFT_CENTERE, 96
	.ident	"clang: Fujitsu C/C++ Compiler 4.1.0 (Dec 23 2019 15:40:15) (based on LLVM 7.1.0)"
	.section	.fj.compile_info, "e"
	.ascii	"C++::clang-libstdc++"
	.section	".note.GNU-stack","",@progbits
	.section	.debug_line,"",@progbits
.Lline_table_start0:
