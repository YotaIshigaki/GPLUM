	.ident	"$Options: Fujitsu C/C++ Compiler Version 1.2.0 P-id: L30000-09 (Sep 17 2014 11:41:29) --preinclude /opt/FJSVtclang/GM-1.2.0-17/bin/../lib/FCC.pre --g++ -D__FUJITSU -Dunix -Dsparc -D__sparc__ -D__unix -D__sparc -D__BUILTIN_VA_ARG_INCR -D__PRAGMA_REDEFINE_EXTNAME -D__FCC_VERSION=700 -D__USER_LABEL_PREFIX__= -D__OPTIMIZE__ -D__HPC_ACE__ -D__ELF__ -D__linux -Asystem(unix) -Dlinux -D__LIBC_6B -D_LP64 -D__LP64__ -D__XOSDEVKIT_PATH__=/opt/FJSVXosDevkit/sparc64fx/target --K=noocl --lp --zmode=64 --sys_include=/opt/FJSVtclang/GM-1.2.0-17/bin/../include/c++/std --sys_include=/opt/FJSVtclang/GM-1.2.0-17/bin/../include/c++ --sys_include=/opt/FJSVtclang/GM-1.2.0-17/bin/../include --sys_include=/opt/FJSVXosDevkit/sparc64fx/target/usr/include --K=opt -D__sparcv9 -D__sparc_v9__ -D__arch64__ --exceptions compile.C -- -ncmdname=FCCpx -Nnoline -Kdalign -zobe=no-static-clump -zobe=cplus -zcfc=8fx -O3 -x- -Kdalign,ns,mfunc,lib,eval,rdconv,prefetch_conditional,fp_contract,fp_relaxed,ilfunc,fast_matmul,omitfp,noswp,nounroll -KLP compile.s $"
	.file	"compile.C"
	.ident	"$Compiler: Fujitsu C/C++ Compiler Version 1.2.0 P-id: L30000-09 (Sep 17 2014 11:41:29) compile.C _ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE $"
	.section	".text"
	.global	_ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE
	.align	64
_ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE:
.LLFB1:
.L5557:

/*     61 */	add	%sp,-2048,%sp


.L5558:

/*     68 */	sub	%o2,1,%g1

/*     68 */	cmp	%o2,%g0

/*     68 */	srl	%g1,31,%g2

/*     68 */	add	%g1,%g2,%g1

/*     68 */	sra	%g1,1,%g1


/*     68 */	ble	.L5570
/*     68 */	add	%g1,1,%g1


.L5559:

/*     68 */	sethi	%h44(.LR0.cnt.7),%g3

/*     68 */	sxar1
/*     68 */	fzero,s	%f180

/*     68 */	sethi	%h44(.LR0.cnt.2),%g4

/*     68 */	sxar1
/*     68 */	fzero	%f64

/*     68 */	or	%g3,%m44(.LR0.cnt.7),%g3

/*     68 */	or	%g4,%m44(.LR0.cnt.2),%g4

/*     68 */	sllx	%g3,12,%g3

/*     68 */	sllx	%g4,12,%g4

/*     68 */	or	%g3,%l44(.LR0.cnt.7),%g3

/*     68 */	or	%g4,%l44(.LR0.cnt.2),%g4

/*     68 */	sxar1
/*     68 */	ldd	[%g3],%f296

/*     68 */	ldd	[%g4],%f32

/*     68 */	sethi	%h44(.LS0.cnt.43),%o0

/*     68 */	sxar1
/*     68 */	sethi	%h44(.LR0.cnt.5),%xg0

/*     68 */	or	%o0,%m44(.LS0.cnt.43),%o0

/*     68 */	sxar1
/*     68 */	sethi	%h44(.LR0.cnt.10),%xg2

/*     68 */	ldd	[%g3],%f40


/*     68 */	sethi	%h44(.LR0.cnt.13),%g5

/*     68 */	sllx	%o0,12,%o0

/*     68 */	or	%o0,%l44(.LS0.cnt.43),%o0

/*     68 */	sxar1
/*     68 */	or	%xg0,%m44(.LR0.cnt.5),%xg0

/*     68 */	or	%g5,%m44(.LR0.cnt.13),%g5

/*     68 */	ldd	[%o0],%f42


/*     68 */	sxar2
/*     68 */	or	%xg2,%m44(.LR0.cnt.10),%xg2
/*     68 */	sllx	%xg0,12,%xg0

/*     68 */	sxar1
/*     68 */	sllx	%xg2,12,%xg2

/*     68 */	fabsd	%f32,%f34

/*     68 */	sllx	%g5,12,%g5


/*     68 */	sxar2
/*     68 */	or	%xg0,%l44(.LR0.cnt.5),%xg0
/*     68 */	sethi	%h44(.LR0.cnt.8),%xg1


/*     68 */	sxar2
/*     68 */	or	%xg2,%l44(.LR0.cnt.10),%xg2
/*     68 */	ldd	[%xg0],%f176


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LR0.cnt.14),%xg3
/*     68 */	sethi	%h44(.LR0.cnt.15),%xg4

/*     68 */	sxar1
/*     68 */	ldd	[%xg0],%f432


/*     68 */	or	%g5,%l44(.LR0.cnt.13),%g5


/*     68 */	sxar2
/*     68 */	sethi	%h44(_ZN6kernel14normalzeFactorE),%xg5
/*     68 */	ldd	[%xg2],%f186


/*     68 */	sxar2
/*     68 */	or	%xg1,%m44(.LR0.cnt.8),%xg1
/*     68 */	or	%xg3,%m44(.LR0.cnt.14),%xg3

/*     68 */	sxar1
/*     68 */	ldd	[%xg2],%f442


/*     68 */	ldd	[%g5],%f44


/*     68 */	sxar2
/*     68 */	or	%xg4,%m44(.LR0.cnt.15),%xg4
/*     68 */	or	%xg5,%m44(_ZN6kernel14normalzeFactorE),%xg5

/*     68 */	fdtox	%f34,%f36


/*     68 */	sxar2
/*     68 */	ldd	[%g5],%f300
/*     68 */	sethi	%h44(.LS0.cnt.34),%xg6



/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.31),%xg7
/*     68 */	sethi	%h44(.LS0.cnt.32),%xg8


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.46),%xg9
/*     68 */	sethi	%h44(.LS0.cnt.33),%xg10


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.28),%xg11
/*     68 */	sethi	%h44(.LS0.cnt.27),%xg12


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.26),%xg13
/*     68 */	sethi	%h44(.LS0.cnt.25),%xg14


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.24),%xg15
/*     68 */	sethi	%h44(.LS0.cnt.44),%xg16

/*     68 */	sxar1
/*     68 */	sethi	%h44(.LS0.cnt.36),%xg17

/*     68 */	fxtod	%f36,%f38


/*     68 */	sxar2
/*     68 */	fpmaddx	%f36,%f42,%f296,%f36
/*     68 */	sethi	%h44(.LR0.cnt.11),%xg18


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.45),%xg19
/*     68 */	sethi	%h44(.LS0.cnt.29),%xg20


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.47),%xg21
/*     68 */	sethi	%h44(.LS0.cnt.30),%xg22


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.41),%xg23
/*     68 */	sethi	%h44(.LS0.cnt.42),%xg24


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.40),%xg25
/*     68 */	sethi	%h44(.LS0.cnt.39),%xg26


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.38),%xg27
/*     68 */	sethi	%h44(.LS0.cnt.37),%xg28


/*     68 */	sxar2
/*     68 */	sethi	%h44(.LS0.cnt.35),%xg29
/*     68 */	sethi	%h44(.LS0.cnt.49),%xg30

/*     68 */	sxar1
/*     68 */	sethi	%h44(.LS0.cnt.50),%xg31

/*     68 */	sethi	%h44(.LS0.cnt.48),%g2

/*     68 */	sethi	%h44(.LR0.cnt.3),%g3


/*     68 */	sxar2
/*     68 */	sllx	%xg1,12,%xg1
/*     68 */	sllx	%xg3,12,%xg3


/*     68 */	sxar2
/*     68 */	sllx	%xg4,12,%xg4
/*     68 */	sllx	%xg5,12,%xg5


/*     68 */	sxar2
/*     68 */	or	%xg6,%m44(.LS0.cnt.34),%xg6
/*     68 */	or	%xg7,%m44(.LS0.cnt.31),%xg7


/*     68 */	sxar2
/*     68 */	or	%xg8,%m44(.LS0.cnt.32),%xg8
/*     68 */	or	%xg9,%m44(.LS0.cnt.46),%xg9


/*     68 */	sxar2
/*     68 */	or	%xg10,%m44(.LS0.cnt.33),%xg10
/*     68 */	or	%xg11,%m44(.LS0.cnt.28),%xg11


/*     68 */	sxar2
/*     68 */	or	%xg12,%m44(.LS0.cnt.27),%xg12
/*     68 */	or	%xg13,%m44(.LS0.cnt.26),%xg13


/*     68 */	sxar2
/*     68 */	or	%xg14,%m44(.LS0.cnt.25),%xg14
/*     68 */	or	%xg15,%m44(.LS0.cnt.24),%xg15


/*     68 */	sxar2
/*     68 */	or	%xg16,%m44(.LS0.cnt.44),%xg16
/*     68 */	or	%xg17,%m44(.LS0.cnt.36),%xg17


/*     68 */	sxar2
/*     68 */	or	%xg18,%m44(.LR0.cnt.11),%xg18
/*     68 */	or	%xg19,%m44(.LS0.cnt.45),%xg19


/*     68 */	sxar2
/*     68 */	or	%xg20,%m44(.LS0.cnt.29),%xg20
/*     68 */	or	%xg21,%m44(.LS0.cnt.47),%xg21


/*     68 */	sxar2
/*     68 */	or	%xg22,%m44(.LS0.cnt.30),%xg22
/*     68 */	or	%xg23,%m44(.LS0.cnt.41),%xg23


/*     68 */	sxar2
/*     68 */	or	%xg24,%m44(.LS0.cnt.42),%xg24
/*     68 */	or	%xg25,%m44(.LS0.cnt.40),%xg25


/*     68 */	sxar2
/*     68 */	or	%xg26,%m44(.LS0.cnt.39),%xg26
/*     68 */	or	%xg27,%m44(.LS0.cnt.38),%xg27


/*     68 */	sxar2
/*     68 */	or	%xg28,%m44(.LS0.cnt.37),%xg28
/*     68 */	or	%xg29,%m44(.LS0.cnt.35),%xg29


/*     68 */	sxar2
/*     68 */	or	%xg30,%m44(.LS0.cnt.49),%xg30
/*     68 */	or	%xg31,%m44(.LS0.cnt.50),%xg31

/*     68 */	or	%g2,%m44(.LS0.cnt.48),%g2

/*     68 */	or	%g3,%m44(.LR0.cnt.3),%g3


/*     68 */	sxar2
/*     68 */	or	%xg1,%l44(.LR0.cnt.8),%xg1
/*     68 */	or	%xg3,%l44(.LR0.cnt.14),%xg3


/*     68 */	sxar2
/*     68 */	or	%xg4,%l44(.LR0.cnt.15),%xg4
/*     68 */	or	%xg5,%l44(_ZN6kernel14normalzeFactorE),%xg5


/*     68 */	sxar2
/*     68 */	ldd	[%xg1],%f178
/*     68 */	sllx	%xg6,12,%xg6


/*     68 */	sxar2
/*     68 */	sllx	%xg7,12,%xg7
/*     68 */	ldd	[%xg1],%f434



/*     68 */	sxar2
/*     68 */	sllx	%xg8,12,%xg8
/*     68 */	sllx	%xg9,12,%xg9


/*     68 */	sxar2
/*     68 */	ldd	[%xg3],%f184
/*     68 */	sllx	%xg10,12,%xg10


/*     68 */	sxar2
/*     68 */	sllx	%xg11,12,%xg11
/*     68 */	ldd	[%xg3],%f440



/*     68 */	sxar2
/*     68 */	sllx	%xg12,12,%xg12
/*     68 */	ldd	[%xg4],%f182


/*     68 */	sxar2
/*     68 */	sllx	%xg13,12,%xg13
/*     68 */	sllx	%xg14,12,%xg14


/*     68 */	sxar2
/*     68 */	ldd	[%xg4],%f438
/*     68 */	sllx	%xg15,12,%xg15



/*     68 */	sxar2
/*     68 */	ldd	[%xg5],%f56
/*     68 */	sllx	%xg16,12,%xg16


/*     68 */	sxar2
/*     68 */	sllx	%xg17,12,%xg17
/*     68 */	ldd	[%xg5],%f312


/*     68 */	sxar2
/*     68 */	sllx	%xg18,12,%xg18
/*     68 */	sllx	%xg19,12,%xg19



/*     68 */	sxar2
/*     68 */	sllx	%xg20,12,%xg20
/*     68 */	sllx	%xg21,12,%xg21


/*     68 */	sxar2
/*     68 */	sllx	%xg22,12,%xg22
/*     68 */	sllx	%xg23,12,%xg23


/*     68 */	sxar2
/*     68 */	sllx	%xg24,12,%xg24
/*     68 */	sllx	%xg25,12,%xg25


/*     68 */	sxar2
/*     68 */	sllx	%xg26,12,%xg26
/*     68 */	sllx	%xg27,12,%xg27


/*     68 */	sxar2
/*     68 */	sllx	%xg28,12,%xg28
/*     68 */	sllx	%xg29,12,%xg29


/*     68 */	sxar2
/*     68 */	sllx	%xg30,12,%xg30
/*     68 */	sllx	%xg31,12,%xg31

/*     68 */	sllx	%g2,12,%g2

/*     68 */	sllx	%g3,12,%g3


/*     68 */	sxar2
/*     68 */	or	%xg6,%l44(.LS0.cnt.34),%xg6
/*     68 */	or	%xg7,%l44(.LS0.cnt.31),%xg7


/*     68 */	sxar2
/*     68 */	or	%xg8,%l44(.LS0.cnt.32),%xg8
/*     68 */	or	%xg9,%l44(.LS0.cnt.46),%xg9


/*     68 */	sxar2
/*     68 */	ldd	[%xg6],%f90
/*     68 */	or	%xg10,%l44(.LS0.cnt.33),%xg10


/*     68 */	sxar2
/*     68 */	or	%xg11,%l44(.LS0.cnt.28),%xg11
/*     68 */	ldd	[%xg7],%f86


/*     68 */	sxar2
/*     68 */	or	%xg12,%l44(.LS0.cnt.27),%xg12
/*     68 */	or	%xg13,%l44(.LS0.cnt.26),%xg13


/*     68 */	sxar2
/*     68 */	ldd	[%xg8],%f94
/*     68 */	or	%xg14,%l44(.LS0.cnt.25),%xg14


/*     68 */	sxar2
/*     68 */	or	%xg15,%l44(.LS0.cnt.24),%xg15
/*     68 */	ldd	[%xg9],%f98


/*     68 */	sxar2
/*     68 */	or	%xg16,%l44(.LS0.cnt.44),%xg16
/*     68 */	or	%xg17,%l44(.LS0.cnt.36),%xg17


/*     68 */	sxar2
/*     68 */	ldd	[%xg10],%f96
/*     68 */	or	%xg18,%l44(.LR0.cnt.11),%xg18


/*     68 */	sxar2
/*     68 */	or	%xg19,%l44(.LS0.cnt.45),%xg19
/*     68 */	ldd	[%xg11],%f100


/*     68 */	sxar2
/*     68 */	or	%xg20,%l44(.LS0.cnt.29),%xg20
/*     68 */	or	%xg21,%l44(.LS0.cnt.47),%xg21


/*     68 */	sxar2
/*     68 */	ldd	[%xg12],%f122
/*     68 */	or	%xg22,%l44(.LS0.cnt.30),%xg22


/*     68 */	sxar2
/*     68 */	or	%xg23,%l44(.LS0.cnt.41),%xg23
/*     68 */	ldd	[%xg13],%f118


/*     68 */	sxar2
/*     68 */	or	%xg24,%l44(.LS0.cnt.42),%xg24
/*     68 */	or	%xg25,%l44(.LS0.cnt.40),%xg25


/*     68 */	sxar2
/*     68 */	ldd	[%xg14],%f124
/*     68 */	or	%xg26,%l44(.LS0.cnt.39),%xg26


/*     68 */	sxar2
/*     68 */	or	%xg27,%l44(.LS0.cnt.38),%xg27
/*     68 */	ldd	[%xg15],%f126


/*     68 */	sxar2
/*     68 */	or	%xg28,%l44(.LS0.cnt.37),%xg28
/*     68 */	or	%xg29,%l44(.LS0.cnt.35),%xg29


/*     68 */	sxar2
/*     68 */	ldd	[%xg16],%f102
/*     68 */	or	%xg30,%l44(.LS0.cnt.49),%xg30


/*     68 */	sxar2
/*     68 */	or	%xg31,%l44(.LS0.cnt.50),%xg31
/*     68 */	ldd	[%xg17],%f130

/*     68 */	or	%g2,%l44(.LS0.cnt.48),%g2

/*     68 */	or	%g3,%l44(.LR0.cnt.3),%g3

/*     68 */	sxar1
/*     68 */	ldd	[%xg18],%f138

/*     68 */	sethi	%hi(-1047553),%g4


/*     68 */	sxar2
/*     68 */	ldd	[%xg19],%f146
/*     68 */	sethi	%h44(__fj_dlogexp1_const_),%xg0

/*     68 */	sxar1
/*     68 */	ldd	[%xg20],%f150

/*     68 */	or	%g4,1023,%g4


/*     68 */	sxar2
/*     68 */	sethi	%h44(__fj_dlogexp2_const_),%xg2
/*     68 */	ldd	[%xg21],%f108

/*     68 */	sxar1
/*     68 */	ldd	[%xg22],%f154

/*     68 */	sethi	%hi(-1),%g5


/*     68 */	sxar2
/*     68 */	or	%xg0,%m44(__fj_dlogexp1_const_),%xg0
/*     68 */	ldd	[%xg23],%f144

/*     68 */	sxar1
/*     68 */	ldd	[%xg24],%f148

/*     68 */	sllx	%g4,32,%g4


/*     68 */	sxar2
/*     68 */	or	%xg2,%m44(__fj_dlogexp2_const_),%xg2
/*     68 */	ldd	[%xg25],%f160

/*     68 */	sxar1
/*     68 */	ldd	[%xg26],%f156

/*     68 */	or	%g5,1023,%g5

/*     68 */	sethi	%hi(1024),%o0


/*     68 */	sxar2
/*     68 */	ldd	[%xg27],%f166
/*     68 */	ldd	[%xg28],%f168


/*     68 */	sxar2
/*     68 */	sllx	%xg0,12,%xg0
/*     68 */	sllx	%xg2,12,%xg2


/*     68 */	sxar2
/*     68 */	ldd	[%xg29],%f134
/*     68 */	ldd	[%xg30],%f70


/*     68 */	sxar2
/*     68 */	mov	1,%xg1
/*     68 */	mov	16,%xg3

/*     68 */	sxar1
/*     68 */	ldd	[%xg31],%f174

/*     68 */	ldd	[%g2],%f62

/*     68 */	sxar1
/*     68 */	mov	%o5,%xg4

/*     68 */	add	%g5,%g4,%g5

/*     68 */	sxar1
/*     68 */	ldd	[%g3],%f188

/*     68 */	sllx	%o0,32,%o0


/*     68 */	sxar2
/*     68 */	or	%xg0,%l44(__fj_dlogexp1_const_),%xg0
/*     68 */	or	%xg2,%l44(__fj_dlogexp2_const_),%xg2

.L5560:

/*     27 */	ldd	[%o1+32],%f46



/*     27 */	sxar2
/*     27 */	ldd	[%o1+80],%f302
/*     27 */	ldd	[%o1],%f248



/*     27 */	sxar2
/*     27 */	ldd	[%o1+48],%f504
/*     27 */	ldd	[%o1+8],%f254

/*     27 */	sxar1
/*     27 */	ldd	[%o1+56],%f510


/*     27 */	ldd	[%o1+16],%f54



/*     26 */	sxar2
/*     26 */	ldd	[%o1+64],%f310
/*     26 */	std,s	%f180,[%sp+2479]


/*    121 */	sxar2
/*    121 */	fmuld,s	%f46,%f44,%f46
/*    121 */	frcpad,s	%f46,%f48


/*     38 */	sxar2
/*     38 */	fnmsubd,s	%f46,%f48,%f40,%f46
/*     38 */	fmuld,s	%f46,%f46,%f50


/*     32 */	sxar2
/*     32 */	faddd,s	%f46,%f50,%f52
/*     32 */	fmaddd,s	%f50,%f50,%f46,%f50


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f52,%f48,%f48,%f52
/*     32 */	fmaddd,s	%f50,%f52,%f48,%f50

.L5561:

/*     78 */	cmp	%o4,%g0

/*     78 */	ble	.L5565
	nop


.L5562:


/*     78 */	sxar2
/*     78 */	mov	%g0,%xg12
/*     78 */	fzero,s	%f78

/*     78 */	sxar1
/*     78 */	mov	%o4,%xg15

.L5563:


/*    268 */	sxar2
/*    268 */	add	%o3,%xg12,%xg13
/*    268 */	add	%xg12,40,%xg12


/*    219 */	sxar2
/*    219 */	ldd,s	[%xg13],%f252
/*    219 */	ldd,s	[%xg13+16],%f46


/*    267 */	sxar2
/*    267 */	subcc	%xg15,1,%xg15
/*    267 */	fmsubd,sc	%f508,%f40,%f248,%f250


/*    267 */	sxar2
/*    267 */	fmsubd,sc	%f46,%f40,%f254,%f42
/*    267 */	fmsubd,sc	%f302,%f40,%f54,%f46


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f250,%f250,%f176,%f250
/*     32 */	fmaddd,s	%f42,%f42,%f250,%f42


/*     83 */	sxar2
/*     83 */	fmaddd,s	%f46,%f46,%f42,%f46
/*     83 */	frsqrtad,s	%f46,%f58


/*     38 */	sxar2
/*     38 */	fmuld,s	%f46,%f178,%f60
/*     38 */	fmuld,s	%f50,%f46,%f46


/*     32 */	sxar2
/*     32 */	fmuld,s	%f58,%f58,%f66
/*     32 */	fnmsubd,s	%f60,%f66,%f178,%f66


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f58,%f66,%f58,%f58
/*     32 */	fmuld,s	%f58,%f58,%f68


/*     32 */	sxar2
/*     32 */	fnmsubd,s	%f60,%f68,%f178,%f68
/*     32 */	fmaddd,s	%f58,%f68,%f58,%f58


/*     38 */	sxar2
/*     38 */	fmuld,s	%f58,%f58,%f72
/*     38 */	fnmsubd,s	%f60,%f72,%f178,%f60


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f58,%f60,%f58,%f58
/*     38 */	fmuld,s	%f46,%f58,%f46


/*     57 */	sxar2
/*     57 */	fsubd,s	%f40,%f46,%f74
/*     57 */	fmaddd,s	%f182,%f46,%f184,%f76


/*     57 */	sxar2
/*     57 */	fmaxd,s	%f180,%f74,%f74
/*     57 */	fmaddd,s	%f46,%f76,%f186,%f76


/*     57 */	sxar2
/*     57 */	fmuld,s	%f74,%f74,%f74
/*     57 */	fmaddd,s	%f46,%f76,%f40,%f46


/*     57 */	sxar2
/*     57 */	fmuld,s	%f74,%f74,%f74
/*     57 */	fmuld,s	%f74,%f74,%f74


/*    226 */	sxar2
/*    226 */	fmuld,s	%f74,%f46,%f74
/*    226 */	fmaddd,sc	%f252,%f74,%f78,%f78

/*    268 */	bne,pt	%icc, .L5563
	nop


.L5564:

/*    268 */	sxar1
/*    268 */	std,s	%f78,[%sp+2479]

.L5565:


/*     38 */	sxar2
/*     38 */	fmuld,s	%f50,%f50,%f54
/*     38 */	fmuld,s	%f50,%f56,%f50


/*    102 */	sxar2
/*    102 */	ldd,s	[%sp+2479],%f72
/*    102 */	cmp	%xg1,%o2

/*    101 */	sxar1
/*    101 */	ldd	[%o1+24],%f76

/*    101 */	fcmpeqd	%f38,%f34,%f58

/*    101 */	fcmpgeed	%f34,%f62,%f60


/*    101 */	sxar2
/*    101 */	fcmpgeed	%f32,%f64,%f66
/*    101 */	fcmpeqd	%f32,%f64,%f68


/*    101 */	sxar2
/*    101 */	fmuld,s	%f54,%f50,%f54
/*    101 */	fselmovd	%f64,%f70,%f66,%f66



/*    101 */	sxar2
/*    101 */	fmuld,s	%f54,%f72,%f54
/*    101 */	frcpad	%f54,%f74


/*    101 */	sxar2
/*    101 */	std	%f54,[%xg4]
/*    101 */	fnmsubd	%f74,%f54,%f296,%f54


/*    101 */	sxar2
/*    101 */	fmuld	%f76,%f74,%f76
/*    101 */	fmaddd	%f54,%f54,%f54,%f78


/*    101 */	sxar2
/*    101 */	fmuld	%f54,%f54,%f80
/*    101 */	fmaddd	%f78,%f76,%f76,%f78


/*    101 */	sxar2
/*    101 */	fmaddd	%f80,%f80,%f54,%f80
/*    101 */	fmaddd	%f80,%f78,%f76,%f80


/*    101 */	sxar2
/*    101 */	fabsd	%f80,%f82
/*    101 */	fcmplted	%f80,%f64,%f84


/*    101 */	sxar2
/*    101 */	fcmpeqd	%f80,%f64,%f80
/*    101 */	fpmaddx	%f82,%f90,%f86,%f88


/*    101 */	sxar2
/*    101 */	std	%f82,[%sp+2247]
/*    101 */	fpmaddxhi	%f88,%f94,%f64,%f92


/*    101 */	sxar2
/*    101 */	fand	%f88,%f96,%f88
/*    101 */	ldx	[%sp+2247],%xg5


/*    101 */	sxar2
/*    101 */	fxtod	%f92,%f92
/*    101 */	for	%f88,%f296,%f88


/*    101 */	sxar2
/*    101 */	andn	%xg5,%g5,%xg5
/*    101 */	add	%xg5,%o0,%xg5


/*    101 */	sxar2
/*    101 */	srax	%xg5,43,%xg5
/*    101 */	sllx	%xg5,5,%xg5


/*    101 */	sxar2
/*    101 */	add	%xg5,16,%xg6
/*    101 */	ldd	[%xg0+%xg5],%f114


/*    101 */	sxar2
/*    101 */	add	%xg5,8,%xg5
/*    101 */	ldd	[%xg0+%xg6],%f128


/*    101 */	sxar2
/*    101 */	ldd	[%xg0+%xg5],%f106
/*    101 */	fsubd	%f92,%f98,%f92


/*    101 */	sxar2
/*    101 */	fsubd	%f88,%f100,%f88
/*    101 */	fmuld	%f106,%f102,%f104


/*    101 */	sxar2
/*    101 */	fand	%f106,%f108,%f110
/*    101 */	fmuld	%f114,%f88,%f112


/*    101 */	sxar2
/*    101 */	fsubd	%f106,%f110,%f106
/*    101 */	fmuld	%f112,%f112,%f116


/*    101 */	sxar2
/*    101 */	fmaddd	%f122,%f112,%f118,%f120
/*    101 */	fmaddd	%f112,%f120,%f124,%f120


/*    101 */	sxar2
/*    101 */	fmaddd	%f112,%f120,%f126,%f112
/*    101 */	fmaddd	%f116,%f112,%f128,%f116


/*    101 */	sxar2
/*    101 */	fmaddd	%f114,%f88,%f116,%f114
/*    101 */	fmaddd	%f114,%f102,%f104,%f104


/*    101 */	sxar2
/*    101 */	faddd	%f106,%f114,%f106
/*    101 */	faddd	%f104,%f92,%f104


/*    101 */	sxar2
/*    101 */	fmuld	%f104,%f32,%f104
/*    101 */	fcmpleed	%f104,%f130,%f132


/*    101 */	sxar2
/*    101 */	fcmpgeed	%f104,%f134,%f136
/*    101 */	fselmovd	%f64,%f104,%f132,%f104


/*    101 */	sxar2
/*    101 */	fmuld	%f104,%f138,%f104
/*    101 */	fdtox	%f104,%f104


/*    101 */	sxar2
/*    101 */	fxtod	%f104,%f140
/*    101 */	fpmaddxhi	%f104,%f144,%f64,%f142


/*    101 */	sxar2
/*    101 */	std	%f104,[%sp+2239]
/*    101 */	fmuld	%f140,%f146,%f140


/*    101 */	sxar2
/*    101 */	fpmaddx	%f142,%f148,%f296,%f142
/*    101 */	ldx	[%sp+2239],%xg7


/*    101 */	sxar2
/*    101 */	fmaddd	%f92,%f32,%f140,%f92
/*    101 */	and	%xg7,255,%xg7


/*    101 */	sxar2
/*    101 */	sllx	%xg7,4,%xg7
/*    101 */	ldd	[%xg2+%xg7],%f164


/*    101 */	sxar2
/*    101 */	add	%xg7,8,%xg7
/*    101 */	ldd	[%xg2+%xg7],%f170


/*    101 */	sxar2
/*    101 */	fmuld	%f92,%f150,%f152
/*    101 */	fmaddd	%f32,%f110,%f152,%f110


/*    101 */	sxar2
/*    101 */	fmaddd	%f32,%f106,%f110,%f106
/*    101 */	fmaddd	%f154,%f92,%f106,%f92


/*    101 */	sxar2
/*    101 */	fmaddd	%f160,%f92,%f156,%f158
/*    101 */	fmuld	%f164,%f92,%f162


/*    101 */	sxar2
/*    101 */	fmaddd	%f158,%f92,%f166,%f158
/*    101 */	fmaddd	%f158,%f92,%f168,%f158


/*    101 */	sxar2
/*    101 */	fmaddd	%f158,%f162,%f170,%f158
/*    101 */	faddd	%f164,%f158,%f164


/*    101 */	sxar2
/*    101 */	fmuld	%f142,%f164,%f142
/*    101 */	fselmovd	%f64,%f142,%f132,%f142


/*    101 */	sxar2
/*    101 */	fselmovd	%f70,%f142,%f136,%f142
/*    101 */	fmuld	%f142,%f36,%f172


/*    101 */	sxar2
/*    101 */	fselmovd	%f172,%f174,%f58,%f172
/*    101 */	fselmovd	%f142,%f172,%f60,%f172


/*    101 */	sxar2
/*    101 */	fselmovd	%f172,%f142,%f84,%f172
/*    101 */	fselmovd	%f66,%f172,%f80,%f66


/*    101 */	sxar2
/*    101 */	fselmovd	%f296,%f66,%f68,%f66
/*    101 */	fmuld	%f66,%f188,%f66

/*    101 */	sxar1
/*    101 */	std	%f66,[%xg4+8]

/*    102 */	bge	.L5567
	nop


.L5566:


/*    104 */	sxar2
/*    104 */	add	%o5,%xg3,%xg8
/*    104 */	ldd,d	[%o1+72],%f198


/*    104 */	sxar2
/*    104 */	fcmpgeed	%f32,%f64,%f190
/*    104 */	fcmpeqd	%f32,%f64,%f192


/*    104 */	sxar2
/*    104 */	std	%f310,[%xg8]
/*    104 */	fselmovd	%f64,%f70,%f190,%f190

/*    104 */	prefetch	[%o1+1352],21

/*    104 */	sxar1
/*    104 */	ldd,d	[%xg4+16],%f196

/*    104 */	prefetch	[%o1+328],0


/*    104 */	sxar2
/*    104 */	prefetch	[%xg4+1296],23
/*    104 */	prefetch	[%xg4+272],2


/*    104 */	sxar2
/*    104 */	frcpad	%f196,%f194
/*    104 */	fnmsubd	%f194,%f196,%f296,%f196


/*    104 */	sxar2
/*    104 */	fmuld	%f198,%f194,%f198
/*    104 */	fmaddd	%f196,%f196,%f196,%f200


/*    104 */	sxar2
/*    104 */	fmuld	%f196,%f196,%f202
/*    104 */	fmaddd	%f200,%f198,%f198,%f200


/*    104 */	sxar2
/*    104 */	fmaddd	%f202,%f202,%f196,%f202
/*    104 */	fmaddd	%f202,%f200,%f198,%f202


/*    104 */	sxar2
/*    104 */	fabsd	%f202,%f204
/*    104 */	fcmplted	%f202,%f64,%f206


/*    104 */	sxar2
/*    104 */	fcmpeqd	%f202,%f64,%f202
/*    104 */	fpmaddx	%f204,%f90,%f86,%f208


/*    104 */	sxar2
/*    104 */	std	%f204,[%sp+2247]
/*    104 */	fpmaddxhi	%f208,%f94,%f64,%f210


/*    104 */	sxar2
/*    104 */	fand	%f208,%f96,%f208
/*    104 */	ldx	[%sp+2247],%xg9


/*    104 */	sxar2
/*    104 */	fxtod	%f210,%f210
/*    104 */	for	%f208,%f296,%f208


/*    104 */	sxar2
/*    104 */	andn	%xg9,%g5,%xg9
/*    104 */	add	%xg9,%o0,%xg9


/*    104 */	sxar2
/*    104 */	srax	%xg9,43,%xg9
/*    104 */	sllx	%xg9,5,%xg9


/*    104 */	sxar2
/*    104 */	add	%xg9,16,%xg10
/*    104 */	ldd	[%xg0+%xg9],%f220


/*    104 */	sxar2
/*    104 */	add	%xg9,8,%xg9
/*    104 */	ldd	[%xg0+%xg10],%f226


/*    104 */	sxar2
/*    104 */	ldd	[%xg0+%xg9],%f214
/*    104 */	fsubd	%f210,%f98,%f210


/*    104 */	sxar2
/*    104 */	fsubd	%f208,%f100,%f208
/*    104 */	fmuld	%f214,%f102,%f212


/*    104 */	sxar2
/*    104 */	fand	%f214,%f108,%f216
/*    104 */	fmuld	%f220,%f208,%f218


/*    104 */	sxar2
/*    104 */	fsubd	%f214,%f216,%f214
/*    104 */	fmuld	%f218,%f218,%f222


/*    104 */	sxar2
/*    104 */	fmaddd	%f122,%f218,%f118,%f224
/*    104 */	fmaddd	%f218,%f224,%f124,%f224


/*    104 */	sxar2
/*    104 */	fmaddd	%f218,%f224,%f126,%f218
/*    104 */	fmaddd	%f222,%f218,%f226,%f222


/*    104 */	sxar2
/*    104 */	fmaddd	%f220,%f208,%f222,%f220
/*    104 */	fmaddd	%f220,%f102,%f212,%f212


/*    104 */	sxar2
/*    104 */	faddd	%f214,%f220,%f214
/*    104 */	faddd	%f212,%f210,%f212


/*    104 */	sxar2
/*    104 */	fmuld	%f212,%f32,%f212
/*    104 */	fcmpleed	%f212,%f130,%f228


/*    104 */	sxar2
/*    104 */	fcmpgeed	%f212,%f134,%f230
/*    104 */	fselmovd	%f64,%f212,%f228,%f212


/*    104 */	sxar2
/*    104 */	fmuld	%f212,%f138,%f212
/*    104 */	fdtox	%f212,%f212


/*    104 */	sxar2
/*    104 */	fxtod	%f212,%f232
/*    104 */	fpmaddxhi	%f212,%f144,%f64,%f234


/*    104 */	sxar2
/*    104 */	std	%f212,[%sp+2239]
/*    104 */	fmuld	%f232,%f146,%f232


/*    104 */	sxar2
/*    104 */	fpmaddx	%f234,%f148,%f296,%f234
/*    104 */	ldx	[%sp+2239],%xg11


/*    104 */	sxar2
/*    104 */	fmaddd	%f210,%f32,%f232,%f210
/*    104 */	and	%xg11,255,%xg11


/*    104 */	sxar2
/*    104 */	sllx	%xg11,4,%xg11
/*    104 */	ldd	[%xg2+%xg11],%f242


/*    104 */	sxar2
/*    104 */	add	%xg11,8,%xg11
/*    104 */	ldd	[%xg2+%xg11],%f244


/*    104 */	sxar2
/*    104 */	fmuld	%f210,%f150,%f236
/*    104 */	fmaddd	%f32,%f216,%f236,%f216


/*    104 */	sxar2
/*    104 */	fmaddd	%f32,%f214,%f216,%f214
/*    104 */	fmaddd	%f154,%f210,%f214,%f210


/*    104 */	sxar2
/*    104 */	fmaddd	%f160,%f210,%f156,%f238
/*    104 */	fmuld	%f242,%f210,%f240


/*    104 */	sxar2
/*    104 */	fmaddd	%f238,%f210,%f166,%f238
/*    104 */	fmaddd	%f238,%f210,%f168,%f238


/*    104 */	sxar2
/*    104 */	fmaddd	%f238,%f240,%f244,%f238
/*    104 */	faddd	%f242,%f238,%f242


/*    104 */	sxar2
/*    104 */	fmuld	%f234,%f242,%f234
/*    104 */	fselmovd	%f64,%f234,%f228,%f234


/*    104 */	sxar2
/*    104 */	fselmovd	%f70,%f234,%f230,%f234
/*    104 */	fmuld	%f234,%f36,%f246


/*    104 */	sxar2
/*    104 */	fselmovd	%f246,%f174,%f58,%f246
/*    104 */	fselmovd	%f234,%f246,%f60,%f246


/*    104 */	sxar2
/*    104 */	fselmovd	%f246,%f234,%f206,%f246
/*    104 */	fselmovd	%f190,%f246,%f202,%f190


/*    104 */	sxar2
/*    104 */	fselmovd	%f296,%f190,%f192,%f190
/*    104 */	fmuld	%f190,%f188,%f190

/*    104 */	sxar1
/*    104 */	std,d	%f190,[%xg4+24]

.L5567:


/*    106 */	sxar2
/*    106 */	add	%xg1,2,%xg1
/*    106 */	add	%xg3,32,%xg3

/*    106 */	sxar1
/*    106 */	add	%xg4,32,%xg4


/*    106 */	subcc	%g1,1,%g1

/*    106 */	bne,pt	%icc, .L5560
/*    106 */	add	%o1,96,%o1


.L5568:


.L5570:

/*    107 */	sub	%sp,-2048,%sp
	retl
	nop


.LLFE1:
	.size	_ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE,.-_ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE
	.type	_ZN11CalcDensityclEPKN3EPI4DensEiPKN3EPJ4DensEiPN6RESULT4DensE,#function
	.ident	"$Compiler: Fujitsu C/C++ Compiler Version 1.2.0 P-id: L30000-09 (Sep 17 2014 11:41:29) compile.C _ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE $"
	.section	".text"
	.global	_ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE
	.align	64
_ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE:
.LLFB2:
.L5571:

/*    112 */	add	%sp,-2784,%sp


.L5572:

/*    119 */	sub	%o2,1,%g1

/*    119 */	cmp	%o2,%g0

/*    119 */	srl	%g1,31,%g2

/*    119 */	add	%g1,%g2,%g1

/*    119 */	sra	%g1,1,%g1


/*    119 */	ble	.L5584
/*    119 */	add	%g1,1,%g1


.L5573:

/*    119 */	sethi	%h44(.LR0.cnt.13),%g3

/*    119 */	sethi	%h44(.LR0.cnt.7),%g4

/*    119 */	sxar1
/*    119 */	fzero,s	%f106

/*    119 */	sethi	%h44(.LR0.cnt.5),%g5

/*    119 */	sethi	%h44(.LR0.cnt.8),%o0


/*    119 */	sxar2
/*    119 */	sethi	%h44(.LR0.cnt.10),%xg0
/*    119 */	sethi	%h44(.LR0.cnt.16),%xg1


/*    119 */	sxar2
/*    119 */	sethi	%h44(.LR0.cnt.17),%xg2
/*    119 */	sethi	%h44(.LR0.cnt.18),%xg3


/*    119 */	sxar2
/*    119 */	sethi	%h44(.LR0.cnt.19),%xg4
/*    119 */	sethi	%h44(.LR0.cnt.11),%xg5

/*    119 */	or	%g3,%m44(.LR0.cnt.13),%g3

/*    119 */	sxar1
/*    119 */	sethi	%h44(_ZN6kernel14normalzeFactorE),%xg6

/*    119 */	or	%g4,%m44(.LR0.cnt.7),%g4

/*    119 */	or	%g5,%m44(.LR0.cnt.5),%g5

/*    119 */	or	%o0,%m44(.LR0.cnt.8),%o0


/*    119 */	sxar2
/*    119 */	or	%xg0,%m44(.LR0.cnt.10),%xg0
/*    119 */	or	%xg1,%m44(.LR0.cnt.16),%xg1


/*    119 */	sxar2
/*    119 */	or	%xg2,%m44(.LR0.cnt.17),%xg2
/*    119 */	or	%xg3,%m44(.LR0.cnt.18),%xg3


/*    119 */	sxar2
/*    119 */	or	%xg4,%m44(.LR0.cnt.19),%xg4
/*    119 */	or	%xg5,%m44(.LR0.cnt.11),%xg5

/*    119 */	sxar1
/*    119 */	or	%xg6,%m44(_ZN6kernel14normalzeFactorE),%xg6

/*    119 */	sllx	%g3,12,%g3

/*    119 */	sllx	%g4,12,%g4

/*    119 */	sllx	%g5,12,%g5

/*    119 */	sllx	%o0,12,%o0


/*    119 */	sxar2
/*    119 */	sllx	%xg0,12,%xg0
/*    119 */	sllx	%xg1,12,%xg1


/*    119 */	sxar2
/*    119 */	sllx	%xg2,12,%xg2
/*    119 */	sllx	%xg3,12,%xg3


/*    119 */	sxar2
/*    119 */	sllx	%xg4,12,%xg4
/*    119 */	sllx	%xg5,12,%xg5

/*    119 */	or	%g3,%l44(.LR0.cnt.13),%g3

/*    119 */	sxar1
/*    119 */	sllx	%xg6,12,%xg6

/*    119 */	or	%g4,%l44(.LR0.cnt.7),%g4

/*    119 */	ldd	[%g3],%f32

/*    119 */	or	%g5,%l44(.LR0.cnt.5),%g5

/*    119 */	or	%o0,%l44(.LR0.cnt.8),%o0


/*    119 */	sxar2
/*    119 */	ldd	[%g3],%f288
/*    119 */	or	%xg0,%l44(.LR0.cnt.10),%xg0


/*    119 */	sxar1
/*    119 */	or	%xg1,%l44(.LR0.cnt.16),%xg1

/*    119 */	ldd	[%g4],%f38


/*    119 */	sxar2
/*    119 */	or	%xg2,%l44(.LR0.cnt.17),%xg2
/*    119 */	or	%xg3,%l44(.LR0.cnt.18),%xg3


/*    119 */	sxar2
/*    119 */	ldd	[%g4],%f294
/*    119 */	or	%xg4,%l44(.LR0.cnt.19),%xg4



/*    119 */	sxar2
/*    119 */	or	%xg5,%l44(.LR0.cnt.11),%xg5
/*    119 */	ldd	[%g5],%f88


/*    119 */	sxar2
/*    119 */	or	%xg6,%l44(_ZN6kernel14normalzeFactorE),%xg6
/*    119 */	mov	1,%xg23


/*    119 */	sxar2
/*    119 */	ldd	[%g5],%f344
/*    119 */	ldd	[%o0],%f94



/*    119 */	sxar2
/*    119 */	add	%o5,32,%xg7
/*    119 */	ldd	[%o0],%f350



/*    119 */	sxar2
/*    119 */	ldd	[%xg0],%f128
/*    119 */	ldd	[%xg0],%f384



/*    119 */	sxar2
/*    119 */	ldd	[%xg1],%f112
/*    119 */	ldd	[%xg1],%f368



/*    119 */	sxar2
/*    119 */	ldd	[%xg2],%f116
/*    119 */	ldd	[%xg2],%f372



/*    119 */	sxar2
/*    119 */	ldd	[%xg3],%f130
/*    119 */	ldd	[%xg3],%f386



/*    119 */	sxar2
/*    119 */	ldd	[%xg4],%f118
/*    119 */	ldd	[%xg4],%f374



/*    119 */	sxar2
/*    119 */	ldd	[%xg5],%f122
/*    119 */	ldd	[%xg5],%f378



/*    119 */	sxar2
/*    119 */	ldd	[%xg6],%f52
/*    119 */	ldd	[%xg6],%f308


.L5574:

/*     27 */	ldd	[%o1+48],%f34



/*     26 */	sxar2
/*     26 */	ldd	[%o1+112],%f290
/*     26 */	std,s	%f106,[%sp+2527]


/*     27 */	sxar2
/*     27 */	ldd	[%o1],%f64
/*     27 */	ldd	[%o1+64],%f320



/*     27 */	sxar2
/*     27 */	ldd	[%o1+8],%f68
/*     27 */	ldd	[%o1+72],%f324



/*     27 */	sxar2
/*     27 */	ldd	[%o1+16],%f74
/*     27 */	ldd	[%o1+80],%f330



/*     38 */	sxar2
/*     38 */	ldd	[%o1+24],%f76
/*     38 */	fmuld,s	%f34,%f32,%f34


/*     27 */	sxar2
/*     27 */	ldd	[%o1+88],%f332
/*     27 */	ldd	[%o1+32],%f82



/*     27 */	sxar2
/*     27 */	std,s	%f106,[%sp+2543]
/*     27 */	ldd	[%o1+96],%f338



/*     27 */	sxar2
/*     27 */	ldd	[%o1+40],%f84
/*     27 */	ldd	[%o1+104],%f340



/*     26 */	sxar2
/*     26 */	frcpad,s	%f34,%f36
/*     26 */	std,s	%f106,[%sp+2559]


/*     26 */	sxar2
/*     26 */	fnmsubd,s	%f34,%f36,%f38,%f34
/*     26 */	std,s	%f106,[%sp+2575]


/*     32 */	sxar2
/*     32 */	fmuld,s	%f34,%f34,%f40
/*     32 */	faddd,s	%f34,%f40,%f42


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f40,%f40,%f34,%f40
/*     32 */	fmaddd,s	%f42,%f36,%f36,%f42

/*     32 */	sxar1
/*     32 */	fmaddd,s	%f40,%f42,%f36,%f40

.L5575:

/*    135 */	cmp	%o4,%g0

/*    135 */	ble	.L5579
	nop


.L5576:


/*    135 */	sxar2
/*    135 */	mov	%g0,%xg8
/*    135 */	fzero,s	%f138


/*    135 */	sxar2
/*    135 */	mov	%o4,%xg14
/*    135 */	fmovd,s	%f138,%f136


/*    135 */	sxar2
/*    135 */	fmovd,s	%f138,%f134
/*    135 */	fmovd,s	%f138,%f132

.L5577:


/*    268 */	sxar2
/*    268 */	add	%o3,%xg8,%xg9
/*    268 */	add	%xg8,64,%xg8


/*    219 */	sxar2
/*    219 */	ldd,s	[%xg9],%f66
/*    219 */	ldd,s	[%xg9+16],%f72


/*    219 */	sxar2
/*    219 */	subcc	%xg14,1,%xg14
/*    219 */	ldd,s	[%xg9+32],%f80


/*    267 */	sxar2
/*    267 */	ldd,s	[%xg9+48],%f86
/*    267 */	fmsubd,sc	%f322,%f38,%f64,%f62


/*    267 */	sxar2
/*    267 */	fmsubd,sc	%f72,%f38,%f68,%f70
/*    267 */	fmsubd,sc	%f328,%f38,%f74,%f72


/*    267 */	sxar2
/*    267 */	fmsubd,sc	%f80,%f38,%f76,%f78
/*    267 */	fmsubd,sc	%f336,%f38,%f82,%f80


/*     32 */	sxar2
/*     32 */	fmsubd,sc	%f86,%f38,%f84,%f86
/*     32 */	fmaddd,s	%f62,%f62,%f88,%f90


/*     53 */	sxar2
/*     53 */	fmuld,s	%f80,%f70,%f124
/*     53 */	fmuld,s	%f80,%f72,%f126


/*     53 */	sxar2
/*     53 */	fmaddd,s	%f70,%f70,%f90,%f90
/*     53 */	fmaddd,s	%f78,%f62,%f124,%f124


/*     53 */	sxar2
/*     53 */	fmsubd,s	%f86,%f70,%f126,%f126
/*     53 */	fmuld,s	%f78,%f70,%f70


/*     53 */	sxar2
/*     53 */	fmaddd,s	%f72,%f72,%f90,%f90
/*     53 */	fmaddd,s	%f86,%f72,%f124,%f124


/*     53 */	sxar2
/*     53 */	fmuld,s	%f86,%f62,%f86
/*     53 */	fmsubd,s	%f80,%f62,%f70,%f80


/*     38 */	sxar2
/*     38 */	frsqrtad,s	%f90,%f92
/*     38 */	fmuld,s	%f90,%f94,%f96


/*     53 */	sxar2
/*     53 */	fmuld,s	%f40,%f90,%f90
/*     53 */	fmsubd,s	%f78,%f72,%f86,%f78


/*     32 */	sxar2
/*     32 */	fmuld,s	%f92,%f92,%f98
/*     32 */	fnmsubd,s	%f96,%f98,%f94,%f98


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f92,%f98,%f92,%f92
/*     32 */	fmuld,s	%f92,%f92,%f100


/*     32 */	sxar2
/*     32 */	fnmsubd,s	%f96,%f100,%f94,%f100
/*     32 */	fmaddd,s	%f92,%f100,%f92,%f92


/*     32 */	sxar2
/*     32 */	fmuld,s	%f92,%f92,%f102
/*     32 */	fnmsubd,s	%f96,%f102,%f94,%f96


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f92,%f96,%f92,%f92
/*     38 */	fmuld,s	%f92,%f90,%f90


/*     38 */	sxar2
/*     38 */	fsubd,s	%f38,%f90,%f104
/*     38 */	fmaddd,s	%f116,%f90,%f112,%f114


/*    192 */	sxar2
/*    192 */	fmaddd,s	%f122,%f90,%f118,%f120
/*    192 */	fmaxd,s	%f106,%f104,%f104


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f90,%f114,%f128,%f114
/*     38 */	fmaddd,s	%f90,%f120,%f130,%f120


/*     38 */	sxar2
/*     38 */	fmuld,s	%f104,%f104,%f108
/*     38 */	fmuld,s	%f104,%f108,%f110


/*     38 */	sxar2
/*     38 */	fmuld,s	%f108,%f108,%f108
/*     38 */	fmsubd,s	%f104,%f114,%f128,%f104


/*     38 */	sxar2
/*     38 */	fmuld,s	%f92,%f110,%f92
/*     38 */	fnmsubd,s	%f90,%f120,%f104,%f90


/*     38 */	sxar2
/*     38 */	fmuld,s	%f108,%f92,%f108
/*     38 */	fmuld,s	%f108,%f90,%f108


/*     53 */	sxar2
/*     53 */	fmaddd,sc	%f66,%f108,%f106,%f66
/*     53 */	fnmsubd,s	%f66,%f124,%f132,%f132


/*     53 */	sxar2
/*     53 */	fmaddd,s	%f66,%f126,%f134,%f134
/*     53 */	fmaddd,s	%f66,%f78,%f136,%f136

/*     53 */	sxar1
/*     53 */	fmaddd,s	%f66,%f80,%f138,%f138

/*    268 */	bne,pt	%icc, .L5577
	nop


.L5578:


/*    268 */	sxar2
/*    268 */	std,s	%f132,[%sp+2527]
/*    268 */	std,s	%f134,[%sp+2543]


/*    268 */	sxar2
/*    268 */	std,s	%f136,[%sp+2559]
/*    268 */	std,s	%f138,[%sp+2575]

.L5579:

/*     27 */	ldd	[%o1+56],%f46


/*     38 */	sxar2
/*     38 */	ldd	[%o1+120],%f302
/*     38 */	fmuld,s	%f40,%f40,%f40



/*     57 */	sxar2
/*     57 */	cmp	%xg23,%o2
/*     57 */	ldd,s	[%sp+2527],%f54


/*     57 */	sxar2
/*     57 */	ldd,s	[%sp+2543],%f56
/*     57 */	ldd,s	[%sp+2559],%f58


/*    121 */	sxar2
/*    121 */	ldd,s	[%sp+2575],%f60
/*    121 */	frcpad,s	%f46,%f44


/*     35 */	sxar2
/*     35 */	fmuld,s	%f40,%f40,%f40
/*     35 */	fnmsubd,s	%f46,%f44,%f38,%f46


/*     32 */	sxar2
/*     32 */	fmuld,s	%f46,%f46,%f48
/*     32 */	faddd,s	%f46,%f48,%f50


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f48,%f48,%f46,%f48
/*     32 */	fmaddd,s	%f50,%f44,%f44,%f50


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f48,%f50,%f44,%f48
/*     38 */	fmuld,s	%f48,%f52,%f48


/*     57 */	sxar2
/*     57 */	fmuld,s	%f48,%f40,%f48
/*     57 */	fmuld,s	%f48,%f54,%f54


/*     57 */	sxar2
/*     57 */	fmuld,s	%f48,%f56,%f56
/*     57 */	fmuld,s	%f48,%f58,%f58

/*     57 */	sxar1
/*     57 */	fmuld,s	%f48,%f60,%f48



/*    149 */	std	%f54,[%o5]

/*    149 */	std	%f56,[%o5+8]



/*    149 */	std	%f58,[%o5+16]


/*    179 */	bge	.L5581
/*    179 */	std	%f48,[%o5+24]


.L5580:


/*    152 */	sxar2
/*    152 */	std	%f310,[%xg7]
/*    152 */	std	%f312,[%xg7+8]


/*    152 */	sxar2
/*    152 */	std	%f314,[%xg7+16]
/*    152 */	std	%f304,[%xg7+24]

.L5581:


/*    185 */	sxar2
/*    185 */	add	%xg23,2,%xg23
/*    185 */	add	%xg7,64,%xg7

/*    185 */	add	%o5,64,%o5


/*    185 */	subcc	%g1,1,%g1

/*    185 */	bne,pt	%icc, .L5574
/*    185 */	add	%o1,128,%o1


.L5582:


.L5584:

/*    186 */	sub	%sp,-2784,%sp
	retl
	nop


.LLFE2:
	.size	_ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE,.-_ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE
	.type	_ZN14CalcDerivativeclEPKN3EPI4DrvtEiPKN3EPJ4DrvtEiPN6RESULT4DrvtE,#function
	.ident	"$Compiler: Fujitsu C/C++ Compiler Version 1.2.0 P-id: L30000-09 (Sep 17 2014 11:41:29) compile.C _ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE $"
	.section	".text"
	.global	_ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE
	.align	64
_ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE:
.LLFB3:
.L5585:

/*    191 */	add	%sp,-3408,%sp


.L5586:

/*    211 */	cmp	%o4,%g0

/*    211 */	ble	.L5591
	nop


.L5587:

/*    211 */	sethi	%h44(.LR0.cnt.7),%g1

/*    211 */	sethi	%h44(.LR0.cnt.1),%g2

/*    211 */	or	%g1,%m44(.LR0.cnt.7),%g1

/*    211 */	or	%g2,%m44(.LR0.cnt.1),%g2

/*    211 */	sllx	%g1,12,%g1

/*    211 */	sllx	%g2,12,%g2

/*    211 */	or	%g1,%l44(.LR0.cnt.7),%g1

/*    211 */	or	%g2,%l44(.LR0.cnt.1),%g2

/*    211 */	sxar1
/*    211 */	mov	%g0,%xg14

/*    211 */	ldd	[%g1],%f44

/*    211 */	ldd	[%g2],%f46

/*    211 */	sxar1
/*    211 */	mov	%o4,%xg17

.L5588:


/*    218 */	sxar2
/*    218 */	add	%xg14,%o3,%xg13
/*    218 */	add	%xg14,104,%xg14


/*    216 */	sxar2
/*    216 */	ldd	[%xg13+48],%f36
/*    216 */	ldd	[%xg13+64],%f38


/*    215 */	sxar2
/*    215 */	subcc	%xg17,1,%xg17
/*    215 */	ldd	[%xg13+72],%f50

/*    215 */	fmuld	%f36,%f36,%f36

/*    216 */	frcpad	%f38,%f40

/*    215 */	frcpad	%f36,%f42

/*    216 */	fnmsubd	%f40,%f38,%f44,%f38

/*    216 */	fmuld	%f40,%f46,%f40

/*    215 */	fnmsubd	%f42,%f36,%f44,%f36

/*    216 */	fmaddd	%f38,%f38,%f38,%f48

/*    215 */	fmuld	%f50,%f42,%f50

/*    216 */	fmuld	%f38,%f38,%f52

/*    215 */	fmaddd	%f36,%f36,%f36,%f54

/*    215 */	fmuld	%f36,%f36,%f56

/*    216 */	fmaddd	%f48,%f40,%f40,%f48

/*    216 */	fmaddd	%f52,%f52,%f38,%f52

/*    215 */	fmaddd	%f54,%f50,%f50,%f54

/*    215 */	fmaddd	%f56,%f56,%f36,%f56

/*    216 */	fmaddd	%f52,%f48,%f40,%f52

/*    215 */	fmaddd	%f56,%f54,%f50,%f56


/*    217 */	sxar2
/*    217 */	std	%f52,[%xg13+64]
/*    217 */	std	%f56,[%xg13+72]

/*    218 */	bne,pt	%icc, .L5588
	nop


.L5589:


.L5591:

/*    220 */	sub	%o2,1,%g3

/*    220 */	cmp	%o2,%g0

/*    220 */	srl	%g3,31,%g4

/*    220 */	add	%g3,%g4,%g3

/*    220 */	sra	%g3,1,%g3


/*    220 */	ble	.L5603
/*    220 */	add	%g3,1,%g3


.L5592:

/*    220 */	sxar1
/*    220 */	fzero,s	%f146

/*    220 */	sethi	%h44(.LR0.cnt.8),%g5

/*    220 */	sethi	%h44(_ZN6kernel14normalzeFactorE),%o0

/*    220 */	or	%g5,%m44(.LR0.cnt.8),%g5

/*    220 */	or	%o0,%m44(_ZN6kernel14normalzeFactorE),%o0

/*    220 */	sllx	%g5,12,%g5

/*    220 */	sllx	%o0,12,%o0

/*    220 */	or	%g5,%l44(.LR0.cnt.8),%g5

/*    220 */	or	%o0,%l44(_ZN6kernel14normalzeFactorE),%o0


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.13),%xg0
/*    220 */	sethi	%h44(.LR0.cnt.7),%xg1

/*    220 */	ldd	[%g5],%f34

/*    220 */	ldd	[%o0],%f32


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.5),%xg2
/*    220 */	sethi	%h44(.LR0.cnt.9),%xg3


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.10),%xg4
/*    220 */	sethi	%h44(.LR0.cnt.16),%xg5


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.17),%xg6
/*    220 */	sethi	%h44(.LR0.cnt.18),%xg7


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.19),%xg8
/*    220 */	sethi	%h44(.LR0.cnt.11),%xg9



/*    220 */	sxar2
/*    220 */	fmovd	%f34,%f290
/*    220 */	sethi	%h44(.LR0.cnt.21),%xg10


/*    220 */	sxar2
/*    220 */	sethi	%h44(.LR0.cnt.4),%xg11
/*    220 */	or	%xg0,%m44(.LR0.cnt.13),%xg0

/*    220 */	sxar1
/*    220 */	or	%xg1,%m44(.LR0.cnt.7),%xg1

/*    220 */	fmuld	%f32,%f34,%f32


/*    220 */	sxar2
/*    220 */	or	%xg2,%m44(.LR0.cnt.5),%xg2
/*    220 */	or	%xg3,%m44(.LR0.cnt.9),%xg3


/*    220 */	sxar2
/*    220 */	or	%xg4,%m44(.LR0.cnt.10),%xg4
/*    220 */	or	%xg5,%m44(.LR0.cnt.16),%xg5


/*    220 */	sxar2
/*    220 */	or	%xg6,%m44(.LR0.cnt.17),%xg6
/*    220 */	or	%xg7,%m44(.LR0.cnt.18),%xg7


/*    220 */	sxar2
/*    220 */	or	%xg8,%m44(.LR0.cnt.19),%xg8
/*    220 */	or	%xg9,%m44(.LR0.cnt.11),%xg9


/*    220 */	sxar2
/*    220 */	or	%xg10,%m44(.LR0.cnt.21),%xg10
/*    220 */	or	%xg11,%m44(.LR0.cnt.4),%xg11


/*    220 */	sxar2
/*    220 */	sllx	%xg0,12,%xg0
/*    220 */	sllx	%xg1,12,%xg1



/*    220 */	sxar2
/*    220 */	fmovd	%f32,%f288
/*    220 */	sllx	%xg2,12,%xg2


/*    220 */	sxar2
/*    220 */	sllx	%xg3,12,%xg3
/*    220 */	sllx	%xg4,12,%xg4


/*    220 */	sxar2
/*    220 */	sllx	%xg5,12,%xg5
/*    220 */	sllx	%xg6,12,%xg6


/*    220 */	sxar2
/*    220 */	sllx	%xg7,12,%xg7
/*    220 */	sllx	%xg8,12,%xg8


/*    220 */	sxar2
/*    220 */	sllx	%xg9,12,%xg9
/*    220 */	sllx	%xg10,12,%xg10


/*    220 */	sxar2
/*    220 */	sllx	%xg11,12,%xg11
/*    220 */	or	%xg0,%l44(.LR0.cnt.13),%xg0


/*    220 */	sxar2
/*    220 */	or	%xg1,%l44(.LR0.cnt.7),%xg1
/*    220 */	ldd	[%xg0],%f58


/*    220 */	sxar2
/*    220 */	or	%xg2,%l44(.LR0.cnt.5),%xg2
/*    220 */	or	%xg3,%l44(.LR0.cnt.9),%xg3


/*    220 */	sxar2
/*    220 */	ldd	[%xg0],%f314
/*    220 */	or	%xg4,%l44(.LR0.cnt.10),%xg4



/*    220 */	sxar2
/*    220 */	or	%xg5,%l44(.LR0.cnt.16),%xg5
/*    220 */	ldd	[%xg1],%f66


/*    220 */	sxar2
/*    220 */	or	%xg6,%l44(.LR0.cnt.17),%xg6
/*    220 */	or	%xg7,%l44(.LR0.cnt.18),%xg7


/*    220 */	sxar2
/*    220 */	ldd	[%xg1],%f322
/*    220 */	or	%xg8,%l44(.LR0.cnt.19),%xg8



/*    220 */	sxar2
/*    220 */	or	%xg9,%l44(.LR0.cnt.11),%xg9
/*    220 */	ldd	[%xg2],%f120


/*    220 */	sxar2
/*    220 */	or	%xg10,%l44(.LR0.cnt.21),%xg10
/*    220 */	or	%xg11,%l44(.LR0.cnt.4),%xg11


/*    220 */	sxar2
/*    220 */	ldd	[%xg2],%f376
/*    220 */	ldd	[%xg3],%f160

/*    220 */	mov	1,%g4



/*    220 */	sxar2
/*    220 */	add	%o5,40,%xg12
/*    220 */	ldd	[%xg3],%f416



/*    220 */	sxar2
/*    220 */	ldd	[%xg4],%f204
/*    220 */	ldd	[%xg4],%f460



/*    220 */	sxar2
/*    220 */	ldd	[%xg5],%f182
/*    220 */	ldd	[%xg5],%f438



/*    220 */	sxar2
/*    220 */	ldd	[%xg6],%f186
/*    220 */	ldd	[%xg6],%f442



/*    220 */	sxar2
/*    220 */	ldd	[%xg7],%f206
/*    220 */	ldd	[%xg7],%f462



/*    220 */	sxar2
/*    220 */	ldd	[%xg8],%f190
/*    220 */	ldd	[%xg8],%f446



/*    220 */	sxar2
/*    220 */	ldd	[%xg9],%f194
/*    220 */	ldd	[%xg9],%f450



/*    220 */	sxar2
/*    220 */	ldd	[%xg10],%f208
/*    220 */	ldd	[%xg10],%f464



/*    220 */	sxar2
/*    220 */	ldd	[%xg11],%f86
/*    220 */	ldd	[%xg11],%f342


.L5593:

/*     27 */	ldd	[%o1+48],%f62



/*     27 */	sxar2
/*     27 */	ldd	[%o1+144],%f318
/*     27 */	ldd	[%o1+56],%f72



/*     26 */	sxar2
/*     26 */	ldd	[%o1+152],%f328
/*     26 */	std,s	%f146,[%sp+2623]


/*     27 */	sxar2
/*     27 */	ldd	[%o1+64],%f84
/*     27 */	ldd	[%o1+160],%f340



/*     27 */	sxar2
/*     27 */	ldd	[%o1],%f106
/*     27 */	ldd	[%o1+96],%f362



/*     38 */	sxar2
/*     38 */	ldd	[%o1+8],%f110
/*     38 */	fmuld,s	%f62,%f58,%f60



/*     38 */	sxar2
/*     38 */	ldd	[%o1+104],%f366
/*     38 */	fmuld,s	%f72,%f72,%f74


/*     27 */	sxar2
/*     27 */	ldd	[%o1+16],%f114
/*     27 */	ldd	[%o1+112],%f370



/*     27 */	sxar2
/*     27 */	std,s	%f146,[%sp+2639]
/*     27 */	ldd	[%o1+24],%f118



/*     27 */	sxar2
/*     27 */	ldd	[%o1+120],%f374
/*     27 */	ldd	[%o1+32],%f124



/*     27 */	sxar2
/*     27 */	ldd	[%o1+128],%f380
/*     27 */	ldd	[%o1+40],%f130


/*     27 */	sxar2
/*     27 */	frcpad,s	%f60,%f64
/*     27 */	ldd	[%o1+136],%f386



/*     27 */	sxar2
/*     27 */	frcpad,s	%f74,%f78
/*     27 */	ldd	[%o1+72],%f150



/*     26 */	sxar2
/*     26 */	ldd	[%o1+168],%f406
/*     26 */	std,s	%f146,[%sp+2655]


/*     27 */	sxar2
/*     27 */	ldd	[%o1+80],%f168
/*     27 */	ldd	[%o1+176],%f424



/*     35 */	sxar2
/*     35 */	fnmsubd,s	%f60,%f64,%f66,%f60
/*     35 */	fnmsubd,s	%f74,%f78,%f66,%f74


/*     38 */	sxar2
/*     38 */	std,s	%f146,[%sp+2671]
/*     38 */	fmuld,s	%f60,%f60,%f68


/*     26 */	sxar2
/*     26 */	fmuld,s	%f74,%f74,%f80
/*     26 */	std,s	%f146,[%sp+2687]


/*     32 */	sxar2
/*     32 */	faddd,s	%f60,%f68,%f70
/*     32 */	fmaddd,s	%f68,%f68,%f60,%f68


/*     32 */	sxar2
/*     32 */	faddd,s	%f74,%f80,%f82
/*     32 */	fmaddd,s	%f80,%f80,%f74,%f80


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f70,%f64,%f64,%f70
/*     32 */	fmaddd,s	%f82,%f78,%f78,%f82


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f68,%f70,%f64,%f68
/*     32 */	fmaddd,s	%f80,%f82,%f78,%f80


/*     38 */	sxar2
/*     38 */	fmuld,s	%f68,%f68,%f76
/*     38 */	fmuld,s	%f84,%f80,%f84

/*     38 */	sxar1
/*     38 */	fmuld,s	%f76,%f76,%f76

.L5594:

/*    244 */	cmp	%o4,%g0

/*    244 */	ble	.L5598
	nop


.L5595:


/*    244 */	sxar2
/*    244 */	mov	%g0,%xg15
/*    244 */	fzero,s	%f212


/*    244 */	sxar2
/*    244 */	mov	%o4,%xg23
/*    244 */	fmovd,s	%f212,%f218


/*    244 */	sxar2
/*    244 */	fmovd,s	%f212,%f216
/*    244 */	fmovd,s	%f212,%f214

/*    244 */	sxar1
/*    244 */	fmovd,s	%f212,%f164

.L5596:


/*    311 */	sxar2
/*    311 */	add	%o3,%xg15,%xg16
/*    311 */	add	%xg15,104,%xg15


/*    219 */	sxar2
/*    219 */	ldd,s	[%xg16],%f108
/*    219 */	ldd,s	[%xg16+16],%f116


/*    219 */	sxar2
/*    219 */	subcc	%xg23,1,%xg23
/*    219 */	ldd,s	[%xg16+32],%f128


/*    231 */	sxar2
/*    231 */	ldd,s	[%xg16+64],%f148
/*    231 */	ldd,s	[%xg16+80],%f154


/*    219 */	sxar2
/*    219 */	ldd,s	[%xg16+48],%f158
/*    219 */	fmsubd,sc	%f108,%f66,%f106,%f104


/*    219 */	sxar2
/*    219 */	fmsubd,sc	%f364,%f66,%f110,%f108
/*    219 */	fmsubd,sc	%f116,%f66,%f114,%f112



/*    267 */	sxar2
/*    267 */	fmsubd,sc	%f128,%f66,%f124,%f126
/*    267 */	fmsubd,sc	%f372,%f66,%f118,%f116


/*    243 */	sxar2
/*    243 */	fmsubd,sc	%f384,%f66,%f130,%f128
/*    243 */	fmovd	%f148,%f196



/*    231 */	sxar2
/*    231 */	fmovd	%f148,%f452
/*    231 */	fmaddd,sc	%f154,%f66,%f150,%f152


/*    282 */	sxar2
/*    282 */	fmaddd,sc	%f158,%f66,%f72,%f156
/*    282 */	fmaddd,sc	%f410,%f66,%f168,%f154


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f104,%f104,%f120,%f122
/*     32 */	fmuld,s	%f126,%f108,%f126


/*    121 */	sxar2
/*    121 */	fmuld,s	%f196,%f196,%f196
/*    121 */	frcpad,s	%f156,%f162


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f108,%f108,%f122,%f122
/*     32 */	fmaddd,s	%f116,%f104,%f126,%f116


/*     32 */	sxar2
/*     32 */	fnmsubd,s	%f156,%f162,%f66,%f156
/*     32 */	fmaddd,s	%f112,%f112,%f122,%f122


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f128,%f112,%f116,%f128
/*     38 */	fmuld,s	%f156,%f156,%f166


/*     38 */	sxar2
/*     38 */	frsqrtad,s	%f122,%f132
/*     38 */	fmuld,s	%f122,%f34,%f134


/*     38 */	sxar2
/*     38 */	faddd,s	%f156,%f166,%f170
/*     38 */	fmaddd,s	%f166,%f166,%f156,%f166


/*     38 */	sxar2
/*     38 */	fmuld,s	%f132,%f132,%f136
/*     38 */	fmaddd,s	%f162,%f170,%f162,%f170


/*     38 */	sxar2
/*     38 */	fnmsubd,s	%f134,%f136,%f34,%f136
/*     38 */	fmaddd,s	%f166,%f170,%f162,%f166


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f132,%f136,%f132,%f132
/*     32 */	fmuld,s	%f132,%f132,%f138


/*     32 */	sxar2
/*     32 */	fnmsubd,s	%f134,%f138,%f34,%f138
/*     32 */	fmaddd,s	%f132,%f138,%f132,%f132


/*     32 */	sxar2
/*     32 */	fmuld,s	%f132,%f132,%f140
/*     32 */	fnmsubd,s	%f134,%f140,%f34,%f134


/*     38 */	sxar2
/*     38 */	fmaddd,s	%f132,%f134,%f132,%f132
/*     38 */	fmuld,s	%f132,%f122,%f122


/*     38 */	sxar2
/*     38 */	fmuld,s	%f132,%f128,%f142
/*     38 */	fmuld,s	%f122,%f68,%f144


/*    195 */	sxar2
/*    195 */	fmaddd,sc	%f148,%f122,%f146,%f122
/*    195 */	fmind,s	%f146,%f142,%f142


/*    192 */	sxar2
/*    192 */	fmaddd,sc	%f404,%f66,%f84,%f148
/*    192 */	fsubd,s	%f66,%f144,%f172


/*     32 */	sxar2
/*     32 */	fsubd,s	%f66,%f122,%f176
/*     32 */	fmaddd,s	%f186,%f144,%f182,%f184


/*     35 */	sxar2
/*     35 */	fmaddd,s	%f186,%f122,%f182,%f188
/*     35 */	fnmsubd,s	%f160,%f142,%f152,%f152


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f194,%f144,%f190,%f192
/*     32 */	fmaddd,s	%f194,%f122,%f190,%f202


/*    192 */	sxar2
/*    192 */	fmaxd,s	%f146,%f172,%f172
/*    192 */	fmaxd,s	%f146,%f176,%f176


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f144,%f184,%f204,%f184
/*     32 */	fmaddd,s	%f122,%f188,%f204,%f188


/*     32 */	sxar2
/*     32 */	fmaxd,s	%f164,%f152,%f164
/*     32 */	fmaddd,s	%f144,%f192,%f206,%f192


/*     32 */	sxar2
/*     32 */	fmuld,s	%f152,%f142,%f152
/*     32 */	fmaddd,s	%f122,%f202,%f206,%f202


/*     38 */	sxar2
/*     38 */	fmuld,s	%f172,%f172,%f174
/*     38 */	fmuld,s	%f176,%f176,%f178


/*     32 */	sxar2
/*     32 */	fmuld,s	%f154,%f152,%f154
/*     32 */	fmuld,s	%f76,%f174,%f180


/*     32 */	sxar2
/*     32 */	fmuld,s	%f178,%f178,%f198
/*     32 */	fmuld,s	%f172,%f174,%f200


/*     32 */	sxar2
/*     32 */	fmuld,s	%f178,%f196,%f178
/*     32 */	fmuld,s	%f196,%f176,%f196


/*     38 */	sxar2
/*     38 */	fmsubd,s	%f172,%f184,%f204,%f172
/*     38 */	fmuld,s	%f154,%f166,%f154


/*     32 */	sxar2
/*     32 */	fmsubd,s	%f176,%f188,%f204,%f176
/*     32 */	fmuld,s	%f174,%f180,%f174


/*     32 */	sxar2
/*     32 */	fmuld,s	%f196,%f198,%f196
/*     32 */	fnmsubd,s	%f144,%f192,%f172,%f144


/*     49 */	sxar2
/*     49 */	fnmsubd,s	%f122,%f202,%f176,%f122
/*     49 */	fnmsubd,s	%f208,%f154,%f84,%f210


/*     32 */	sxar2
/*     32 */	fnmsubd,s	%f34,%f154,%f148,%f154
/*     32 */	fmuld,s	%f200,%f174,%f200


/*     32 */	sxar2
/*     32 */	fmuld,s	%f178,%f196,%f178
/*     32 */	fmuld,s	%f200,%f144,%f200


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f178,%f122,%f200,%f178
/*     32 */	fmuld,s	%f178,%f132,%f178


/*     49 */	sxar2
/*     49 */	fmaddd,sc	%f414,%f178,%f146,%f158
/*     49 */	fmuld,s	%f158,%f128,%f128


/*     49 */	sxar2
/*     49 */	fmuld,s	%f158,%f154,%f158
/*     49 */	fmaddd,s	%f128,%f210,%f212,%f212


/*     49 */	sxar2
/*     49 */	fmaddd,s	%f158,%f104,%f214,%f214
/*     49 */	fmaddd,s	%f158,%f108,%f216,%f216

/*     49 */	sxar1
/*     49 */	fmaddd,s	%f158,%f112,%f218,%f218

/*    311 */	bne,pt	%icc, .L5596
	nop


.L5597:


/*    311 */	sxar2
/*    311 */	std,s	%f164,[%sp+2623]
/*    311 */	std,s	%f214,[%sp+2639]


/*    311 */	sxar2
/*    311 */	std,s	%f216,[%sp+2655]
/*    311 */	std,s	%f218,[%sp+2671]

/*    311 */	sxar1
/*    311 */	std,s	%f212,[%sp+2687]

.L5598:


/*     38 */	sxar2
/*     38 */	ldd,s	[%sp+2639],%f88
/*     38 */	fmuld,s	%f86,%f62,%f62

/*    120 */	sxar1
/*    120 */	ldd,s	[%sp+2623],%f92

/*    324 */	cmp	%g4,%o2


/*     57 */	sxar2
/*     57 */	ldd,s	[%sp+2655],%f94
/*     57 */	ldd,s	[%sp+2671],%f96


/*     57 */	sxar2
/*     57 */	ldd,s	[%sp+2687],%f98
/*     57 */	fmuld,s	%f32,%f88,%f88


/*     57 */	sxar2
/*     57 */	frcpad,s	%f92,%f90
/*     57 */	fmuld,s	%f32,%f94,%f94


/*     57 */	sxar2
/*     57 */	fmuld,s	%f32,%f96,%f96
/*     57 */	fmuld,s	%f32,%f98,%f98





/*    149 */	sxar2
/*    149 */	fnmsubd,s	%f92,%f90,%f66,%f92
/*    149 */	std	%f88,[%o5]



/*    149 */	sxar2
/*    149 */	std	%f94,[%o5+8]
/*    149 */	std	%f96,[%o5+16]


/*     38 */	sxar2
/*     38 */	std	%f98,[%o5+24]
/*     38 */	fmuld,s	%f92,%f92,%f100


/*     32 */	sxar2
/*     32 */	faddd,s	%f92,%f100,%f102
/*     32 */	fmaddd,s	%f100,%f100,%f92,%f100


/*     32 */	sxar2
/*     32 */	fmaddd,s	%f102,%f90,%f90,%f102
/*     32 */	fmaddd,s	%f100,%f102,%f90,%f100

/*     38 */	sxar1
/*     38 */	fmuld,s	%f62,%f100,%f62



/*    324 */	bge	.L5600
/*    324 */	std	%f62,[%o5+32]


.L5599:


/*    152 */	sxar2
/*    152 */	std	%f344,[%xg12]
/*    152 */	std	%f350,[%xg12+8]


/*    152 */	sxar2
/*    152 */	std	%f352,[%xg12+16]
/*    152 */	std	%f354,[%xg12+24]

/*    152 */	sxar1
/*    152 */	std	%f318,[%xg12+32]

.L5600:

/*    331 */	add	%g4,2,%g4

/*    331 */	sxar1
/*    331 */	add	%xg12,80,%xg12

/*    331 */	add	%o5,80,%o5


/*    331 */	subcc	%g3,1,%g3

/*    331 */	bne,pt	%icc, .L5593
/*    331 */	add	%o1,192,%o1


.L5601:


.L5603:

/*    332 */	sub	%sp,-3408,%sp
	retl
	nop


.LLFE3:
	.size	_ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE,.-_ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE
	.type	_ZN14CalcHydroForceclEPKN3EPI5HydroEiPKN3EPJ5HydroEiPN6RESULT5HydroE,#function
	.ident	"$Compiler: Fujitsu C/C++ Compiler Version 1.2.0 P-id: L30000-09 (Sep 17 2014 11:41:29) compile.C __sti___9_compile_C_76f9cf6c $"
	.section	".text"
	.global	__sti___9_compile_C_76f9cf6c
	.align	64
__sti___9_compile_C_76f9cf6c:
.LLFB4:
.L219:

/*      0 */	save	%sp,-208,%sp
.LLCFI0:


.L220:

/*      4 */	sethi	%h44(.LR0.cnt.7),%g1

/*      4 */	or	%g1,%m44(.LR0.cnt.7),%g1

/*      4 */	sllx	%g1,12,%g1

/*      4 */	or	%g1,%l44(.LR0.cnt.7),%g1

/*      4 */	ldd	[%g1],%f0


/*      4 */	call	atan
/*    ??? */	std	%f0,[%fp+2039]
/*      4 */	sethi	%h44(.LR0.cnt.12),%g2

/*      5 */	fzero	%f36
/*    ??? */	ldd	[%fp+2039],%f52

/*      4 */	or	%g2,%m44(.LR0.cnt.12),%g2
/*      4 */	sethi	%h44(_ZN4math2piE),%g4
/*      4 */	sllx	%g2,12,%g2

/*     17 */	sethi	%h44(_ZN6kernel2piE),%o0

/*      4 */	or	%g2,%l44(.LR0.cnt.12),%g2
/*      4 */	or	%g4,%m44(_ZN4math2piE),%g4
/*      4 */	ldd	[%g2],%f32
/*      4 */	sllx	%g4,12,%g4

/*     17 */	or	%o0,%m44(_ZN6kernel2piE),%o0

/*      4 */	or	%g4,%l44(_ZN4math2piE),%g4

/*     17 */	sllx	%o0,12,%o0

/*     19 */	sethi	%h44(.LR0.cnt.22),%g3

/*     17 */	or	%o0,%l44(_ZN6kernel2piE),%o0

/*      5 */	frcpad	%f36,%f34

/*     19 */	or	%g3,%m44(.LR0.cnt.22),%g3

/*      5 */	sethi	%h44(_ZN4math3NaNE),%g5

/*     19 */	sllx	%g3,12,%g3
/*     19 */	sethi	%h44(_ZN6kernel14normalzeFactorE),%o1
/*     19 */	or	%g3,%l44(.LR0.cnt.22),%g3

/*      5 */	or	%g5,%m44(_ZN4math3NaNE),%g5

/*      4 */	fmuld	%f0,%f32,%f0

/*      5 */	sllx	%g5,12,%g5

/*     19 */	ldd	[%g3],%f46
/*     19 */	or	%o1,%m44(_ZN6kernel14normalzeFactorE),%o1

/*      5 */	or	%g5,%l44(_ZN4math3NaNE),%g5

/*     19 */	sllx	%o1,12,%o1
/*     19 */	or	%o1,%l44(_ZN6kernel14normalzeFactorE),%o1

/*      5 */	fnmsubd	%f34,%f36,%f52,%f38
/*      5 */	fmuld	%f36,%f34,%f36

/*     19 */	frcpad	%f0,%f40

/*      4 */	std	%f0,[%g4]

/*     17 */	std	%f0,[%o0]

/*      5 */	fmaddd	%f38,%f38,%f38,%f42
/*      5 */	fmuld	%f38,%f38,%f44

/*     19 */	fnmsubd	%f40,%f0,%f52,%f0
/*     19 */	fmuld	%f46,%f40,%f46

/*      5 */	fmaddd	%f42,%f36,%f36,%f42
/*      5 */	fmaddd	%f44,%f44,%f38,%f44

/*     19 */	fmaddd	%f0,%f0,%f0,%f48
/*     19 */	fmuld	%f0,%f0,%f50

/*      5 */	fmaddd	%f44,%f42,%f36,%f44

/*     19 */	fmaddd	%f48,%f46,%f46,%f48
/*     19 */	fmaddd	%f50,%f50,%f0,%f50

/*      5 */	std	%f44,[%g5]

/*     19 */	fmaddd	%f50,%f48,%f46,%f50
/*     19 */	std	%f50,[%o1]
/*     19 */	ret
	restore


.L221:


.LLFE4:
	.size	__sti___9_compile_C_76f9cf6c,.-__sti___9_compile_C_76f9cf6c
	.type	__sti___9_compile_C_76f9cf6c,#function
	.section	".ctors",#alloc,#write
	.align	8
	.xword		__sti___9_compile_C_76f9cf6c
	.section	".rodata"
	.align	128
__fj_dlogexp1_const_:
	.word	0X3FF00000,0
	.word	0,0
	.word	0,0
	.word	0,0
	.word	0X3FEFF007,0XFC01FF00
	.word	0X3F5FF802,0XA9A00000
	.word	0X3D4621CC,0XF14F1D0B
	.word	0,0
	.word	0X3FEFE01F,0XE01FE020
	.word	0X3F6FF00A,0XA2B00000
	.word	0X3D20BC04,0XA086B56A
	.word	0,0
	.word	0X3FEFD047,0X94A10E6A
	.word	0X3F77EE11,0XEBD80000
	.word	0X3D0749D3,0XC2D23A07
	.word	0,0
	.word	0X3FEFC07F,0X1FC07F0
	.word	0X3F7FE02A,0X6B100000
	.word	0X3D19E23F,0XDDA40E4
	.word	0,0
	.word	0X3FEFB0C6,0X10D5E939
	.word	0X3F83E729,0X5D240000
	.word	0X3D4A7D8F,0X803597BB
	.word	0,0
	.word	0X3FEFA11C,0XAA01FA12
	.word	0X3F87DC47,0X5F800000
	.word	0X3D40A76D,0XD2512F06
	.word	0,0
	.word	0X3FEF9182,0XB6813BAF
	.word	0X3F8BCF71,0X2C740000
	.word	0X3D1C25E0,0X97BD9771
	.word	0,0
	.word	0X3FEF81F8,0X1F81F820
	.word	0X3F8FC0A8,0XB0FC0000
	.word	0X3CDF1E7C,0XF6D3A69C
	.word	0,0
	.word	0X3FEF727C,0XCE5F530A
	.word	0X3F91D7F7,0XEB9E0000
	.word	0X3D4D7CD8,0XAF80670B
	.word	0,0
	.word	0X3FEF6310,0XACA0DBB5
	.word	0X3F93CEA4,0X43460000
	.word	0X3D44AE9D,0XE694ADFB
	.word	0,0
	.word	0X3FEF53B3,0XA3FA204E
	.word	0X3F95C45A,0X51B80000
	.word	0X3D4A7112,0X77A49E0A
	.word	0,0
	.word	0X3FEF4465,0X9E4A4271
	.word	0X3F97B91B,0X7D60000
	.word	0XBD33B955,0XB602ACE4
	.word	0,0
	.word	0X3FEF3526,0X859B8CEC
	.word	0X3F99ACE7,0X551C0000
	.word	0X3D48A289,0XA04E0EFC
	.word	0,0
	.word	0X3FEF25F6,0X44230AB5
	.word	0X3F9B9FC0,0X27B00000
	.word	0XBD3B9A01,0XAE6922A
	.word	0,0
	.word	0X3FEF16D4,0XC4401F17
	.word	0X3F9D91A6,0X6C540000
	.word	0X3D2E61F1,0X658CFB9A
	.word	0,0
	.word	0X3FEF07C1,0XF07C1F08
	.word	0X3F9F829B,0XE780000
	.word	0X3D298026,0X7C7E09E4
	.word	0,0
	.word	0X3FEEF8BD,0XB389EBAD
	.word	0X3FA0B94F,0X7C190000
	.word	0X3D485D9D,0XA43F761F
	.word	0,0
	.word	0X3FEEE9C7,0XF8458E02
	.word	0X3FA1B0D9,0X89240000
	.word	0XBD33401E,0X9AE889BB
	.word	0,0
	.word	0X3FEEDAE0,0XA9B3D3A5
	.word	0X3FA2A7EC,0X22150000
	.word	0XBD278CE7,0X7A9163FE
	.word	0,0
	.word	0X3FEECC07,0XB301ECC0
	.word	0X3FA39E87,0XB9FF0000
	.word	0XBD40A815,0XBFA937F5
	.word	0,0
	.word	0X3FEEBD3C,0XFF850B0C
	.word	0X3FA494AC,0XC34E0000
	.word	0XBD4BB8E1,0XD6A40B6E
	.word	0,0
	.word	0X3FEEAE80,0X7ABA01EB
	.word	0X3FA58A5B,0XAFC90000
	.word	0XBD2B2B73,0X9570AD39
	.word	0,0
	.word	0X3FEE9FD2,0X1044E799
	.word	0X3FA67F94,0XF0950000
	.word	0XBD4099F0,0X60C0D895
	.word	0,0
	.word	0X3FEE9131,0XABF0B767
	.word	0X3FA77458,0XF6330000
	.word	0XBD3181DC,0XE586AF09
	.word	0,0
	.word	0X3FEE829F,0X39AEF509
	.word	0X3FA868A8,0X30840000
	.word	0XBD12623A,0X134AC693
	.word	0,0
	.word	0X3FEE741A,0XA59750E4
	.word	0X3FA95C83,0XEC90000
	.word	0XBD2C1482,0X97C5FEB8
	.word	0,0
	.word	0X3FEE65A3,0XDBE74D6B
	.word	0X3FAA4FE9,0XFFA40000
	.word	0XBD36E584,0XA0402925
	.word	0,0
	.word	0X3FEE573A,0XC901E574
	.word	0X3FAB42DD,0X71190000
	.word	0X3D4C6FB0,0XA34531F6
	.word	0,0
	.word	0X3FEE48DF,0X596F3394
	.word	0X3FAC355D,0XD0920000
	.word	0X3D2F2CCC,0X9ABF8388
	.word	0,0
	.word	0X3FEE3A91,0X79DC1A73
	.word	0X3FAD276B,0X8ADB0000
	.word	0X3D16A423,0XC78A64B0
	.word	0,0
	.word	0X3FEE2C51,0X1719EE16
	.word	0X3FAE1907,0XC270000
	.word	0X3D48056E,0X66E7585F
	.word	0,0
	.word	0X3FEE1E1E,0X1E1E1E1E
	.word	0X3FAF0A30,0XC0110000
	.word	0X3D48A998,0X5F325C5C
	.word	0,0
	.word	0X3FEE0FF8,0X7C01E100
	.word	0X3FAFFAE9,0X119C0000
	.word	0XBD4B3F32,0X2F674EAB
	.word	0,0
	.word	0X3FEE01E0,0X1E01E01E
	.word	0X3FB07598,0X35990000
	.word	0XBD3B8ECF,0XE4B59987
	.word	0,0
	.word	0X3FEDF3D4,0XF17DE4DB
	.word	0X3FB0ED83,0X9B550000
	.word	0X3D437F05,0XC95BAA26
	.word	0,0
	.word	0X3FEDE5D6,0XE3F8868A
	.word	0X3FB16536,0XEEA38000
	.word	0XBD147C5E,0X768FA309
	.word	0,0
	.word	0X3FEDD7E5,0XE316D94C
	.word	0X3FB1DCB2,0X63DB0000
	.word	0X3D39444F,0X5E9E8981
	.word	0,0
	.word	0X3FEDCA01,0XDCA01DCA
	.word	0X3FB253F6,0X2F0A0000
	.word	0X3D3416F8,0XFB69A701
	.word	0,0
	.word	0X3FEDBC2A,0XBE7D71D4
	.word	0X3FB2CB02,0X83F60000
	.word	0XBD40F08E,0X9ACD47EF
	.word	0,0
	.word	0X3FEDAE60,0X76B981DB
	.word	0X3FB341D7,0X961C0000
	.word	0XBD4717B6,0XB33E44F8
	.word	0,0
	.word	0X3FEDA0A2,0XF3803B41
	.word	0X3FB3B875,0X98B18000
	.word	0X3D4B76FA,0X9AD4D73B
	.word	0,0
	.word	0X3FED92F2,0X231E7F8A
	.word	0X3FB42EDC,0XBEA68000
	.word	0XBD4C87E2,0X22B06CA6
	.word	0,0
	.word	0X3FED854D,0XF401D855
	.word	0X3FB4A50D,0X3AA18000
	.word	0X3D482007,0XB3DDA307
	.word	0,0
	.word	0X3FED77B6,0X54B82C34
	.word	0X3FB51B07,0X3F060000
	.word	0X3D383F69,0X278E686A
	.word	0,0
	.word	0X3FED6A2B,0X33EF7448
	.word	0X3FB590CA,0XFDF00000
	.word	0X3D3C284F,0X5722ABAA
	.word	0,0
	.word	0X3FED5CAC,0X807572B2
	.word	0X3FB60658,0XA9378000
	.word	0XBD479E27,0X108B1D84
	.word	0,0
	.word	0X3FED4F3A,0X293769CA
	.word	0X3FB67BB0,0X726F0000
	.word	0XBD4F8236,0XD258429C
	.word	0,0
	.word	0X3FED41D4,0X1D41D41D
	.word	0X3FB6F0D2,0X8AE58000
	.word	0XBD34B464,0X1B664613
	.word	0,0
	.word	0X3FED347A,0X4BC01D34
	.word	0X3FB765BF,0X23A68000
	.word	0X3D4F09A1,0XFE51DED5
	.word	0,0
	.word	0X3FED272C,0XA3FC5B1A
	.word	0X3FB7DA76,0X6D7B0000
	.word	0X3D32CC84,0X4480C89B
	.word	0,0
	.word	0X3FED19EB,0X155F08A4
	.word	0X3FB84EF8,0X98E80000
	.word	0X3D441519,0X6DCB441C
	.word	0,0
	.word	0X3FED0CB5,0X8F6EC074
	.word	0X3FB8C345,0XD6318000
	.word	0X3D3B20F5,0XACB42A66
	.word	0,0
	.word	0X3FECFF8C,0X1CFF8C0
	.word	0X3FB9375E,0X55598000
	.word	0XBD40911E,0X463F9E4E
	.word	0,0
	.word	0X3FECF26E,0X5C44BFC6
	.word	0X3FB9AB42,0X46200000
	.word	0X3D49D66D,0XF661E3E8
	.word	0,0
	.word	0X3FECE55C,0X8EAC7900
	.word	0X3FBA1EF1,0XD8060000
	.word	0X3D3CD417,0X6DF97BCB
	.word	0,0
	.word	0X3FECD856,0X89039B0B
	.word	0X3FBA926D,0X3A4B0000
	.word	0XBD454E4D,0X7A16EAB2
	.word	0,0
	.word	0X3FECCB5C,0X3B636E3A
	.word	0X3FBB05B4,0X9BEE8000
	.word	0XBD4E00DD,0X3E707B5A
	.word	0,0
	.word	0X3FECBE6D,0X9601CBE7
	.word	0X3FBB78C8,0X2BB10000
	.word	0XBD325EF7,0XBC3987E7
	.word	0,0
	.word	0X3FECB18A,0X8930DE60
	.word	0X3FBBEBA8,0X18148000
	.word	0XBD389B78,0XB6DF1F57
	.word	0,0
	.word	0X3FECA4B3,0X55EE191
	.word	0X3FBC5E54,0X8F5C0000
	.word	0XBD4C5E75,0X14F4083F
	.word	0,0
	.word	0X3FEC97E6,0XFB15E44D
	.word	0X3FBCD0CD,0XBF8C0000
	.word	0X3D33E14D,0XB50DD743
	.word	0,0
	.word	0X3FEC8B26,0X5AFB8A42
	.word	0X3FBD4313,0XD66C8000
	.word	0X3D49AEAF,0X21BB2A3B
	.word	0,0
	.word	0X3FEC7E71,0X15D0CE95
	.word	0X3FBDB527,0X1880000
	.word	0XBD436C43,0XD4A8F3F0
	.word	0,0
	.word	0X3FEC71C7,0X1C71C71C
	.word	0X3FBE2707,0X6E2B0000
	.word	0XBD2A342C,0X2AF0003C
	.word	0,0
	.word	0X3FEC6528,0X5FD56843
	.word	0X3FBE98B5,0X49670000
	.word	0X3D346774,0X89C50E97
	.word	0,0
	.word	0X3FEC5894,0XD10D4986
	.word	0X3FBF0A30,0XC0118000
	.word	0XBD3D599E,0X83368E91
	.word	0,0
	.word	0X3FEC4C0C,0X61456A8E
	.word	0X3FBF7B79,0XFEC38000
	.word	0XBD010987,0XE897ED01
	.word	0,0
	.word	0X3FEC3F8F,0X1C3F8F0
	.word	0X3FBFEC91,0X31DC0000
	.word	0XBD354555,0XD1AE6607
	.word	0,0
	.word	0X3FEC331C,0XA3E91679
	.word	0X3FC02EBB,0X42BF4000
	.word	0XBD15A8FA,0X5CE00E5D
	.word	0,0
	.word	0X3FEC26B5,0X392EA01C
	.word	0X3FC06715,0X12CA4000
	.word	0X3D496E2A,0X18C8FD71
	.word	0,0
	.word	0X3FEC1A58,0XB327F576
	.word	0X3FC09F56,0X1EE70000
	.word	0X3D49C33E,0XA3AA0B96
	.word	0,0
	.word	0X3FEC0E07,0X381C0E0
	.word	0X3FC0D77E,0X7CD08000
	.word	0X3D3CB2CD,0X2EE2F482
	.word	0,0
	.word	0X3FEC01C0,0X1C01C01C
	.word	0X3FC10F8E,0X42254000
	.word	0XBD293B38,0X43396307
	.word	0,0
	.word	0X3FEBF583,0XEE868D8B
	.word	0X3FC14785,0X84674000
	.word	0X3D156345,0X1027C750
	.word	0,0
	.word	0X3FEBE952,0X6D0769FA
	.word	0X3FC17F64,0X58FCC000
	.word	0XBD49EF01,0X4BDB0DC9
	.word	0,0
	.word	0X3FEBDD2B,0X899406F7
	.word	0X3FC1B72A,0XD52F8000
	.word	0XBD485FD6,0XF9FB971A
	.word	0,0
	.word	0X3FEBD10F,0X365451B6
	.word	0X3FC1EED9,0XE2DC000
	.word	0X3D161563,0X7097648F
	.word	0,0
	.word	0X3FEBC4FD,0X65883E7B
	.word	0X3FC2266F,0X190A4000
	.word	0X3D4ACB7D,0X51EFC602
	.word	0,0
	.word	0X3FEBB8F6,0X9879493
	.word	0X3FC25DED,0XABC8000
	.word	0XBD452E3D,0X58744DDC
	.word	0,0
	.word	0X3FEBACF9,0X14C1BAD0
	.word	0X3FC29552,0XF8200000
	.word	0XBD35B967,0XF4471DFC
	.word	0,0
	.word	0X3FEBA106,0X79BD8488
	.word	0X3FC2CCA0,0XF5F60000
	.word	0XBD3B5EF1,0X91AFF120
	.word	0,0
	.word	0X3FEB951E,0X2B18FF23
	.word	0X3FC303D7,0X18E48000
	.word	0XBCD680B5,0XCE3ECB05
	.word	0,0
	.word	0X3FEB8940,0X1B89401C
	.word	0X3FC33AF5,0X75770000
	.word	0X3D3C9ECC,0XA2FE72A5
	.word	0,0
	.word	0X3FEB7D6C,0X3DDA338B
	.word	0X3FC371FC,0X201E8000
	.word	0X3D3EE877,0X9B2D8ABC
	.word	0,0
	.word	0X3FEB71A2,0X84EE6B34
	.word	0X3FC3A8EB,0X2D31C000
	.word	0XBD4C8A09,0X105455F8
	.word	0,0
	.word	0X3FEB65E2,0XE3BEEE05
	.word	0X3FC3DFC2,0XB0ECC000
	.word	0X3D28A72A,0X62B8C13F
	.word	0,0
	.word	0X3FEB5A2D,0X4D5B081F
	.word	0X3FC41682,0XBF728000
	.word	0XBD210047,0X81F849D
	.word	0,0
	.word	0X3FEB4E81,0XB4E81B4F
	.word	0X3FC44D2B,0X6CCB8000
	.word	0XBD170CC1,0X6135783C
	.word	0,0
	.word	0X3FEB42E0,0XDA17007
	.word	0X3FC483BC,0XCCE70000
	.word	0XBD4C22B5,0XB1B81393
	.word	0,0
	.word	0X3FEB3748,0X4AD806CE
	.word	0X3FC4BA36,0XF39A4000
	.word	0X3D45E55A,0X2606F30E
	.word	0,0
	.word	0X3FEB2BBA,0X5FF26A23
	.word	0X3FC4F099,0XF4A24000
	.word	0XBD3E9BF2,0XFAFEAF27
	.word	0,0
	.word	0X3FEB2036,0X406C80D9
	.word	0X3FC526E5,0XE3A1C000
	.word	0XBD3790BA,0X37FC5238
	.word	0,0
	.word	0X3FEB14BB,0XDFD760E6
	.word	0X3FC55D1A,0XD4234000
	.word	0XBD429135,0X912CDC0C
	.word	0,0
	.word	0X3FEB094B,0X31D922A4
	.word	0X3FC59338,0XD9984000
	.word	0XBD4F7A2C,0XBA455516
	.word	0,0
	.word	0X3FEAFDE4,0X2A2CB482
	.word	0X3FC5C940,0X7598000
	.word	0XBD3A8D94,0X8CD23322
	.word	0,0
	.word	0X3FEAF286,0XBCA1AF28
	.word	0X3FC5FF30,0X70A78000
	.word	0X3D43D3C8,0X73E20A07
	.word	0,0
	.word	0X3FEAE732,0XDD1C2A09
	.word	0X3FC6350A,0X28AAC000
	.word	0XBD48A84F,0XD51FA714
	.word	0,0
	.word	0X3FEADBE8,0X7F94905E
	.word	0X3FC66ACD,0X4272C000
	.word	0XBD42AF21,0X201C9C3D
	.word	0,0
	.word	0X3FEAD0A7,0X98177693
	.word	0X3FC6A079,0XD0F7C000
	.word	0XBD452E03,0XDDB97585
	.word	0,0
	.word	0X3FEAC570,0X1AC5701B
	.word	0X3FC6D60F,0XE719C000
	.word	0X3D421C8D,0X54765C4D
	.word	0,0
	.word	0X3FEABA41,0XFBD2E5B1
	.word	0X3FC70B8F,0X97A1C000
	.word	0XBD458ACD,0X84AB5209
	.word	0,0
	.word	0X3FEAAF1D,0X2F87EBFD
	.word	0X3FC740F8,0XF5404000
	.word	0XBD30B66C,0X99018AA1
	.word	0,0
	.word	0X3FEAA401,0XAA401AA4
	.word	0X3FC7764C,0X128F4000
	.word	0XBD4ED8B6,0XFCB861C3
	.word	0,0
	.word	0X3FEA98EF,0X606A63BE
	.word	0X3FC7AB89,0X210C000
	.word	0X3D49091B,0XE36B2D6A
	.word	0,0
	.word	0X3FEA8DE6,0X4688EBAB
	.word	0X3FC7E0AF,0XD630C000
	.word	0X3D139E7C,0X1D8F1034
	.word	0,0
	.word	0X3FEA82E6,0X5130E159
	.word	0X3FC815C0,0XA1434000
	.word	0X3D47EAD6,0X836FF18C
	.word	0,0
	.word	0X3FEA77EF,0X750A56DA
	.word	0X3FC84ABB,0X75864000
	.word	0X3D41392A,0X9058EA17
	.word	0,0
	.word	0X3FEA6D01,0XA6D01A6D
	.word	0X3FC87FA0,0X6520C000
	.word	0X3D322120,0X401202FC
	.word	0,0
	.word	0X3FEA621C,0XDB4F8FDF
	.word	0X3FC8B46F,0X82238000
	.word	0XBD4DA4C1,0XBDFA4507
	.word	0,0
	.word	0X3FEA5741,0X7688A4A
	.word	0X3FC8E928,0XDE888000
	.word	0XBD42BF55,0XA7614696
	.word	0,0
	.word	0X3FEA4C6E,0X200D2637
	.word	0X3FC91DCC,0X8C340000
	.word	0X3D37BC6A,0XBDDEFF46
	.word	0,0
	.word	0X3FEA41A4,0X1A41A41A
	.word	0X3FC9525A,0X9CF44000
	.word	0X3D46B476,0X41307539
	.word	0,0
	.word	0X3FEA36E2,0XEB1C432D
	.word	0X3FC986D3,0X22818000
	.word	0X3CF93B56,0X4DD44000
	.word	0,0
	.word	0X3FEA2C2A,0X87C51CA0
	.word	0X3FC9BB36,0X2E7E0000
	.word	0XBD21F2A8,0XA1CE0FFC
	.word	0,0
	.word	0X3FEA217A,0XE575FF2F
	.word	0X3FC9EF83,0XD2768000
	.word	0X3D4A33D2,0X13018C80
	.word	0,0
	.word	0X3FEA16D3,0XF97A4B02
	.word	0X3FCA23BC,0X1FE2C000
	.word	0XBD3539CD,0X91DC9F0B
	.word	0,0
	.word	0X3FEA0C35,0XB92ECDF1
	.word	0X3FCA57DF,0X28244000
	.word	0X3D3B99C8,0XCA1D9ABB
	.word	0,0
	.word	0X3FEA01A0,0X1A01A01A
	.word	0X3FCA8BEC,0XFC884000
	.word	0XBD40E73D,0X186F2318
	.word	0,0
	.word	0X3FE9F713,0X117200D0
	.word	0X3FCABFE5,0XAE460000
	.word	0X3D424B85,0X63507632
	.word	0,0
	.word	0X3FE9EC8E,0X951033D9
	.word	0X3FCAF3C9,0X4E80C000
	.word	0XBCBA4E63,0X3FCD9066
	.word	0,0
	.word	0X3FE9E212,0X9A7D5F0A
	.word	0X3FCB2797,0XEE464000
	.word	0XBD3BE88A,0X906D00A9
	.word	0,0
	.word	0X3FE9D79F,0X176B682D
	.word	0X3FCB5B51,0X9E8FC000
	.word	0XBD34B722,0XEC011F31
	.word	0,0
	.word	0X3FE9CD34,0X19CD340
	.word	0X3FCB8EF6,0X70420000
	.word	0X3D387533,0X321788E0
	.word	0,0
	.word	0X3FE9C2D1,0X4EE4A102
	.word	0X3FCBC286,0X742D8000
	.word	0X3D39AC53,0XF39D121C
	.word	0,0
	.word	0X3FE9B876,0XF5262DD1
	.word	0X3FCBF601,0XBB0E4000
	.word	0X3D2386A9,0X47C378B5
	.word	0,0
	.word	0X3FE9AE24,0XEA5510DA
	.word	0X3FCC2968,0X558C0000
	.word	0X3D48C0A3,0X8471D70
	.word	0,0
	.word	0X3FE9A3DB,0X2474FB98
	.word	0X3FCC5CBA,0X543B0000
	.word	0XBD4BDB58,0X84D2EAC1
	.word	0,0
	.word	0X3FE99999,0X9999999A
	.word	0X3FCC8FF7,0XC79A8000
	.word	0X3D4A21AC,0X25D81EF3
	.word	0,0
	.word	0X3FE98F60,0X3FE670A0
	.word	0X3FCCC320,0XC0178000
	.word	0XBD4AFDBF,0X1966B21B
	.word	0,0
	.word	0X3FE9852F,0XD8EC0FF
	.word	0X3FCCF635,0X4E09C000
	.word	0X3D277123,0X9A07D55B
	.word	0,0
	.word	0X3FE97B05,0XF8D56652
	.word	0X3FCD2935,0X81B6C000
	.word	0XBD383270,0X128AAA5F
	.word	0,0
	.word	0X3FE970E4,0XF80CB872
	.word	0X3FCD5C21,0X6B4FC000
	.word	0XBD21BA91,0XBBCA681B
	.word	0,0
	.word	0X3FE966CC,0X1966CC0
	.word	0X3FCD8EF9,0X1AF30000
	.word	0X3D4D5DF4,0X7936710
	.word	0,0
	.word	0X3FE95CBB,0XBE377AE
	.word	0X3FCDC1BC,0XA0AC0000
	.word	0XBD43829F,0X2CEB999D
	.word	0,0
	.word	0X3FE952B2,0XD73EE97
	.word	0X3FCDF46C,0XC724000
	.word	0XBD42D0BE,0XA7A437E3
	.word	0,0
	.word	0X3FE948B0,0XFCD6E9E0
	.word	0X3FCE2707,0X6E2B0000
	.word	0XBD3A342C,0X2AF0003C
	.word	0,0
	.word	0X3FE93EB7,0XD0AA6759
	.word	0X3FCE598E,0XD5A88000
	.word	0XBD0D134B,0XCF1E98A1
	.word	0,0
	.word	0X3FE934C6,0X7F9B2CE6
	.word	0X3FCE8C02,0X52AA4000
	.word	0X3D4A5FE9,0X1FC5C640
	.word	0,0
	.word	0X3FE92ADD,0X64AB74
	.word	0X3FCEBE61,0XF4DD8000
	.word	0XBD23D453,0X30FDCA4D
	.word	0,0
	.word	0X3FE920FB,0X49D0E229
	.word	0X3FCEF0AD,0XCBDC4000
	.word	0X3D493652,0X18DE5437
	.word	0,0
	.word	0X3FE91721,0X52B841DD
	.word	0X3FCF22E5,0XE72F0000
	.word	0X3D405D58,0X75F40DF0
	.word	0,0
	.word	0X3FE90D4F,0X120190D5
	.word	0X3FCF550A,0X564B8000
	.word	0XBD2323E3,0XA09202FE
	.word	0,0
	.word	0X3FE90384,0X7EA1CEC1
	.word	0X3FCF871B,0X28954000
	.word	0X3D404501,0X4AD6C881
	.word	0,0
	.word	0X3FE8F9C1,0X8F9C18FA
	.word	0X3FCFB918,0X6D5E4000
	.word	0XBD0D572A,0XAB993C87
	.word	0,0
	.word	0X3FE8F006,0X3C018F00
	.word	0X3FCFEB02,0X33E60000
	.word	0X3D2F316E,0X32D5E8C7
	.word	0,0
	.word	0X3FE8E652,0X7AF1373F
	.word	0X3FD00E6C,0X45AD6000
	.word	0XBD4FC672,0XE55A3FDC
	.word	0,0
	.word	0X3FE8DCA6,0X4397E408
	.word	0X3FD0274D,0XC16C2000
	.word	0X3D2979E8,0X9CF835C2
	.word	0,0
	.word	0X3FE8D301,0X8D3018D3
	.word	0X3FD04025,0X94B4E000
	.word	0XBD4F7E4A,0X3B085E94
	.word	0,0
	.word	0X3FE8C964,0X4F01EFBC
	.word	0X3FD058F3,0XC703E000
	.word	0X3D478BCC,0XA196E4A9
	.word	0,0
	.word	0X3FE8BFCE,0X8062FF3A
	.word	0X3FD071B8,0X5FCD6000
	.word	0XBD3BCB8B,0XA3E01A11
	.word	0,0
	.word	0X3FE8B640,0X18B64019
	.word	0X3FD08A73,0X667C6000
	.word	0XBD40A1F1,0X5F9D2E6C
	.word	0,0
	.word	0X3FE8ACB9,0XF6BF3AA
	.word	0X3FD0A324,0XE273A000
	.word	0XBD4E3941,0X1810BFCF
	.word	0,0
	.word	0X3FE8A339,0X5C018A34
	.word	0X3FD0BBCC,0XDB0D2000
	.word	0X3D32F32C,0XCC5DCDFB
	.word	0,0
	.word	0X3FE899C0,0XF601899C
	.word	0X3FD0D46B,0X579AC000
	.word	0XBD4169BF,0X4DF8F0D
	.word	0,0
	.word	0X3FE8904F,0XD503744B
	.word	0X3FD0ED00,0X5F658000
	.word	0XBD22DC75,0X285AA803
	.word	0,0
	.word	0X3FE886E5,0XF0ABB04A
	.word	0X3FD1058B,0XF9AE4000
	.word	0X3D45AA31,0X3F415699
	.word	0,0
	.word	0X3FE87D83,0X40AB6E97
	.word	0X3FD11E0E,0X2DADA000
	.word	0XBD2A47F8,0X8FCCE5BA
	.word	0,0
	.word	0X3FE87427,0XBCC092B9
	.word	0X3FD13687,0X293A000
	.word	0X3D4160BD,0XB314C76F
	.word	0,0
	.word	0X3FE86AD3,0X5CB59A84
	.word	0X3FD14EF6,0X7F886000
	.word	0X3D40B42C,0XD101B436
	.word	0,0
	.word	0X3FE86186,0X18618618
	.word	0X3FD1675C,0XABABA000
	.word	0X3D38380E,0X731F55C4
	.word	0,0
	.word	0X3FE8583F,0XE7A7C018
	.word	0X3FD17FB9,0X8E150000
	.word	0X3D42BABA,0X2C5AECBE
	.word	0,0
	.word	0X3FE84F00,0XC2780614
	.word	0X3FD1980D,0X2DD42000
	.word	0X3D2B7B3A,0X7A361C9A
	.word	0,0
	.word	0X3FE845C8,0XA0CE5129
	.word	0X3FD1B057,0X91F08000
	.word	0XBD32DD46,0X6DC55E2D
	.word	0,0
	.word	0X3FE83C97,0X7AB2BEDD
	.word	0X3FD1C898,0XC169A000
	.word	0XBD381410,0XE5C62AFF
	.word	0,0
	.word	0X3FE8336D,0X48397A24
	.word	0X3FD1E0D0,0XC3372000
	.word	0XBD428386,0XAB2795B0
	.word	0,0
	.word	0X3FE82A4A,0X182A4A0
	.word	0X3FD1F8FF,0X9E48A000
	.word	0X3D27946C,0X40CBE77
	.word	0,0
	.word	0X3FE8212D,0X9EBA4018
	.word	0X3FD21125,0X59862000
	.word	0XBD43E8C3,0XA2EB579E
	.word	0,0
	.word	0X3FE81818,0X18181818
	.word	0X3FD22941,0XFBCF8000
	.word	0XBD3A6976,0XF5EB0963
	.word	0,0
	.word	0X3FE80F09,0X65DFABCB
	.word	0X3FD24155,0X8BFD2000
	.word	0XBD47F800,0X66EB81A9
	.word	0,0
	.word	0X3FE80601,0X80601806
	.word	0X3FD25960,0X10DF8000
	.word	0XBD438C21,0XEED8AE0F
	.word	0,0
	.word	0X3FE7FD00,0X5FF40180
	.word	0X3FD27161,0X913F8000
	.word	0X3D34F4F1,0XF61564B4
	.word	0,0
	.word	0X3FE7F405,0XFD017F40
	.word	0X3FD2895A,0X13DE8000
	.word	0X3D3A8D7A,0XD24C13F0
	.word	0,0
	.word	0X3FE7EB12,0X4FFA053B
	.word	0X3FD2A149,0X9F762000
	.word	0X3D479220,0X57062A92
	.word	0,0
	.word	0X3FE7E225,0X515A4F1D
	.word	0X3FD2B930,0X3AB8A000
	.word	0XBD26DB12,0XD6BFB0A5
	.word	0,0
	.word	0X3FE7D93E,0XF9AA4B46
	.word	0X3FD2D10D,0XEC508000
	.word	0X3D360C61,0XF7088353
	.word	0,0
	.word	0X3FE7D05F,0X417D05F4
	.word	0X3FD2E8E2,0XBAE12000
	.word	0XBD267B1E,0X99B72BD8
	.word	0,0
	.word	0X3FE7C786,0X2170949F
	.word	0X3FD300AE,0XAD064000
	.word	0XBD45E854,0XBA4501AA
	.word	0,0
	.word	0X3FE7BEB3,0X922E017C
	.word	0X3FD31871,0XC9544000
	.word	0X3D184FAB,0X94CECFD9
	.word	0,0
	.word	0X3FE7B5E7,0X8C693733
	.word	0X3FD3302C,0X16586000
	.word	0X3D36217D,0XC2A3E08B
	.word	0,0
	.word	0X3FE7AD22,0X8E0ECC3
	.word	0X3FD347DD,0X9A988000
	.word	0XBD25594D,0XD4C58092
	.word	0,0
	.word	0X3FE7A463,0X5E918C
	.word	0X3FD35F86,0X5C932000
	.word	0X3D427C10,0XD8AF2D5B
	.word	0,0
	.word	0X3FE79BAA,0X6BB6398B
	.word	0X3FD37726,0X62BFE000
	.word	0XBD3E9436,0XAC53B023
	.word	0,0
	.word	0X3FE792F8,0X43C689C3
	.word	0X3FD38EBD,0XB38EE000
	.word	0XBD49BE9F,0X4660ACD8
	.word	0,0
	.word	0X3FE78A4C,0X8178A4C8
	.word	0X3FD3A64C,0X55694000
	.word	0X3D37A71C,0XBCD735D0
	.word	0,0
	.word	0X3FE781A7,0X1DC01782
	.word	0X3FD3BDD2,0X4EB14000
	.word	0X3D46D425,0XB478C893
	.word	0,0
	.word	0X3FE77908,0X119AC60D
	.word	0X3FD3D54F,0XA5C20000
	.word	0XBD41E0F1,0X932E350E
	.word	0,0
	.word	0X3FE7706F,0X5610D8D0
	.word	0X3FD3ECC4,0X60EF6000
	.word	0XBD060286,0X27C1300F
	.word	0,0
	.word	0X3FE767DC,0XE434A9B1
	.word	0X3FD40430,0X8686A000
	.word	0X3D3F8EF4,0X3049F7D3
	.word	0,0
	.word	0X3FE75F50,0XB522B17C
	.word	0X3FD41B94,0X1CCE0000
	.word	0X3D47DCB7,0XF60DE01C
	.word	0,0
	.word	0X3FE756CA,0XC201756D
	.word	0X3FD432EF,0X2A04E000
	.word	0X3D40276B,0X3674752A
	.word	0,0
	.word	0X3FE74E4B,0X40174E5
	.word	0X3FD44A41,0XB463C000
	.word	0X3D31EE28,0XF37CF612
	.word	0,0
	.word	0X3FE745D1,0X745D1746
	.word	0X3FD4618B,0XC21C6000
	.word	0XBD13D82F,0X484C84CC
	.word	0,0
	.word	0X3FE73D5E,0XC5899F7
	.word	0X3FD478CD,0X5959C000
	.word	0XBD484DD9,0X3CDCBDA
	.word	0,0
	.word	0X3FE734F0,0XC541FE8D
	.word	0X3FD49006,0X80400000
	.word	0X3D43A198,0XF2F83A
	.word	0,0
	.word	0X3FE72C89,0X9870F91F
	.word	0X3FD4A737,0X3CED0000
	.word	0XBD39A234,0XEBF35449
	.word	0,0
	.word	0X3FE72428,0X7F46DEBC
	.word	0X3FD4BE5F,0X95778000
	.word	0XBD3D7C92,0XCD9AD824
	.word	0,0
	.word	0X3FE71BCD,0X732E940A
	.word	0X3FD4D57F,0X8FEFE000
	.word	0X3D23F926,0X7FD06868
	.word	0,0
	.word	0X3FE71378,0X6D9C7C09
	.word	0X3FD4EC97,0X32600000
	.word	0X3D234D7A,0XAF04D104
	.word	0,0
	.word	0X3FE70B29,0X680E66FA
	.word	0X3FD503A6,0X82CB2000
	.word	0XBD2A68C8,0XF16F9B5D
	.word	0,0
	.word	0X3FE702E0,0X5C0B8170
	.word	0X3FD51AAD,0X872E0000
	.word	0XBD3F4BD8,0XDB0A7CC1
	.word	0,0
	.word	0X3FE6FA9D,0X43244380
	.word	0X3FD531AC,0X457EE000
	.word	0X3D3DF83B,0X7D931501
	.word	0,0
	.word	0X3FE6F260,0X16F26017
	.word	0X3FD548A2,0XC3ADE000
	.word	0XBD4B3A60,0X673DF8C2
	.word	0,0
	.word	0X3FE6EA28,0XD118B474
	.word	0X3FD55F91,0X7A44000
	.word	0XBD11E647,0X78DF4A62
	.word	0,0
	.word	0X3FE6E1F7,0X6B4337C7
	.word	0X3FD57677,0X17456000
	.word	0XBD364EAD,0X9524D7CA
	.word	0,0
	.word	0X3FE6D9CB,0XDF26EAEF
	.word	0X3FD58D54,0XF86E0000
	.word	0X3D2791F3,0XA795215
	.word	0,0
	.word	0X3FE6D1A6,0X2681C861
	.word	0X3FD5A42A,0XB0F4C000
	.word	0X3D4FC338,0XA1A4108B
	.word	0,0
	.word	0X3FE6C986,0X3B1AB429
	.word	0X3FD5BAF8,0X46AA2000
	.word	0XBD339AE8,0XF873FA41
	.word	0,0
	.word	0X3FE6C16C,0X16C16C17
	.word	0X3FD5D1BD,0XBF580000
	.word	0X3D4394A1,0X1B1C1EE4
	.word	0,0
	.word	0X3FE6B957,0XB34E7803
	.word	0X3FD5E87B,0X20C2A000
	.word	0XBD456C17,0X38446383
	.word	0,0
	.word	0X3FE6B149,0XAA31A3D
	.word	0X3FD5FF30,0X70A7A000
	.word	0XBD48586F,0X183BEBF2
	.word	0,0
	.word	0X3FE6A940,0X16A94017
	.word	0X3FD615DD,0XB4BEC000
	.word	0X3D13C7CA,0X90BC04B2
	.word	0,0
	.word	0X3FE6A13C,0XD1537290
	.word	0X3FD62C82,0XF2B9C000
	.word	0X3D3E54BD,0XBD7C8A98
	.word	0,0
	.word	0X3FE6993F,0X349CC726
	.word	0X3FD64320,0X30444000
	.word	0X3D3EFE02,0X7A01D7DF
	.word	0,0
	.word	0X3FE69147,0X3A88D0C0
	.word	0X3FD659B5,0X7303E000
	.word	0X3D1F281D,0XB0AF8EFC
	.word	0,0
	.word	0X3FE68954,0XDD2390BA
	.word	0X3FD67042,0XC0984000
	.word	0XBD1CF5B9,0X2118779C
	.word	0,0
	.word	0X3FE68168,0X16816817
	.word	0X3FD686C8,0X1E9B2000
	.word	0XBD46A277,0X7A83DFD6
	.word	0,0
	.word	0X3FE67980,0XE0BF08C7
	.word	0X3FD69D45,0X92A04000
	.word	0XBD43A721,0X2034002D
	.word	0,0
	.word	0X3FE6719F,0X3601671A
	.word	0X3FD6B3BB,0X2235A000
	.word	0XBD4784ED,0X42B666CC
	.word	0,0
	.word	0X3FE669C3,0X1075AB40
	.word	0X3FD6CA28,0XD2E34000
	.word	0X3D430AD0,0XB8C49166
	.word	0,0
	.word	0X3FE661EC,0X6A5122F9
	.word	0X3FD6E08E,0XAA2BA000
	.word	0X3D1E38C1,0X39318D71
	.word	0,0
	.word	0X3FE65A1B,0X3DD13357
	.word	0X3FD6F6EC,0XAD8B2000
	.word	0X3D249058,0XFDF08376
	.word	0,0
	.word	0X3FE6524F,0X853B4AA3
	.word	0X3FD70D42,0XE278A000
	.word	0XBD4B9454,0XB320475E
	.word	0,0
	.word	0X3FE64A89,0X3ADCD25F
	.word	0X3FD72391,0X4E650000
	.word	0X3D0C1D52,0XBDC87D8A
	.word	0,0
	.word	0X3FE642C8,0X590B2164
	.word	0X3FD739D7,0XF6BBE000
	.word	0XBD4FF2C6,0X3B67580A
	.word	0,0
	.word	0X3FE63B0C,0XDA236E1C
	.word	0X3FD75016,0XE0E2C000
	.word	0XBD3677E8,0XB799D03C
	.word	0,0
	.word	0X3FE63356,0XB88AC0DE
	.word	0X3FD7664E,0X1239E000
	.word	0XBD30C4FB,0X6AEB27AF
	.word	0,0
	.word	0X3FE62BA5,0XEEADE65E
	.word	0X3FD77C7D,0X901BC000
	.word	0XBD3BAFC1,0X943804E0
	.word	0,0
	.word	0X3FE623FA,0X77016240
	.word	0X3FD792A5,0X5FDD4000
	.word	0X3D3E89F0,0X57691FEA
	.word	0,0
	.word	0X3FE61C54,0X4C0161C5
	.word	0X3FD7A8C5,0X86CE0000
	.word	0XBD4577C1,0X4AED80C0
	.word	0,0
	.word	0X3FE614B3,0X6831AE94
	.word	0X3FD7BEDE,0XA37A000
	.word	0X3D4F7F3C,0X3E1A33FF
	.word	0,0
	.word	0X3FE60D17,0XC61DA198
	.word	0X3FD7D4EE,0XEF5EE000
	.word	0X3D48DBF9,0XC1C3609C
	.word	0,0
	.word	0X3FE60581,0X60581606
	.word	0X3FD7EAF8,0X3B82A000
	.word	0X3D4F86C9,0X674BCF69
	.word	0,0
	.word	0X3FE5FDF0,0X317B5C6F
	.word	0X3FD800F9,0XF3DCA000
	.word	0XBD466A18,0XC911BAF
	.word	0,0
	.word	0X3FE5F664,0X34292DFC
	.word	0X3FD816F4,0X1DA0E000
	.word	0XBD46D494,0X91E5025C
	.word	0,0
	.word	0X3FE5EEDD,0X630A9FB3
	.word	0X3FD82CE6,0XBDFE4000
	.word	0X3D4B39B3,0XDD334022
	.word	0,0
	.word	0X3FE5E75B,0XB8D015E7
	.word	0X3FD842D1,0XDA1E8000
	.word	0X3D462E92,0X7628CBC2
	.word	0,0
	.word	0X3FE5DFDF,0X303137B6
	.word	0X3FD858B5,0X7725C000
	.word	0X3D48855B,0X2ACDA048
	.word	0,0
	.word	0X3FE5D867,0XC3ECE2A5
	.word	0X3FD86E91,0X9A330000
	.word	0X3D474013,0XF9B16FEB
	.word	0,0
	.word	0X3FE5D0F5,0X6EC91E57
	.word	0X3FD88466,0X48600000
	.word	0X3D388FE9,0XDF324F6A
	.word	0,0
	.word	0X3FE5C988,0X2B931057
	.word	0X3FD89A33,0X86C14000
	.word	0X3D22D5AD,0X38C40882
	.word	0,0
	.word	0X3FE5C21F,0XF51EF005
	.word	0X3FD8AFF9,0X5A662000
	.word	0XBD41002E,0X2C0D76A1
	.word	0,0
	.word	0X3FE5BABC,0XC647FA91
	.word	0X3FD8C5B7,0XC858C000
	.word	0XBD46EAF0,0X55A7EFD0
	.word	0,0
	.word	0X3FE5B35E,0X99F06714
	.word	0X3FD8DB6E,0XD59E2000
	.word	0X3D3CB29E,0X13B3EBAE
	.word	0,0
	.word	0X3FE5AC05,0X6B015AC0
	.word	0X3FD8F11E,0X87366000
	.word	0X3D263BF0,0XBB4EAB4C
	.word	0,0
	.word	0X3FE5A4B1,0X346ADD2B
	.word	0X3FD906C6,0XE21C4000
	.word	0X3D3D4E79,0X8E962C0C
	.word	0,0
	.word	0X3FE59D61,0XF123CCAA
	.word	0X3FD91C67,0XEB45A000
	.word	0X3D407B0F,0X8FA8EE5B
	.word	0,0
	.word	0X3FE59617,0X9C29D2CE
	.word	0X3FD93201,0XA7A36000
	.word	0XBD4992BC,0X9C47D591
	.word	0,0
	.word	0X3FE58ED2,0X308158ED
	.word	0X3FD94794,0X1C212000
	.word	0XBD420A8B,0X6645D706
	.word	0,0
	.word	0X3FE58791,0XA9357CCE
	.word	0X3FD95D1F,0X4DA5C000
	.word	0X3D44143F,0XF0110E3
	.word	0,0
	.word	0X3FE58056,0X1580560
	.word	0X3FD972A3,0X41136000
	.word	0XBD4D4F2D,0X1FB16DA4
	.word	0,0
	.word	0X3FE5791F,0X34015792
	.word	0X3FD9881F,0XFB46A000
	.word	0X3D3BC032,0XB84E5463
	.word	0,0
	.word	0X3FE571ED,0X3C506B3A
	.word	0X3FD99D95,0X8117E000
	.word	0X3D015975,0X25DD88F0
	.word	0,0
	.word	0X3FE56AC0,0X156AC015
	.word	0X3FD9B303,0XD75A4000
	.word	0XBD43C1C6,0X85B184C0
	.word	0,0
	.word	0X3FE56397,0XBA7C52E2
	.word	0X3FD9C86B,0X2DC0000
	.word	0X3D40C537,0X408A4B11
	.word	0,0
	.word	0X3FE55C74,0X26B79286
	.word	0X3FD9DDCB,0X866E000
	.word	0X3D3D066F,0XCA3E126B
	.word	0,0
	.word	0X3FE55555,0X55555555
	.word	0XBFD26962,0X1134E000
	.word	0X3D31B61F,0X10522625
	.word	0,0
	.word	0X3FE54E3B,0X4194CE66
	.word	0XBFD25410,0X494E5000
	.word	0XBD3B1D7A,0XC0EF77F2
	.word	0,0
	.word	0X3FE54725,0XE6BB82FE
	.word	0XBFD23EC5,0X991EC000
	.word	0X3D36DBE4,0X48A2E522
	.word	0,0
	.word	0X3FE54015,0X40154015
	.word	0XBFD22981,0XFBEF8000
	.word	0X3D3A1421,0X609580DA
	.word	0,0
	.word	0X3FE53909,0X48F40FEB
	.word	0XBFD21445,0X6D0EC000
	.word	0X3D3CAF04,0X28B728A3
	.word	0,0
	.word	0X3FE53201,0XFCB02FB1
	.word	0XBFD1FF0F,0XE7CF4000
	.word	0XBD3E9D5B,0X513FF0C1
	.word	0,0
	.word	0X3FE52AFF,0X56A8054B
	.word	0XBFD1E9E1,0X6788A000
	.word	0X3D382EAE,0XD3C8B65E
	.word	0,0
	.word	0X3FE52401,0X52401524
	.word	0XBFD1D4B9,0XE796C000
	.word	0XBD222A66,0X7C42E56D
	.word	0,0
	.word	0X3FE51D07,0XEAE2F815
	.word	0XBFD1BF99,0X635A7000
	.word	0X3D31AC89,0X575C2125
	.word	0,0
	.word	0X3FE51613,0X1C015161
	.word	0XBFD1AA7F,0XD638D000
	.word	0XBD29F60A,0X9616F7A0
	.word	0,0
	.word	0X3FE50F22,0XE111C4C5
	.word	0XBFD1956D,0X3B9BC000
	.word	0XBD27D2F7,0X3AD1AA14
	.word	0,0
	.word	0X3FE50837,0X3590EC9C
	.word	0XBFD18061,0X8EF19000
	.word	0X3D3482FF,0XC86D38E5
	.word	0,0
	.word	0X3FE50150,0X15015015
	.word	0XBFD16B5C,0XCBAD0000
	.word	0X3D323299,0X42D74BF
	.word	0,0
	.word	0X3FE4FA6D,0X7AEB597C
	.word	0XBFD1565E,0XED456000
	.word	0X3CEE75AD,0XFB6ABA25
	.word	0,0
	.word	0X3FE4F38F,0X62DD4C9B
	.word	0XBFD14167,0XEF367000
	.word	0XBD3E0C07,0X824DAAF5
	.word	0,0
	.word	0X3FE4ECB5,0XC86B3D24
	.word	0XBFD12C77,0XCD007000
	.word	0XBD13B294,0X8A11F797
	.word	0,0
	.word	0X3FE4E5E0,0XA72F0539
	.word	0XBFD1178E,0X8227E000
	.word	0XBD31EF78,0XCE2D07F2
	.word	0,0
	.word	0X3FE4DF0F,0XFAC83C01
	.word	0XBFD102AC,0XA35D000
	.word	0X3D2F1FBD,0XDFDFD686
	.word	0,0
	.word	0X3FE4D843,0XBEDC2C4C
	.word	0XBFD0EDD0,0X60B78000
	.word	0XBD0019B5,0X2D8435F5
	.word	0,0
	.word	0X3FE4D17B,0XEF15CB4E
	.word	0XBFD0D8FB,0X813EB000
	.word	0XBD1EE8C8,0X8753FA35
	.word	0,0
	.word	0X3FE4CAB8,0X8725AF6E
	.word	0XBFD0C42D,0X67616000
	.word	0XBD27188B,0X163CEAE9
	.word	0,0
	.word	0X3FE4C3F9,0X82C20723
	.word	0XBFD0AF66,0XEB9E000
	.word	0XBD23C7C3,0XF528D80A
	.word	0,0
	.word	0X3FE4BD3E,0XDDA68FE1
	.word	0XBFD09AA5,0X72E6C000
	.word	0XBD3B50A1,0XE1734342
	.word	0,0
	.word	0X3FE4B688,0X93948D1C
	.word	0XBFD085EB,0X8F8AE000
	.word	0XBD3E5D51,0X3F45FE7B
	.word	0,0
	.word	0X3FE4AFD6,0XA052BF5B
	.word	0XBFD07138,0X604D6000
	.word	0X3D3E7632,0X4E912B17
	.word	0,0
	.word	0X3FE4A928,0XFFAD5B5C
	.word	0XBFD05C8B,0XE0D96000
	.word	0XBD2AD0F1,0XC77CCB58
	.word	0,0
	.word	0X3FE4A27F,0XAD76014A
	.word	0XBFD047E6,0XCDE8000
	.word	0XBD2DBDF1,0XD397F3C
	.word	0,0
	.word	0X3FE49BDA,0XA583B401
	.word	0XBFD03346,0XE0106000
	.word	0XBCF89FF8,0XA966395C
	.word	0,0
	.word	0X3FE49539,0XE3B2D067
	.word	0XBFD01EAE,0X5626C000
	.word	0XBD3A43DC,0XFADE85AE
	.word	0,0
	.word	0X3FE48E9D,0X63E504D1
	.word	0XBFD00A1C,0X6ADDA000
	.word	0XBD31CD8D,0X688B9E18
	.word	0,0
	.word	0X3FE48805,0X22014880
	.word	0XBFCFEB22,0X33EA0000
	.word	0XBD2F3418,0XDE00938B
	.word	0,0
	.word	0X3FE48171,0X19F3D325
	.word	0XBFCFC218,0XBE620000
	.word	0XBD34BBA4,0X6F1CF6A0
	.word	0,0
	.word	0X3FE47AE1,0X47AE147B
	.word	0XBFCF991C,0X6CB3C000
	.word	0X3D390D04,0XCD7CC834
	.word	0,0
	.word	0X3FE47455,0XA726ABF2
	.word	0XBFCF702D,0X36778000
	.word	0X3D108195,0X16673E23
	.word	0,0
	.word	0X3FE46DCE,0X34596066
	.word	0XBFCF474B,0X134E0000
	.word	0X3D3BAE49,0XF1DF7B5E
	.word	0,0
	.word	0X3FE4674A,0XEB4717E9
	.word	0XBFCF1E75,0XFADFA000
	.word	0X3D20862B,0X25D83F6D
	.word	0,0
	.word	0X3FE460CB,0XC7F5CF9A
	.word	0XBFCEF5AD,0XE4DD0000
	.word	0X3CCA2115,0X65BB8E11
	.word	0,0
	.word	0X3FE45A50,0XC670938F
	.word	0XBFCECCF2,0XC8FEA000
	.word	0X3D3BEC63,0XA3E75640
	.word	0,0
	.word	0X3FE453D9,0XE2C776CA
	.word	0XBFCEA444,0X9F04A000
	.word	0XBD35E916,0X63732A36
	.word	0,0
	.word	0X3FE44D67,0X190F8B43
	.word	0XBFCE7BA3,0X5EB78000
	.word	0X3D0D5EEE,0X23793649
	.word	0,0
	.word	0X3FE446F8,0X6562D9FB
	.word	0XBFCE530E,0XFFE72000
	.word	0X3D3FDBDB,0XB13F7C18
	.word	0,0
	.word	0X3FE4408D,0XC3E05B22
	.word	0XBFCE2A87,0X7A6B2000
	.word	0XBD382381,0X7787081A
	.word	0,0
	.word	0X3FE43A27,0X30ABEE4D
	.word	0XBFCE020C,0XC6236000
	.word	0X3D252B00,0XADB91424
	.word	0,0
	.word	0X3FE433C4,0XA7EE52B4
	.word	0XBFCDD99E,0XDAF6E000
	.word	0X3D302EC6,0X69C756EB
	.word	0,0
	.word	0X3FE42D66,0X25D51F87
	.word	0XBFCDB13D,0XB0D48000
	.word	0XBD32806A,0X847527E6
	.word	0,0
	.word	0X3FE4270B,0XA692BC4D
	.word	0XBFCD88E9,0X3FB30000
	.word	0X3D375F28,0X234BF51
	.word	0,0
	.word	0X3FE420B5,0X265E5951
	.word	0XBFCD60A1,0X7F904000
	.word	0X3D35D6E0,0X6FC20D39
	.word	0,0
	.word	0X3FE41A62,0XA173E821
	.word	0XBFCD3866,0X68720000
	.word	0X3D373650,0XB38932BC
	.word	0,0
	.word	0X3FE41414,0X14141414
	.word	0XBFCD1037,0XF2656000
	.word	0X3D084A7E,0X75B6F6E4
	.word	0,0
	.word	0X3FE40DC9,0X7A843AE8
	.word	0XBFCCE816,0X157F2000
	.word	0X3D29E0AB,0XA2099515
	.word	0,0
	.word	0X3FE40782,0XD10E6566
	.word	0XBFCCC000,0XC9DB4000
	.word	0X3D1D6D58,0X5D57AFF9
	.word	0,0
	.word	0X3FE40140,0X14014014
	.word	0XBFCC97F8,0X79D4000
	.word	0XBD23B161,0XA8C6E6C5
	.word	0,0
	.word	0X3FE3FB01,0X3FB013FB
	.word	0XBFCC6FFB,0XC6F00000
	.word	0XBD3EE138,0XD3A69D43
	.word	0,0
	.word	0X3FE3F4C6,0X5072BF74
	.word	0XBFCC480C,0X5C000
	.word	0XBD39A294,0XD5E44E76
	.word	0,0
	.word	0X3FE3EE8F,0X42A5AF07
	.word	0XBFCC2028,0XAB180000
	.word	0X3D292E0E,0XE55C7AC6
	.word	0,0
	.word	0X3FE3E85C,0X12A9D651
	.word	0XBFCBF851,0XC0676000
	.word	0X3D35420E,0X4C0854AD
	.word	0,0
	.word	0X3FE3E22C,0XBCE4A902
	.word	0XBFCBD087,0X383BE000
	.word	0X3D2D4BC4,0X595412B6
	.word	0,0
	.word	0X3FE3DC01,0X3DC013DC
	.word	0XBFCBA8C9,0XAE4A000
	.word	0XBD3A32E7,0XF44432DA
	.word	0,0
	.word	0X3FE3D5D9,0X91AA75C6
	.word	0XBFCB8117,0X30B82000
	.word	0XBD1E9068,0X3B9CD768
	.word	0,0
	.word	0X3FE3CFB5,0XB51698EB
	.word	0XBFCB5971,0XA213A000
	.word	0XBD39B50E,0X83AA91DF
	.word	0,0
	.word	0X3FE3C995,0XA47BABE7
	.word	0XBFCB31D8,0X575BC000
	.word	0XBD3C794E,0X562A63CB
	.word	0,0
	.word	0X3FE3C379,0X5C553AFB
	.word	0XBFCB0A4B,0X48FC2000
	.word	0X3D22E72D,0X5C3998ED
	.word	0,0
	.word	0X3FE3BD60,0XD9232955
	.word	0XBFCAE2CA,0X6F672000
	.word	0XBD37A8D5,0XAE54F550
	.word	0,0
	.word	0X3FE3B74C,0X1769AA5C
	.word	0XBFCABB55,0XC316A000
	.word	0X3D38A65A,0XCAF14CD8
	.word	0,0
	.word	0X3FE3B13B,0X13B13B14
	.word	0XBFCA93ED,0X3C8AE000
	.word	0X3D287243,0X50562169
	.word	0,0
	.word	0X3FE3AB2D,0XCA869B81
	.word	0XBFCA6C90,0XD44B8000
	.word	0X3D3F63B7,0XF037B0C6
	.word	0,0
	.word	0X3FE3A524,0X387AC822
	.word	0XBFCA4540,0X82E6A000
	.word	0XBD360A77,0XC81F7171
	.word	0,0
	.word	0X3FE39F1E,0X5A22F36E
	.word	0XBFCA1DFC,0X40F1C000
	.word	0X3D301E0F,0X4F3781
	.word	0,0
	.word	0X3FE3991C,0X2C187F63
	.word	0XBFC9F6C4,0X708A000
	.word	0X3D3337D9,0X4BCD3F43
	.word	0,0
	.word	0X3FE3931D,0XAAF8F721
	.word	0XBFC9CF97,0XCDCE0000
	.word	0XBD3D862F,0X10C414E3
	.word	0,0
	.word	0X3FE38D22,0XD366088E
	.word	0XBFC9A877,0X8DEBA000
	.word	0XBD3470FA,0X3EFEC390
	.word	0,0
	.word	0X3FE3872B,0XA2057E04
	.word	0XBFC98163,0X4011A000
	.word	0XBD34EADD,0X9E9045E2
	.word	0,0
	.word	0X3FE38138,0X13813814
	.word	0XBFC95A5A,0XDCF70000
	.word	0XBD07F228,0X58A0FF6F
	.word	0,0
	.word	0X3FE37B48,0X24872744
	.word	0XBFC9335E,0X5D594000
	.word	0XBD33115C,0X3ABD47DA
	.word	0,0
	.word	0X3FE3755B,0XD1C945EE
	.word	0XBFC90C6D,0XB9FCC000
	.word	0X3D1935F5,0X7718D7CA
	.word	0,0
	.word	0X3FE36F73,0X17FD9212
	.word	0XBFC8E588,0XEBAC2000
	.word	0XBD3B7D5C,0XAB2D1140
	.word	0,0
	.word	0X3FE3698D,0XF3DE0748
	.word	0XBFC8BEAF,0XEB390000
	.word	0X3D073D54,0XAAE92CD1
	.word	0,0
	.word	0X3FE363AC,0X622898B1
	.word	0XBFC897E2,0XB17B2000
	.word	0X3D296B37,0X380CBE9E
	.word	0,0
	.word	0X3FE35DCE,0X5F9F2AF8
	.word	0XBFC87121,0X3750E000
	.word	0XBD3328EB,0X42F9AF75
	.word	0,0
	.word	0X3FE357F3,0XE9078E5B
	.word	0XBFC84A6B,0X759F6000
	.word	0X3D3DA280,0X2ADF8609
	.word	0,0
	.word	0X3FE3521C,0XFB2B78C1
	.word	0XBFC823C1,0X6551A000
	.word	0XBD1E0DDB,0X9A631E83
	.word	0,0
	.word	0X3FE34C49,0X92D87FD9
	.word	0XBFC7FD22,0XFF59A000
	.word	0X3D158BEB,0XF457B7D2
	.word	0,0
	.word	0X3FE34679,0XACE01346
	.word	0XBFC7D690,0X3CAF6000
	.word	0X3D24C06B,0X17C301D7
	.word	0,0
	.word	0X3FE340AD,0X461776D3
	.word	0XBFC7B009,0X16516000
	.word	0X3D3AE75F,0XCB067E57
	.word	0,0
	.word	0X3FE33AE4,0X5B57BCB2
	.word	0XBFC7898D,0X85444000
	.word	0XBD38E67B,0XE3DBAF3F
	.word	0,0
	.word	0X3FE3351E,0XE97DBFC6
	.word	0XBFC7631D,0X82936000
	.word	0X3D25E77D,0XC7C5F3E1
	.word	0,0
	.word	0X3FE32F5C,0XED6A1DFA
	.word	0XBFC73CB9,0X74FE000
	.word	0X3D3D66A9,0XD0005A6
	.word	0,0
	.word	0X3FE3299E,0X6401329A
	.word	0XBFC71660,0XC914000
	.word	0XBCE51B15,0X7CEC3838
	.word	0,0
	.word	0X3FE323E3,0X4A2B10BF
	.word	0XBFC6F012,0X8B756000
	.word	0XBD357739,0XD31EF0F
	.word	0,0
	.word	0X3FE31E2B,0X9CD37DC2
	.word	0XBFC6C9D0,0X7D204000
	.word	0X3CDC73FA,0XFD9B2DCA
	.word	0,0
	.word	0X3FE31877,0X58E9EBB6
	.word	0XBFC6A399,0XDABBE000
	.word	0X3D38F934,0XE66A15A6
	.word	0,0
	.word	0X3FE312C6,0X7B6173EE
	.word	0XBFC67D6E,0X9D786000
	.word	0X3D311E88,0X30A706D3
	.word	0,0
	.word	0X3FE30D19,0X130D190
	.word	0XBFC6574E,0XBE8C2000
	.word	0X3D398C1D,0X34F0F462
	.word	0,0
	.word	0X3FE3076E,0XE7525C2C
	.word	0XBFC6313A,0X37336000
	.word	0X3D144DF5,0X4F21EA6D
	.word	0,0
	.word	0X3FE301C8,0X2AC40260
	.word	0XBFC60B31,0XB0A000
	.word	0X3D371456,0XC988F814
	.word	0,0
	.word	0X3FE2FC24,0XC8874486
	.word	0XBFC5E533,0X144C2000
	.word	0X3D31CE0B,0XF3B290EA
	.word	0,0
	.word	0X3FE2F684,0XBDA12F68
	.word	0XBFC5BF40,0X6B544000
	.word	0X3D127023,0XEB68981C
	.word	0,0
	.word	0X3FE2F0E8,0X71A5703
	.word	0XBFC59958,0XFF1D6000
	.word	0X3D3A1D05,0X9769CA05
	.word	0,0
	.word	0X3FE2EB4E,0XA1FED14B
	.word	0XBFC5737C,0XC9018000
	.word	0XBD39BAA7,0XA6B887F6
	.word	0,0
	.word	0X3FE2E5B8,0X8B5E3104
	.word	0XBFC54DAB,0XC2610000
	.word	0XBD2746FE,0XE5C8D0D8
	.word	0,0
	.word	0X3FE2E025,0XC04B8097
	.word	0XBFC527E5,0XE4A1C000
	.word	0X3D34E60B,0X8D4B411D
	.word	0,0
	.word	0X3FE2DA96,0X3DDD3CFB
	.word	0XBFC5022B,0X292F6000
	.word	0XBD348A05,0XFF36A25B
	.word	0,0
	.word	0X3FE2D50A,0X12D50A0
	.word	0XBFC4DC7B,0X897BC000
	.word	0XBD0C79B6,0XAE1FF0F
	.word	0,0
	.word	0X3FE2CF81,0X7590E67
	.word	0XBFC4B6D6,0XFEFE2000
	.word	0XBD1522EC,0XF56E7952
	.word	0,0
	.word	0X3FE2C9FB,0X4D812CA0
	.word	0XBFC4913D,0X8333C000
	.word	0X3D353E43,0X558124C4
	.word	0,0
	.word	0X3FE2C478,0XD0C9C013
	.word	0XBFC46BAF,0XF9F6000
	.word	0X3D1249CD,0X790841A
	.word	0,0
	.word	0X3FE2BEF9,0X8E5A3711
	.word	0XBFC4462B,0X9DC9C000
	.word	0X3D384858,0XA711B062
	.word	0,0
	.word	0X3FE2B97D,0X835D548E
	.word	0XBFC420B3,0X27410000
	.word	0X3D116282,0XC85A0884
	.word	0,0
	.word	0X3FE2B404,0XAD012B40
	.word	0XBFC3FB45,0XA5992000
	.word	0XBD319713,0XC0CAE559
	.word	0,0
	.word	0X3FE2AE8F,0X87718D0
	.word	0XBFC3D5E3,0X126BC000
	.word	0XBD13FB2F,0X85096C4B
	.word	0,0
	.word	0X3FE2A91C,0X92F3C105
	.word	0XBFC3B08B,0X67580000
	.word	0X3D3AADE8,0XF29320FB
	.word	0,0
	.word	0X3FE2A3AD,0X49AF0907
	.word	0XBFC38B3E,0X9E028000
	.word	0X3D370EF0,0X545C17F9
	.word	0,0
	.word	0X3FE29E41,0X29E4129E
	.word	0XBFC365FC,0XB015A000
	.word	0X3D3FD3A0,0XAFB9691B
	.word	0,0
	.word	0X3FE298D8,0X30D13780
	.word	0XBFC340C5,0X97412000
	.word	0X3D37A3DC,0XF7D9D386
	.word	0,0
	.word	0X3FE29372,0X5BB804A5
	.word	0XBFC31B99,0X4D3A4000
	.word	0XBD3F098E,0XE3A50810
	.word	0,0
	.word	0X3FE28E0F,0XA7DD35A3
	.word	0XBFC2F677,0XCBBC0000
	.word	0XBD352B30,0X2160F40D
	.word	0,0
	.word	0X3FE288B0,0X1288B013
	.word	0XBFC2D161,0XC868000
	.word	0XBD039D6C,0XCB81B4A1
	.word	0,0
	.word	0X3FE28353,0X99057EFD
	.word	0XBFC2AC55,0X95F6000
	.word	0X3D1D3466,0XD0C6C8A8
	.word	0,0
	.word	0X3FE27DFA,0X38A1CE4D
	.word	0XBFC28753,0XBC11A000
	.word	0XBD37494E,0X359302E6
	.word	0,0
	.word	0X3FE278A3,0XEEAEE650
	.word	0XBFC2625D,0X1E6DE000
	.word	0X3CF52962,0XF09E3D82
	.word	0,0
	.word	0X3FE27350,0XB8812735
	.word	0XBFC23D71,0X2A49C000
	.word	0XBD100D23,0X8FD3DF5C
	.word	0,0
	.word	0X3FE26E00,0X9370049C
	.word	0XBFC2188F,0XD9808000
	.word	0X3D3B3A1E,0X7F50C701
	.word	0,0
	.word	0X3FE268B3,0X7CD60127
	.word	0XBFC1F3B9,0X25F26000
	.word	0X3D15F74E,0X9B083633
	.word	0,0
	.word	0X3FE26369,0X7210AA18
	.word	0XBFC1CEED,0X9854000
	.word	0X3D315C1C,0X39192AF9
	.word	0,0
	.word	0X3FE25E22,0X708092F1
	.word	0XBFC1AA2B,0X7E240000
	.word	0X3D31AC38,0XDDE3B366
	.word	0,0
	.word	0X3FE258DE,0X75895121
	.word	0XBFC18574,0X7DBEC000
	.word	0XBD3E6744,0X45BD9B49
	.word	0,0
	.word	0X3FE2539D,0X7E9177B2
	.word	0XBFC160C8,0X24B2000
	.word	0XBD2EC2D2,0XA9009E3D
	.word	0,0
	.word	0X3FE24E5F,0X89029305
	.word	0XBFC13C26,0X5C3A000
	.word	0X3D2CF5FD,0XD94F6509
	.word	0,0
	.word	0X3FE24924,0X92492492
	.word	0XBFC1178E,0X8227E000
	.word	0XBD21EF78,0XCE2D07F2
	.word	0,0
	.word	0X3FE243EC,0X97D49EAE
	.word	0XBFC0F301,0X717D0000
	.word	0X3D3E09B4,0X41AE86C5
	.word	0,0
	.word	0X3FE23EB7,0X9717605B
	.word	0XBFC0CE7E,0XCDCCC000
	.word	0XBD14652D,0XABFF5447
	.word	0,0
	.word	0X3FE23985,0X8D86B11F
	.word	0XBFC0AA06,0X91268000
	.word	0X3D345519,0XD7032129
	.word	0,0
	.word	0X3FE23456,0X789ABCDF
	.word	0XBFC08598,0XB59E4000
	.word	0X3D27E5DD,0X7009902C
	.word	0,0
	.word	0X3FE22F2A,0X55CE8FC5
	.word	0XBFC06135,0X354D4000
	.word	0XBD363046,0X28340EE9
	.word	0,0
	.word	0X3FE22A01,0X22A0122A
	.word	0XBFC03CDC,0XA51E000
	.word	0XBD381A9C,0XF169FC5C
	.word	0,0
	.word	0X3FE224DA,0XDC900489
	.word	0XBFC0188D,0X2ECF6000
	.word	0XBD03F965,0X1CFF9DFE
	.word	0,0
	.word	0X3FE21FB7,0X8121FB78
	.word	0XBFBFE891,0X39DBC000
	.word	0XBD356594,0XD82F7A82
	.word	0,0
	.word	0X3FE21A97,0XDDC5BA7
	.word	0XBFBFA01C,0X9DB58000
	.word	0X3D08F351,0XFA48A730
	.word	0,0
	.word	0X3FE21579,0X804855E6
	.word	0XBFBF57BC,0X7D900000
	.word	0XBD176A6C,0X9EA8B04E
	.word	0,0
	.word	0X3FE2105E,0XD5F1E336
	.word	0XBFBF0F70,0XCDD98000
	.word	0XBD32E31F,0X6C272C1E
	.word	0,0
	.word	0X3FE20B47,0XC67C0D9
	.word	0XBFBEC739,0X830A0000
	.word	0XBD311FCB,0XA80CDD10
	.word	0,0
	.word	0X3FE20632,0X213B6C6D
	.word	0XBFBE7F16,0X91A34000
	.word	0X3D32C1C5,0X9BC77BFA
	.word	0,0
	.word	0X3FE20120,0X12012012
	.word	0XBFBE3707,0XEE304000
	.word	0XBD20F684,0XE6766ABD
	.word	0,0
	.word	0X3FE1FC10,0XDC4FCE8B
	.word	0XBFBDEF0D,0X8D468000
	.word	0X3D324750,0X412E9A74
	.word	0,0
	.word	0X3FE1F704,0X7DC11F70
	.word	0XBFBDA727,0X63844000
	.word	0XBD1A8940,0X1FA71733
	.word	0,0
	.word	0X3FE1F1FA,0XF3F16B64
	.word	0XBFBD5F55,0X65920000
	.word	0XBD30E239,0XCC185469
	.word	0,0
	.word	0X3FE1ECF4,0X3C7FB84C
	.word	0XBFBD1797,0X88218000
	.word	0XBD336433,0XB5EFBEED
	.word	0,0
	.word	0X3FE1E7F0,0X550DB594
	.word	0XBFBCCFED,0XBFEE0000
	.word	0XBD33A823,0X2FE71256
	.word	0,0
	.word	0X3FE1E2EF,0X3B3FB874
	.word	0XBFBC8858,0X1BC4000
	.word	0XBD2646D1,0XC65AACD3
	.word	0,0
	.word	0X3FE1DDF0,0XECBCB841
	.word	0XBFBC40D6,0X425A4000
	.word	0XBD3CB112,0X1D1930DD
	.word	0,0
	.word	0X3FE1D8F5,0X672E4ABD
	.word	0XBFBBF968,0X769FC000
	.word	0XBD24218C,0X8D824283
	.word	0,0
	.word	0X3FE1D3FC,0XA840A074
	.word	0XBFBBB20E,0X936D8000
	.word	0X3D368BA8,0X35459B8E
	.word	0,0
	.word	0X3FE1CF06,0XADA2811D
	.word	0XBFBB6AC8,0X8DAD4000
	.word	0XBD3B1BDF,0XF50225C7
	.word	0,0
	.word	0X3FE1CA13,0X750547FE
	.word	0XBFBB2396,0X5A530000
	.word	0X3CEFF64E,0XEA137079
	.word	0,0
	.word	0X3FE1C522,0XFC1CE059
	.word	0XBFBADC77,0XEE5B0000
	.word	0X3D3573B2,0X9C31904
	.word	0,0
	.word	0X3FE1C035,0X409FC1DF
	.word	0XBFBA956D,0X3ECAC000
	.word	0XBD3E6379,0X4C02C4AF
	.word	0,0
	.word	0X3FE1BB4A,0X4046ED29
	.word	0XBFBA4E76,0X40B1C000
	.word	0X3D0E42B6,0XB94407C8
	.word	0,0
	.word	0X3FE1B661,0XF8CDE833
	.word	0XBFBA0792,0XE9278000
	.word	0X3D0A9CE6,0XC9AD51BF
	.word	0,0
	.word	0X3FE1B17C,0X67F2BAE3
	.word	0XBFB9C0C3,0X2D4D4000
	.word	0X3D3AB7C0,0X9E838668
	.word	0,0
	.word	0X3FE1AC99,0X8B75EB90
	.word	0XBFB97A07,0X24CC000
	.word	0X3CF8BCC1,0X732093CE
	.word	0,0
	.word	0X3FE1A7B9,0X611A7B96
	.word	0XBFB9335E,0X5D594000
	.word	0XBD23115C,0X3ABD47DA
	.word	0,0
	.word	0X3FE1A2DB,0XE6A5E3E4
	.word	0XBFB8ECC9,0X33AEC000
	.word	0X3D222F39,0XBE67F7AA
	.word	0,0
	.word	0X3FE19E01,0X19E0119E
	.word	0XBFB8A647,0X7A91C000
	.word	0XBD3C28C0,0XAF9BD6DF
	.word	0,0
	.word	0X3FE19928,0XF89362B7
	.word	0XBFB85FD9,0X27508000
	.word	0X3D35B818,0X19970C1C
	.word	0,0
	.word	0X3FE19453,0X808CA29C
	.word	0XBFB8197E,0X2F410000
	.word	0X3D3C0FE4,0X60D20041
	.word	0,0
	.word	0X3FE18F80,0XAF9B06DC
	.word	0XBFB7D336,0X87C28000
	.word	0XBD33C88C,0X3E706706
	.word	0,0
	.word	0X3FE18AB0,0X83902BDB
	.word	0XBFB78D02,0X263D8000
	.word	0XBD069B57,0X94B69FB7
	.word	0,0
	.word	0X3FE185E2,0XFA401186
	.word	0XBFB746E1,0X228000
	.word	0X3D3126D1,0X6E1E21D2
	.word	0,0
	.word	0X3FE18118,0X11811812
	.word	0XBFB700D3,0XAEAC000
	.word	0XBCEC1E8D,0XA99DED32
	.word	0,0
	.word	0X3FE17C4F,0XC72BFCB9
	.word	0XBFB6BAD8,0X3C188000
	.word	0XBD0DAF3C,0XC08926AE
	.word	0,0
	.word	0X3FE1778A,0X191BD684
	.word	0XBFB674F0,0X89364000
	.word	0XBD3A7999,0X4C9D3302
	.word	0,0
	.word	0X3FE172C7,0X52E1316
	.word	0XBFB62F1B,0XE7D78000
	.word	0X3D217995,0X7ED63C4E
	.word	0,0
	.word	0X3FE16E06,0X89427379
	.word	0XBFB5E95A,0X4D978000
	.word	0XBD31CB7C,0XE1D17171
	.word	0,0
	.word	0X3FE16948,0XA33B08FA
	.word	0XBFB5A3AB,0XB01AC000
	.word	0XBD3E2574,0X9E6AFA18
	.word	0,0
	.word	0X3FE1648D,0X50FC3201
	.word	0XBFB55E10,0X50E0000
	.word	0XBD0C1D74,0XC53C72E
	.word	0,0
	.word	0X3FE15FD4,0X906C96F1
	.word	0XBFB51887,0X42260000
	.word	0XBD330A1D,0X96258B3E
	.word	0,0
	.word	0X3FE15B1E,0X5F75270D
	.word	0XBFB4D311,0X5D208000
	.word	0X3CF53A25,0X82F4E1EF
	.word	0,0
	.word	0X3FE1566A,0XBC011567
	.word	0XBFB48DAE,0X4BC30000
	.word	0XBD30185B,0X208C200C
	.word	0,0
	.word	0X3FE151B9,0XA3FDD5C9
	.word	0XBFB4485E,0X3DBC000
	.word	0XBD3FAD46,0XE8D26AB7
	.word	0,0
	.word	0X3FE14D0B,0X155B19AE
	.word	0XBFB40320,0X7B414000
	.word	0XBD26FD84,0XAA8157C0
	.word	0,0
	.word	0X3FE1485F,0XE0ACD3B
	.word	0XBFB3BDF5,0XA7D20000
	.word	0X3D319BD0,0XAD125895
	.word	0,0
	.word	0X3FE143B5,0X8C01143B
	.word	0XBFB378DD,0X7F748000
	.word	0XBD371411,0X28F1FACA
	.word	0,0
	.word	0X3FE13F0E,0X8D344724
	.word	0XBFB333D7,0XF8184000
	.word	0X3CE692B6,0XA81B8848
	.word	0,0
	.word	0X3FE13A6A,0XF9CF01E
	.word	0XBFB2EEE5,0X7B40000
	.word	0XBD08081E,0XDD77C860
	.word	0,0
	.word	0X3FE135C8,0X1135C811
	.word	0XBFB2AA04,0XA4470000
	.word	0XBD37A48B,0XA8B1CB41
	.word	0,0
	.word	0X3FE13128,0X8FFBB3B6
	.word	0XBFB26536,0XC3D8C000
	.word	0XBD0B4BAC,0X97C5BA3
	.word	0,0
	.word	0X3FE12C8B,0X89EDC0AC
	.word	0XBFB2207B,0X5C784000
	.word	0XBD349D8C,0XFC10C7BF
	.word	0,0
	.word	0X3FE127F0,0XFD0D2295
	.word	0XBFB1DBD2,0X643D0000
	.word	0XBD390B24,0XD977C494
	.word	0,0
	.word	0X3FE12358,0XE75D3033
	.word	0XBFB1973B,0XD1464000
	.word	0XBD3566D1,0X54F930B3
	.word	0,0
	.word	0X3FE11EC3,0X46E36092
	.word	0XBFB152B7,0X99BB4000
	.word	0X3D09BB29,0X7030829
	.word	0,0
	.word	0X3FE11A30,0X19A74826
	.word	0XBFB10E45,0XB3CB0000
	.word	0X3D37CF69,0X284A3465
	.word	0,0
	.word	0X3FE1159F,0X5DB29606
	.word	0XBFB0C9E6,0X15AC4000
	.word	0XBD2C2DA8,0X974D976
	.word	0,0
	.word	0X3FE11111,0X11111111
	.word	0XBFB08598,0XB59E4000
	.word	0X3D17E5DD,0X7009902C
	.word	0,0
	.word	0X3FE10C85,0X31D0952E
	.word	0XBFB0415D,0X89E74000
	.word	0XBD1111C0,0X5CF1D753
	.word	0,0
	.word	0X3FE107FB,0XBE011080
	.word	0XBFAFFA69,0X11AB8000
	.word	0XBD23008C,0X98381A8F
	.word	0,0
	.word	0X3FE10374,0XB3B480AA
	.word	0XBFAF723B,0X51800000
	.word	0X3D3D6EB0,0XDD5610D3
	.word	0,0
	.word	0X3FE0FEF0,0X10FEF011
	.word	0XBFAEEA31,0XC0068000
	.word	0XBD3C3DD8,0X3606D891
	.word	0,0
	.word	0X3FE0FA6D,0XD3F67322
	.word	0XBFAE624C,0X4A0B8000
	.word	0X3D30F25C,0X74676689
	.word	0,0
	.word	0X3FE0F5ED,0XFAB325A2
	.word	0XBFADDA8A,0XDC680000
	.word	0X3D21B1AC,0X64D9E42F
	.word	0,0
	.word	0X3FE0F170,0X834F27FA
	.word	0XBFAD52ED,0X64060000
	.word	0X3D33C85D,0X2A29BBD6
	.word	0,0
	.word	0X3FE0ECF5,0X6BE69C90
	.word	0XBFACCB73,0XCDDD8000
	.word	0XBD3965C3,0X6E09F5FE
	.word	0,0
	.word	0X3FE0E87C,0XB297A51E
	.word	0XBFAC441E,0X6F70000
	.word	0XBD354F1F,0X49850D15
	.word	0,0
	.word	0X3FE0E406,0X55826011
	.word	0XBFABBCEB,0XFC690000
	.word	0X3D17BF86,0X8C317C2A
	.word	0,0
	.word	0X3FE0DF92,0X52C8E5E6
	.word	0XBFAB35DD,0X9B588000
	.word	0XBD3D5674,0XD6CF558E
	.word	0,0
	.word	0X3FE0DB20,0XA88F4696
	.word	0XBFAAAEF2,0XD0FB0000
	.word	0XBD20FC1A,0X353BB42E
	.word	0,0
	.word	0X3FE0D6B1,0X54FB86F9
	.word	0XBFAA282B,0X8A938000
	.word	0X3D2E8F59,0X80EFC8E3
	.word	0,0
	.word	0X3FE0D244,0X56359E3A
	.word	0XBFA9A187,0XB5740000
	.word	0X3D30C22E,0X4EC4D90D
	.word	0,0
	.word	0X3FE0CDD9,0XAA677344
	.word	0XBFA91B07,0X3EFD8000
	.word	0X3D19D7C5,0X3F76CA96
	.word	0,0
	.word	0X3FE0C971,0X4FBCDA3B
	.word	0XBFA894AA,0X149F8000
	.word	0XBD39A19A,0X8BE97661
	.word	0,0
	.word	0X3FE0C50B,0X446391F3
	.word	0XBFA80E70,0X23D90000
	.word	0X3D399DC1,0X6F28BF45
	.word	0,0
	.word	0X3FE0C0A7,0X868B4171
	.word	0XBFA78859,0X5A358000
	.word	0X3D108B0D,0X83B3A4C
	.word	0,0
	.word	0X3FE0BC46,0X14657569
	.word	0XBFA70265,0XA5510000
	.word	0X3D2888DF,0X11FD5CE7
	.word	0,0
	.word	0X3FE0B7E6,0XEC259DC8
	.word	0XBFA67C94,0XF2D48000
	.word	0XBD3DAC20,0X827CCA0C
	.word	0,0
	.word	0X3FE0B38A,0XC010B39
	.word	0XBFA5F6E7,0X30790000
	.word	0X3D20485A,0X8012494C
	.word	0,0
	.word	0X3FE0AF2F,0X722EECB5
	.word	0XBFA5715C,0X4C040000
	.word	0X3D38888D,0XDFC47628
	.word	0,0
	.word	0X3FE0AAD7,0X1CE84D16
	.word	0XBFA4EBF4,0X334A0000
	.word	0X3D2D9150,0XF73BE773
	.word	0,0
	.word	0X3FE0A681,0XA6810A7
	.word	0XBFA466AE,0XD42E0000
	.word	0X3D2C1673,0X75BDFD28
	.word	0,0
	.word	0X3FE0A22D,0X38EAF2BF
	.word	0XBFA3E18C,0X1CA08000
	.word	0XBD3748ED,0X3F6E378E
	.word	0,0
	.word	0X3FE09DDB,0XA6AF8360
	.word	0XBFA35C8B,0XFAA10000
	.word	0XBD38357D,0X5EF9EB35
	.word	0,0
	.word	0X3FE0998C,0X51F624D5
	.word	0XBFA2D7AE,0X5C3C8000
	.word	0X3D322939,0X459DA66D
	.word	0,0
	.word	0X3FE0953F,0X39010954
	.word	0XBFA252F3,0X2F8D0000
	.word	0XBD283E9A,0XE021B67B
	.word	0,0
	.word	0X3FE090F4,0X5A1430AA
	.word	0XBFA1CE5A,0X62BC0000
	.word	0XBD3A9CC7,0X8D8DF999
	.word	0,0
	.word	0X3FE08CAB,0XB37565E2
	.word	0XBFA149E3,0XE4008000
	.word	0X3D32B98A,0X9A4168FD
	.word	0,0
	.word	0X3FE08865,0X436C3CF7
	.word	0XBFA0C58F,0XA19E0000
	.word	0X3D0559D1,0X58B17913
	.word	0,0
	.word	0X3FE08421,0X8421084
	.word	0XBFA0415D,0X89E78000
	.word	0X3D3DDDC7,0XF461C516
	.word	0,0
	.word	0X3FE07FDF,0X41FF7C
	.word	0XBF9F7A9B,0X16780000
	.word	0XBD242AD9,0X271BE7D7
	.word	0,0
	.word	0X3FE07B9F,0X29B8EAE2
	.word	0XBF9E72BF,0X28140000
	.word	0X3D28D751,0X49774D47
	.word	0,0
	.word	0X3FE07761,0X82F57386
	.word	0XBF9D6B27,0X25980000
	.word	0X3D39FF7B,0X50D1B838
	.word	0,0
	.word	0X3FE07326,0XA47F7C6
	.word	0XBF9C63D2,0XEC150000
	.word	0X3D35439C,0XE030A687
	.word	0,0
	.word	0X3FE06EEC,0XBE029155
	.word	0XBF9B5CC2,0X58B70000
	.word	0XBD18E611,0XB8AFBFE8
	.word	0,0
	.word	0X3FE06AB5,0X9C7912FB
	.word	0XBF9A55F5,0X48C60000
	.word	0X3D2DE070,0X9F2D03C9
	.word	0,0
	.word	0X3FE06680,0XA4010668
	.word	0XBF994F6B,0X99A20000
	.word	0XBD311D5E,0XF96CF7F5
	.word	0,0
	.word	0X3FE0624D,0XD2F1A9FC
	.word	0XBF984925,0X28C90000
	.word	0X3D2AA0BA,0X325A0C34
	.word	0,0
	.word	0X3FE05E1D,0X27A3EE9C
	.word	0XBF974321,0XD3D00000
	.word	0XBCFB4A69,0XFE94778
	.word	0,0
	.word	0X3FE059EE,0XA0727586
	.word	0XBF963D61,0X78690000
	.word	0XBD07ABF3,0X89596542
	.word	0,0
	.word	0X3FE055C2,0X3BB98E2A
	.word	0XBF9537E3,0XF45F0000
	.word	0XBD2AB259,0XD2D7F253
	.word	0,0
	.word	0X3FE05197,0XF7D73404
	.word	0XBF9432A9,0X25980000
	.word	0XBD098139,0X928637FE
	.word	0,0
	.word	0X3FE04D6F,0XD32B0C7B
	.word	0XBF932DB0,0XEA130000
	.word	0XBD2710CB,0X130895FC
	.word	0,0
	.word	0X3FE04949,0XCC1664C5
	.word	0XBF9228FB,0X1FEA0000
	.word	0XBD2713E3,0X284991FE
	.word	0,0
	.word	0X3FE04525,0XE0FC2FCB
	.word	0XBF912487,0XA5500000
	.word	0XBD3FDBE5,0XFED4B393
	.word	0,0
	.word	0X3FE04104,0X10410410
	.word	0XBF902056,0X58930000
	.word	0XBD3611D2,0X7C8E8417
	.word	0,0
	.word	0X3FE03CE4,0X584B19A0
	.word	0XBF8E38CE,0X30340000
	.word	0X3D39DE88,0XA3DA281A
	.word	0,0
	.word	0X3FE038C6,0XB78247FC
	.word	0XBF8C3173,0X84C80000
	.word	0X3D341F33,0XFCEFB9FE
	.word	0,0
	.word	0X3FE034AB,0X2C50040D
	.word	0XBF8A2A9C,0X6C180000
	.word	0X3D3F73BC,0X4D6D3472
	.word	0,0
	.word	0X3FE03091,0XB51F5E1A
	.word	0XBF882448,0XA3880000
	.word	0XBD345544,0X12C584E0
	.word	0,0
	.word	0X3FE02C7A,0X505CFFBF
	.word	0XBF861E77,0XE8B60000
	.word	0X3D38073E,0XEAF8EAF3
	.word	0,0
	.word	0X3FE02864,0XFC7729E9
	.word	0XBF841929,0XF9680000
	.word	0XBD1977C7,0X55D01368
	.word	0,0
	.word	0X3FE02451,0XB7DDB2D2
	.word	0XBF82145E,0X939E0000
	.word	0XBD3E3D12,0X38C4EA00
	.word	0,0
	.word	0X3FE02040,0X81020408
	.word	0XBF801015,0X75880000
	.word	0XBD3BCE25,0X1998B506
	.word	0,0
	.word	0X3FE01C31,0X5657186B
	.word	0XBF7C189C,0XBB100000
	.word	0X3D3D8055,0X12588560
	.word	0,0
	.word	0X3FE01824,0X36517A37
	.word	0XBF781212,0X14580000
	.word	0XBD1AD503,0X82973F27
	.word	0,0
	.word	0X3FE01419,0X1F674111
	.word	0XBF740C8A,0X74780000
	.word	0XBD1E3871,0XDF070002
	.word	0,0
	.word	0X3FE01010,0X10101010
	.word	0XBF700805,0X59580000
	.word	0XBD2166AF,0XCB31C67B
	.word	0,0
	.word	0X3FE00C09,0X6C513CF
	.word	0XBF680904,0X82880000
	.word	0XBD285C06,0X96A70C0C
	.word	0,0
	.word	0X3FE00804,0X2010080
	.word	0XBF600401,0X55D80000
	.word	0X3D33BB10,0XC7CC7089
	.word	0,0
	.word	0X3FE00401,0X401004
	.word	0XBF500200,0X55600000
	.word	0XBD356224,0XCD5F35F8
	.word	0,0
	.word	0X3FE00000,0
	.word	0,0
	.word	0,0
	.word	0,0
	.align	8
__fj_dlogexp1_const_adr_:
	.xword		__fj_dlogexp1_const_
	.section	".rodata"
	.align	128
__fj_dlogexp2_const_:
	.word	0X3FF00000,0
	.word	0,0
	.word	0X3FF00B1A,0XFA5ABCBF
	.word	0XBC84F6B2,0XA7609F71
	.word	0X3FF0163D,0XA9FB3335
	.word	0X3C9B6129,0X9AB8CDB7
	.word	0X3FF02168,0X143B0281
	.word	0XBC82BF31,0XFC54EB6
	.word	0X3FF02C9A,0X3E778061
	.word	0XBC719083,0X535B085D
	.word	0X3FF037D4,0X2E11BBCC
	.word	0X3C656811,0XEEADE11A
	.word	0X3FF04315,0XE86E7F85
	.word	0XBC90A31C,0X1977C96E
	.word	0X3FF04E5F,0X72F654B1
	.word	0X3C84C379,0X3AA0D08C
	.word	0X3FF059B0,0XD3158574
	.word	0X3C8D73E2,0XA475B465
	.word	0X3FF0650A,0XE3C1F89
	.word	0XBC95CB7B,0X5799C397
	.word	0X3FF0706B,0X29DDF6DE
	.word	0XBC8C91DF,0XE2B13C27
	.word	0X3FF07BD4,0X2B72A836
	.word	0X3C832334,0X54458700
	.word	0X3FF08745,0X18759BC8
	.word	0X3C6186BE,0X4BB284FF
	.word	0X3FF092BD,0XF66607E0
	.word	0XBC968063,0X800A3FD1
	.word	0X3FF09E3E,0XCAC6F383
	.word	0X3C914878,0X18316136
	.word	0X3FF0A9C7,0X9B1F3919
	.word	0X3C85D16C,0X873D1D37
	.word	0X3FF0B558,0X6CF9890F
	.word	0X3C98A62E,0X4ADC610B
	.word	0X3FF0C0F1,0X45E46C85
	.word	0X3C94F989,0X6D21CEF
	.word	0X3FF0CC92,0X2B7247F7
	.word	0X3C901EDC,0X16E24F71
	.word	0X3FF0D83B,0X23395DEC
	.word	0XBC9BC14D,0XE43F316A
	.word	0X3FF0E3EC,0X32D3D1A2
	.word	0X3C403A17,0X27C57B53
	.word	0X3FF0EFA5,0X5FDFA9C5
	.word	0XBC949DB9,0XBC54021B
	.word	0X3FF0FB66,0XAFFED31B
	.word	0XBC6B9BED,0XC44EBD7B
	.word	0X3FF10730,0X28D7233E
	.word	0X3C8D46EB,0X1692FDD5
	.word	0X3FF11301,0XD0125B51
	.word	0XBC96C510,0X39449B3A
	.word	0X3FF11EDB,0XAB5E2AB6
	.word	0XBC9CA454,0XF703FB72
	.word	0X3FF12ABD,0XC06C31CC
	.word	0XBC51B514,0XB36CA5C7
	.word	0X3FF136A8,0X14F204AB
	.word	0XBC67108F,0XBA48DCEF
	.word	0X3FF1429A,0XAEA92DE0
	.word	0XBC932FBF,0X9AF1369E
	.word	0X3FF14E95,0X934F312E
	.word	0XBC8B91E8,0X39BF44AB
	.word	0X3FF15A98,0XC8A58E51
	.word	0X3C82406A,0XB9EEAB0A
	.word	0X3FF166A4,0X5471C3C2
	.word	0X3C58F23B,0X82EA1A32
	.word	0X3FF172B8,0X3C7D517B
	.word	0XBC819041,0XB9D78A76
	.word	0X3FF17ED4,0X8695BBC0
	.word	0X3C709E3F,0XE2AC5A64
	.word	0X3FF18AF9,0X388C8DEA
	.word	0XBC911023,0XD1970F6C
	.word	0X3FF19726,0X58375D2F
	.word	0X3C94AADD,0X85F17E08
	.word	0X3FF1A35B,0XEB6FCB75
	.word	0X3C8E5B4C,0X7B4968E4
	.word	0X3FF1AF99,0XF8138A1C
	.word	0X3C97BF85,0XA4B69280
	.word	0X3FF1BBE0,0X84045CD4
	.word	0XBC995386,0X352EF607
	.word	0X3FF1C82F,0X95281C6B
	.word	0X3C900977,0X8010F8C9
	.word	0X3FF1D487,0X3168B9AA
	.word	0X3C9E016E,0XA2643C
	.word	0X3FF1E0E7,0X5EB44027
	.word	0XBC96FDD8,0X88CB6DE
	.word	0X3FF1ED50,0X22FCD91D
	.word	0XBC91DF98,0X27BB78C
	.word	0X3FF1F9C1,0X8438CE4D
	.word	0XBC9BF524,0XA097AF5C
	.word	0X3FF2063B,0X88628CD6
	.word	0X3C8DC775,0X814A8495
	.word	0X3FF212BE,0X3578A819
	.word	0X3C93592D,0X2CFCAAC9
	.word	0X3FF21F49,0X917DDC96
	.word	0X3C82A97E,0X9494A5EE
	.word	0X3FF22BDD,0XA27912D1
	.word	0X3C8D34FB,0X5577D69F
	.word	0X3FF2387A,0X6E756238
	.word	0X3C99B07E,0XB6C70573
	.word	0X3FF2451F,0XFB82140A
	.word	0X3C8ACFCC,0X911CA996
	.word	0X3FF251CE,0X4FB2A63F
	.word	0X3C8AC155,0XBEF4F4A4
	.word	0X3FF25E85,0X711ECE75
	.word	0X3C93E1A2,0X4AC31B2C
	.word	0X3FF26B45,0X65E27CDD
	.word	0X3C82BD33,0X9940E9D9
	.word	0X3FF2780E,0X341DDF29
	.word	0X3C9E067C,0X5F9E76C
	.word	0X3FF284DF,0XE1F56381
	.word	0XBC9A4C3A,0X8C3F0D7E
	.word	0X3FF291BA,0X7591BB70
	.word	0XBC82CC72,0X28401CBD
	.word	0X3FF29E9D,0XF51FDEE1
	.word	0X3C8612E8,0XAFAD1255
	.word	0X3FF2AB8A,0X66D10F13
	.word	0XBC995743,0X191690A7
	.word	0X3FF2B87F,0XD0DAD990
	.word	0XBC410ADC,0XD6381AA4
	.word	0X3FF2C57E,0X39771B2F
	.word	0XBC950145,0XA6EB5125
	.word	0X3FF2D285,0XA6E4030B
	.word	0X3C900247,0X54DB41D5
	.word	0X3FF2DF96,0X1F641589
	.word	0X3C9D16CF,0XFBBCE198
	.word	0X3FF2ECAF,0XA93E2F56
	.word	0X3C71CA0F,0X45D52383
	.word	0X3FF2F9D2,0X4ABD886B
	.word	0XBC653C55,0X532BDA93
	.word	0X3FF306FE,0XA31B715
	.word	0X3C86F46A,0XD23182E4
	.word	0X3FF31432,0XEDEEB2FD
	.word	0X3C8959A3,0XF3F3FCD1
	.word	0X3FF32170,0XFC4CD831
	.word	0X3C8A9CE7,0X8E18047C
	.word	0X3FF32EB8,0X3BA8EA31
	.word	0X3CA1DD0B,0XE1A58674
	.word	0X3FF33C08,0XB26416FF
	.word	0X3C932721,0X843659A6
	.word	0X3FF34962,0X66E3FA2D
	.word	0XBC835A75,0X930881A4
	.word	0X3FF356C5,0X5F929FF1
	.word	0XBC8B5CEE,0X5C4E4628
	.word	0X3FF36431,0XA2DE883B
	.word	0XBC8C3144,0XA06CB85D
	.word	0X3FF371A7,0X373AA9CB
	.word	0XBC963AEA,0XBF42EAE2
	.word	0X3FF37F26,0X231E754A
	.word	0XBC99F5CA,0X9ECEB23C
	.word	0X3FF38CAE,0X6D05D866
	.word	0XBC9E958D,0X3C9904BD
	.word	0X3FF39A40,0X1B7140EF
	.word	0XBC99A9A5,0XFC8E2934
	.word	0X3FF3A7DB,0X34E59FF7
	.word	0XBC75E436,0XD661F5E3
	.word	0X3FF3B57F,0XBFEC6CF4
	.word	0X3C954C66,0XE26FFF18
	.word	0X3FF3C32D,0XC313A8E5
	.word	0XBC9EFFF8,0X375D29C4
	.word	0X3FF3D0E5,0X44EDE173
	.word	0X3C7FE8D0,0X8C284C71
	.word	0X3FF3DEA6,0X4C123422
	.word	0X3C8ADA09,0X11F09EBC
	.word	0X3FF3EC70,0XDF1C5175
	.word	0XBC8AF663,0X7B8C9BCA
	.word	0X3FF3FA45,0X4AC801C
	.word	0XBC97D023,0XF956F9F3
	.word	0X3FF40822,0XC367A024
	.word	0X3C8BDDF8,0XB6F4D048
	.word	0X3FF4160A,0X21F72E2A
	.word	0XBC5EF369,0X1C309278
	.word	0X3FF423FB,0X2709468A
	.word	0XBC98462D,0XC0B314DD
	.word	0X3FF431F5,0XD950A897
	.word	0XBC81C7DD,0XE35F7999
	.word	0X3FF43FFA,0X3F84B9D4
	.word	0X3C8880BE,0X9704C003
	.word	0X3FF44E08,0X6061892D
	.word	0X3C489B7A,0X4EF80D0
	.word	0X3FF45C20,0X42A7D232
	.word	0XBC686419,0X82FB1F8E
	.word	0X3FF46A41,0XED1D0057
	.word	0X3C9C944B,0XD1648A76
	.word	0X3FF4786D,0X668B3237
	.word	0XBC9C20F0,0XED445733
	.word	0X3FF486A2,0XB5C13CD0
	.word	0X3C73C1A3,0XB69062F0
	.word	0X3FF494E1,0XE192AED2
	.word	0XBC83B289,0X5E499EA0
	.word	0X3FF4A32A,0XF0D7D3DF
	.word	0XBCA31A4E,0X861720D5
	.word	0X3FF4B17D,0XEA6DB7D7
	.word	0XBC8125B8,0X7F2897F0
	.word	0X3FF4BFDA,0XD5362A27
	.word	0X3C7D4397,0XAFEC42E2
	.word	0X3FF4CE41,0XB817C114
	.word	0X3C905E29,0X690ABD5D
	.word	0X3FF4DCB2,0X99FDDD0D
	.word	0X3C98ECDB,0XBC6A7833
	.word	0X3FF4EB2D,0X81D8ABFF
	.word	0XBC95257D,0X2E5D7A52
	.word	0X3FF4F9B2,0X769D2CA7
	.word	0XBC94B309,0XD25957E3
	.word	0X3FF50841,0X7F4531EE
	.word	0X3C7A249B,0X49B7465F
	.word	0X3FF516DA,0XA2CF6642
	.word	0XBC8F7685,0X69BD93EF
	.word	0X3FF5257D,0XE83F4EEF
	.word	0XBC7C998D,0X43EFEF71
	.word	0X3FF5342B,0X569D4F82
	.word	0XBC807ABE,0X1DB13CAD
	.word	0X3FF542E2,0XF4F6AD27
	.word	0X3C87926D,0X192D5F7E
	.word	0X3FF551A4,0XCA5D920F
	.word	0XBC8D689C,0XEFEDE59B
	.word	0X3FF56070,0XDDE910D2
	.word	0XBC90FB6E,0X168EEBF0
	.word	0X3FF56F47,0X36B527DA
	.word	0X3C99BB2C,0X11D93AD
	.word	0X3FF57E27,0XDBE2C4CF
	.word	0XBC90B98C,0X8A57B9C4
	.word	0X3FF58D12,0XD497C7FD
	.word	0X3C8295E1,0X5B9A1DE8
	.word	0X3FF59C08,0X27FF07CB
	.word	0X3CA40E98,0X8DCC0CF9
	.word	0X3FF5AB07,0XDD485429
	.word	0X3C96324C,0X54647AD
	.word	0X3FF5BA11,0XFBA87A03
	.word	0XBC9B77A1,0X4C233E1A
	.word	0X3FF5C926,0X8A5946B7
	.word	0X3C3C4B1B,0X816986A2
	.word	0X3FF5D845,0X90998B93
	.word	0XBC9CD6A7,0XA8B45643
	.word	0X3FF5E76F,0X15AD2148
	.word	0X3C9BA6F9,0X3080E65E
	.word	0X3FF5F6A3,0X20DCEB71
	.word	0XBC89EADD,0XE3CDCF92
	.word	0X3FF605E1,0XB976DC09
	.word	0XBC93E242,0X9B56DE47
	.word	0X3FF6152A,0XE6CDF6F4
	.word	0X3C9E4B3E,0X4AB84C27
	.word	0X3FF6247E,0XB03A5585
	.word	0XBC9383C1,0X7E40B497
	.word	0X3FF633DD,0X1D1929FD
	.word	0X3C984710,0XBEB964E5
	.word	0X3FF64346,0X34CCC320
	.word	0XBC8C483C,0X759D8933
	.word	0X3FF652B9,0XFEBC8FB7
	.word	0XBC9AE3D5,0XC9A73E09
	.word	0X3FF66238,0X82552225
	.word	0XBC9BB609,0X87591C34
	.word	0X3FF671C1,0XC70833F6
	.word	0XBC8E8732,0X586C6133
	.word	0X3FF68155,0XD44CA973
	.word	0X3C6038AE,0X44F73E65
	.word	0X3FF690F4,0XB19E9538
	.word	0X3C8804BD,0X9AEB445D
	.word	0X3FF6A09E,0X667F3BCC
	.word	0X3CA21165,0XF626CDD5
	.word	0X3FF6B052,0XFA75173E
	.word	0X3C7A38F5,0X2C9A9D0E
	.word	0X3FF6C012,0X750BDABF
	.word	0XBC728956,0X67FF0B0D
	.word	0X3FF6CFDC,0XDDD47646
	.word	0XBCA1C2AB,0X2487467C
	.word	0X3FF6DFB2,0X3C651A2F
	.word	0XBC6BBE3A,0X683C88AB
	.word	0X3FF6EF92,0X98593AE5
	.word	0XBC90B974,0X9E1AC8B2
	.word	0X3FF6FF7D,0XF9519484
	.word	0XBC883C0F,0X25860EF6
	.word	0X3FF70F74,0X66F42E87
	.word	0X3C59D644,0XD45AA65F
	.word	0X3FF71F75,0XE8EC5F74
	.word	0XBC816E47,0X86887A99
	.word	0X3FF72F82,0X86EAD08A
	.word	0XBC920AA0,0X2CD62C72
	.word	0X3FF73F9A,0X48A58174
	.word	0XBC90A8D9,0X6C65D53C
	.word	0X3FF74FBD,0X35D7CBFD
	.word	0X3C9047FD,0X618A6E1C
	.word	0X3FF75FEB,0X564267C9
	.word	0XBC902459,0X57316DD3
	.word	0X3FF77024,0XB1AB6E0A
	.word	0XBCA243C4,0X74B75C04
	.word	0X3FF78069,0X4FDE5D3F
	.word	0X3C9866B8,0XA02162D
	.word	0X3FF790B9,0X38AC1CF6
	.word	0X3C9349A8,0X62AADD3E
	.word	0X3FF7A114,0X73EB0187
	.word	0XBC841577,0XEE04992F
	.word	0X3FF7B17B,0X976CFDB
	.word	0XBC9BEBB5,0X8468DC88
	.word	0X3FF7C1ED,0X130C133
	.word	0XBCA076D9,0X9774D915
	.word	0X3FF7D26A,0X62FF86F0
	.word	0X3C91BDDB,0XFB72B8B4
	.word	0X3FF7E2F3,0X36CF4E62
	.word	0X3C705D02,0XBA15797E
	.word	0X3FF7F387,0X8491C491
	.word	0XBC807F11,0XCF9311AE
	.word	0X3FF80427,0X543E1A12
	.word	0XBC927C86,0X626D972B
	.word	0X3FF814D2,0XADD106D9
	.word	0X3C946437,0XD151D4D
	.word	0X3FF82589,0X994CCE13
	.word	0XBC9D4C1D,0XD41532D8
	.word	0X3FF8364C,0X1EB941F7
	.word	0X3C999B9A,0X31DF2BD5
	.word	0X3FF8471A,0X4623C7AD
	.word	0XBC88D684,0XA341CDFB
	.word	0X3FF857F4,0X179F5B21
	.word	0XBC5BA748,0XF8B216CF
	.word	0X3FF868D9,0X9B4492ED
	.word	0XBC9FC6F8,0X9BD4F6BA
	.word	0X3FF879CA,0XD931A436
	.word	0X3C85D2D7,0XD2DB47BD
	.word	0X3FF88AC7,0XD98A6699
	.word	0X3C9994C2,0XF37CB53A
	.word	0X3FF89BD0,0XA478580F
	.word	0X3C9D5395,0X4475202B
	.word	0X3FF8ACE5,0X422AA0DB
	.word	0X3C96E9F1,0X56864B27
	.word	0X3FF8BE05,0XBAD61779
	.word	0XBCA09A50,0X81DE5DC9
	.word	0X3FF8CF32,0X16B5448C
	.word	0XBC70D55E,0X32E9E3AA
	.word	0X3FF8E06A,0X5E0866D9
	.word	0XBC97114A,0X6FC9B2E6
	.word	0X3FF8F1AE,0X99157736
	.word	0X3C85CC13,0XA2E3976C
	.word	0X3FF902FE,0XD0282C8A
	.word	0X3C9592CA,0X85FE3FD2
	.word	0X3FF9145B,0XB91FFC5
	.word	0X3CA114C3,0X68D3ED6E
	.word	0X3FF925C3,0X53AA2FE2
	.word	0XBC83455F,0XA639DB7F
	.word	0X3FF93737,0XB0CDC5E5
	.word	0XBC675FC7,0X81B57EBB
	.word	0X3FF948B8,0X2B5F98E5
	.word	0XBC8DC3D6,0X797D2D99
	.word	0X3FF95A44,0XCBC8520F
	.word	0XBC764B7C,0X96A5F039
	.word	0X3FF96BDD,0X9A7670B3
	.word	0XBC5BA596,0X7F19C895
	.word	0X3FF97D82,0X9FDE4E4F
	.word	0X3CA173D2,0X41F23D18
	.word	0X3FF98F33,0XE47A22A2
	.word	0X3C7CABDA,0XA24C78ED
	.word	0X3FF9A0F1,0X70CA07BA
	.word	0XBC9173BD,0X91CEE632
	.word	0X3FF9B2BB,0X4D53FE0C
	.word	0X3CA113D8,0XD9049574
	.word	0X3FF9C491,0X82A3F090
	.word	0X3C7C7C46,0XB071F2BE
	.word	0X3FF9D674,0X194BB8D5
	.word	0XBC9516BE,0XA3DD8233
	.word	0X3FF9E863,0X19E32323
	.word	0X3C7824CA,0X78E64C6E
	.word	0X3FF9FA5E,0X8D07F29E
	.word	0XBC84A9CE,0XAAF1FACE
	.word	0X3FFA0C66,0X7B5DE565
	.word	0XBC935949,0X5D1CD533
	.word	0X3FFA1E7A,0XED8EB8BB
	.word	0X3C9C6618,0XEE8BE70E
	.word	0X3FFA309B,0XEC4A2D33
	.word	0X3C96305C,0X7DDC36AB
	.word	0X3FFA42C9,0X80460AD8
	.word	0XBC9AA780,0X589FB120
	.word	0X3FFA5503,0XB23E255D
	.word	0XBC9D2F6E,0XDB8D41E1
	.word	0X3FFA674A,0X8AF46052
	.word	0X3C650F56,0X30670366
	.word	0X3FFA799E,0X1330B359
	.word	0XBCA21A40,0X9A9D4E1D
	.word	0X3FFA8BFE,0X53C12E59
	.word	0XBC94F867,0XB2BA15A9
	.word	0X3FFA9E6B,0X5579FDC0
	.word	0XBCA7829B,0X78840167
	.word	0X3FFAB0E5,0X21356EBB
	.word	0XBCA9D8F3,0X8945AEAF
	.word	0X3FFAC36B,0XBFD3F37A
	.word	0XBC8F9234,0XCAE76CD0
	.word	0X3FFAD5FF,0X3A3C2774
	.word	0X3C97EF3B,0XB6B1B8E5
	.word	0X3FFAE89F,0X995AD3AD
	.word	0X3C97A1CD,0X345DCC81
	.word	0X3FFAFB4C,0XE622F2FE
	.word	0X3CA5A681,0XF8675099
	.word	0X3FFB0E07,0X298DB665
	.word	0X3CA21085,0X59BF8DEE
	.word	0X3FFB20CE,0X6C9A8952
	.word	0X3C94DD02,0X4A0756CC
	.word	0X3FFB33A2,0XB84F15FB
	.word	0XBC62805E,0X3084D708
	.word	0X3FFB4684,0X15B749B1
	.word	0XBC7F763D,0XE9DF7C90
	.word	0X3FFB5972,0X8DE5593A
	.word	0XBC9C71DF,0XBBBA6DE3
	.word	0X3FFB6C6E,0X29F1C52A
	.word	0X3C92A8F3,0X52883F6E
	.word	0X3FFB7F76,0XF2FB5E47
	.word	0XBC75584F,0X7E54AC3B
	.word	0X3FFB928C,0XF22749E4
	.word	0XBC9B7216,0X54CB65C6
	.word	0X3FFBA5B0,0X30A1064A
	.word	0XBC9EFCD3,0XE54292E
	.word	0X3FFBB8E0,0XB79A6F1F
	.word	0XBC3F52D1,0XC9696205
	.word	0X3FFBCC1E,0X904BC1D2
	.word	0X3C823DD0,0X7A2D9E84
	.word	0X3FFBDF69,0XC3F3A207
	.word	0XBC3C2623,0X60EA5B53
	.word	0X3FFBF2C2,0X5BD71E08
	.word	0X3CA0811A,0XE04A31C6
	.word	0X3FFC0628,0X6141B33C
	.word	0X3CA89D69,0X57810D73
	.word	0X3FFC199B,0XDD85529C
	.word	0X3C811065,0X895048DD
	.word	0X3FFC2D1C,0XD9FA652B
	.word	0X3CA48D74,0XF41BAD14
	.word	0X3FFC40AB,0X5FFFD07A
	.word	0X3C9B4537,0XE083C60A
	.word	0X3FFC5447,0X78FAFB22
	.word	0X3C912F07,0X2493B5AF
	.word	0X3FFC67F1,0X2E57D14B
	.word	0X3C92884D,0XFF483CAD
	.word	0X3FFC7BA8,0X8988C933
	.word	0XBC8E76BB,0XBE255559
	.word	0X3FFC8F6D,0X9406E7B5
	.word	0X3C71ACBC,0X48805C44
	.word	0X3FFCA340,0X5751C4DB
	.word	0XBC87F2BE,0XD10D08F5
	.word	0X3FFCB720,0XDCEF9069
	.word	0X3C7503CB,0XD1E949DB
	.word	0X3FFCCB0F,0X2E6D1675
	.word	0XBC7D220F,0X86009093
	.word	0X3FFCDF0B,0X555DC3FA
	.word	0XBC8DD83B,0X53829D72
	.word	0X3FFCF315,0X5B5BAB74
	.word	0XBC9A08E9,0XB86DFF57
	.word	0X3FFD072D,0X4A07897C
	.word	0XBC9CBC37,0X43797A9C
	.word	0X3FFD1B53,0X2B08C969
	.word	0XBCA554E4,0XEF32E489
	.word	0X3FFD2F87,0X80D89F1
	.word	0X3CA15BC2,0X47313D44
	.word	0X3FFD43C8,0XEACAA1D6
	.word	0X3C93DB53,0XBF5A1614
	.word	0X3FFD5818,0XDCFBA487
	.word	0X3C82ED02,0XD75B3707
	.word	0X3FFD6C76,0XE862E6D3
	.word	0X3C5FE87A,0X4A8165A0
	.word	0X3FFD80E3,0X16C98398
	.word	0XBC911EC1,0X8BEDDFE8
	.word	0X3FFD955D,0X71FF6075
	.word	0X3C9A052D,0XBB9AF6BE
	.word	0X3FFDA9E6,0X3DB3285
	.word	0X3C9C2300,0X696DB532
	.word	0X3FFDBE7C,0XD63A8315
	.word	0XBC9B76F1,0X926B8BE4
	.word	0X3FFDD321,0XF301B460
	.word	0X3C92DA57,0X78F018C3
	.word	0X3FFDE7D5,0X641C0658
	.word	0XBC9CA552,0X8E79BA8F
	.word	0X3FFDFC97,0X337B9B5F
	.word	0XBC91A5CD,0X4F184B5C
	.word	0X3FFE1167,0X6B197D17
	.word	0XBC72B529,0XBD5C7F44
	.word	0X3FFE2646,0X14F5A129
	.word	0XBC97B627,0X817A1496
	.word	0X3FFE3B33,0X3B16EE12
	.word	0XBC99F4A4,0X31FDC68B
	.word	0X3FFE502E,0XE78B3FF7
	.word	0XBCAB185D,0X9FD58CDC
	.word	0X3FFE6539,0X24676D75
	.word	0X3CAA7001,0XE2B75233
	.word	0X3FFE7A51,0XFBC74C83
	.word	0X3C92D522,0XCA0C8DE2
	.word	0X3FFE8F79,0X77CDB73F
	.word	0X3CA77BB5,0XBFA7D5A8
	.word	0X3FFEA4AF,0XA2A490D9
	.word	0X3CA0B1EE,0X7431EBB6
	.word	0X3FFEB9F4,0X867CCA6E
	.word	0X3C94832F,0X2293E4F2
	.word	0X3FFECF48,0X2D8E67F0
	.word	0X3CA1B606,0X25F7293A
	.word	0X3FFEE4AA,0XA2188510
	.word	0X3C91C68D,0XA487568D
	.word	0X3FFEFA1B,0XEE615A27
	.word	0X3C9DC7F4,0X86A4B6B0
	.word	0X3FFF0F9C,0X1CB6412A
	.word	0XBC932200,0X65181D45
	.word	0X3FFF252B,0X376BBA97
	.word	0X3C93A1A5,0XBF0D8E43
	.word	0X3FFF3AC9,0X48DD7274
	.word	0XBC795A5A,0X3ED837DE
	.word	0X3FFF5076,0X5B6E4541
	.word	0XBCA3160F,0X6913AF3A
	.word	0X3FFF6632,0X798844F9
	.word	0XBCA02E42,0X656365E1
	.word	0X3FFF7BFD,0XAD9CBE14
	.word	0XBC9DBB12,0XD006350A
	.word	0X3FFF91D8,0X2243C89
	.word	0XBC612EA8,0XA779F689
	.word	0X3FFFA7C1,0X819E90D8
	.word	0X3C874853,0XF3A5931E
	.word	0X3FFFBDBA,0X3692D513
	.word	0X3CACD311,0X9D5ECE29
	.word	0X3FFFD3C2,0X2B8F71F1
	.word	0X3C62EB74,0X966579E7
	.word	0X3FFFE9D9,0X6B2A23D9
	.word	0X3C74A603,0X7442FDE3
	.align	8
__fj_dlogexp2_const_adr_:
	.xword		__fj_dlogexp2_const_
	.global	__gxx_personality_v0
	.section	".eh_frame",#alloc
.LLframe1:
	.uaword	.LLECIE1-.LLSCIE1	! CIE Length
.LLSCIE1:
	.uaword	0x0	! CIE ID
	.byte	0x1	! CIE Version
	.asciz	"zPR"	! CIE Augmentation
	.uleb128	0x1	! CIE Code Alignment Factor
	.sleb128	-8	! CIE Data Alignment Factor
	.byte	0xf
	.uleb128	0xa	! CIE Augmentation Section Length 
	.byte	0x0	! Personality Routine Encoding Specifier ( absptr )
	.uaxword	__gxx_personality_v0	! Personality Routine Name
	.byte	0x1b	! FDE Code Encoding Specifier ( pcrel | sdata4 )
	.byte	0xc	! DW_CFA_def_cfa
	.uleb128	0xe
	.uleb128	0x7ff
	.align	8	! CIE Padding
.LLECIE1:
.LLSFDE7:
	.uaword	.LLEFDE7-.LLASFDE7	! FDE Length
.LLASFDE7:
	.uaword	.LLASFDE7-.LLframe1	! FDE CIE Pointer
	.uaword	%r_disp32(.LLFB4)	! FDE Initial Location
	.uaword	.LLFE4-.LLFB4	! FDE Address Range
	.uleb128	0x0	! FDE Augmentation Section Length 
	.byte	0x4	! DW_CFA_advance_loc4
	.uaword	.LLCFI0-.LLFB4
	.byte	0xd	! DW_CFA_def_cfa_register
	.uleb128	0x1e
	.byte	0x2d	! DW_CFA_GNU_window_save
	.byte	0x9	! DW_CFA_register
	.uleb128	0xf
	.uleb128	0x1f
	.align	8	! FDE Padding
.LLEFDE7:
	.section	".rodata"
	.align	8
.LR0.cnt.1:
	.word	0X3FD24924,0X92492492
	.type	.LR0.cnt.1,#object
	.size	.LR0.cnt.1,.-.LR0.cnt.1
	.section	".rodata"
	.align	8
.LR0.cnt.2:
	.word	0X3FD55555,0X55555555
	.type	.LR0.cnt.2,#object
	.size	.LR0.cnt.2,.-.LR0.cnt.2
	.section	".rodata"
	.align	8
.LR0.cnt.3:
	.word	0X3FF33333,0X33333333
	.type	.LR0.cnt.3,#object
	.size	.LR0.cnt.3,.-.LR0.cnt.3
	.section	".rodata"
	.align	8
.LR0.cnt.4:
	.word	0X3FE33333,0X33333333
	.type	.LR0.cnt.4,#object
	.size	.LR0.cnt.4,.-.LR0.cnt.4
	.section	".rodata"
	.align	8
.LR0.cnt.5:
	.word	0X3C9CD2B2,0X97D889BC
	.type	.LR0.cnt.5,#object
	.size	.LR0.cnt.5,.-.LR0.cnt.5
	.section	".rodata"
	.align	8
.LR0.cnt.6:
	.word	0,0
	.type	.LR0.cnt.6,#object
	.size	.LR0.cnt.6,.-.LR0.cnt.6
	.section	".rodata"
	.align	8
.LR0.cnt.7:
	.word	0X3FF00000,0
	.type	.LR0.cnt.7,#object
	.size	.LR0.cnt.7,.-.LR0.cnt.7
	.section	".rodata"
	.align	8
.LR0.cnt.8:
	.word	0X3FE00000,0
	.type	.LR0.cnt.8,#object
	.size	.LR0.cnt.8,.-.LR0.cnt.8
	.section	".rodata"
	.align	8
.LR0.cnt.9:
	.word	0X40080000,0
	.type	.LR0.cnt.9,#object
	.size	.LR0.cnt.9,.-.LR0.cnt.9
	.section	".rodata"
	.align	8
.LR0.cnt.10:
	.word	0X40200000,0
	.type	.LR0.cnt.10,#object
	.size	.LR0.cnt.10,.-.LR0.cnt.10
	.section	".rodata"
	.align	8
.LR0.cnt.11:
	.word	0X40700000,0
	.type	.LR0.cnt.11,#object
	.size	.LR0.cnt.11,.-.LR0.cnt.11
	.section	".rodata"
	.align	8
.LR0.cnt.12:
	.word	0X40100000,0
	.type	.LR0.cnt.12,#object
	.size	.LR0.cnt.12,.-.LR0.cnt.12
	.section	".rodata"
	.align	8
.LR0.cnt.13:
	.word	0X400C0000,0
	.type	.LR0.cnt.13,#object
	.size	.LR0.cnt.13,.-.LR0.cnt.13
	.section	".rodata"
	.align	8
.LR0.cnt.14:
	.word	0X40390000,0
	.type	.LR0.cnt.14,#object
	.size	.LR0.cnt.14,.-.LR0.cnt.14
	.section	".rodata"
	.align	8
.LR0.cnt.15:
	.word	0X40400000,0
	.type	.LR0.cnt.15,#object
	.size	.LR0.cnt.15,.-.LR0.cnt.15
	.section	".rodata"
	.align	8
.LR0.cnt.16:
	.word	0X40490000,0
	.type	.LR0.cnt.16,#object
	.size	.LR0.cnt.16,.-.LR0.cnt.16
	.section	".rodata"
	.align	8
.LR0.cnt.17:
	.word	0X40580000,0
	.type	.LR0.cnt.17,#object
	.size	.LR0.cnt.17,.-.LR0.cnt.17
	.section	".rodata"
	.align	8
.LR0.cnt.18:
	.word	0X40500000,0
	.type	.LR0.cnt.18,#object
	.size	.LR0.cnt.18,.-.LR0.cnt.18
	.section	".rodata"
	.align	8
.LR0.cnt.19:
	.word	0X40690000,0
	.type	.LR0.cnt.19,#object
	.size	.LR0.cnt.19,.-.LR0.cnt.19
	.section	".rodata"
	.align	8
.LR0.cnt.20:
	.word	0X40955400,0
	.type	.LR0.cnt.20,#object
	.size	.LR0.cnt.20,.-.LR0.cnt.20
	.section	".rodata"
	.align	8
.LR0.cnt.21:
	.word	0X3FD00000,0
	.type	.LR0.cnt.21,#object
	.size	.LR0.cnt.21,.-.LR0.cnt.21
	.section	".rodata"
	.align	8
.LR0.cnt.22:
	.word	0X40355400,0
	.type	.LR0.cnt.22,#object
	.size	.LR0.cnt.22,.-.LR0.cnt.22
	.section	".data"
	.align	8
_ZN6kernel14normalzeFactorE:
	.skip	8
	.type	_ZN6kernel14normalzeFactorE,#object
	.size	_ZN6kernel14normalzeFactorE,.-_ZN6kernel14normalzeFactorE
	.section	".data"
	.align	8
_ZN6kernel2piE:
	.skip	8
	.type	_ZN6kernel2piE,#object
	.size	_ZN6kernel2piE,.-_ZN6kernel2piE
	.section	".data"
	.align	8
_ZN4math3NaNE:
	.skip	8
	.type	_ZN4math3NaNE,#object
	.size	_ZN4math3NaNE,.-_ZN4math3NaNE
	.section	".data"
	.align	8
_ZN4math2piE:
	.skip	8
	.type	_ZN4math2piE,#object
	.size	_ZN4math2piE,.-_ZN4math2piE
	.section	".data"
	.align	16
.LS0:
	.align	8
.LS0.cnt.23:
	.word	1065353216
	.type	.LS0.cnt.23,#object
	.size	.LS0.cnt.23,.-.LS0.cnt.23
	.align	8
.LS0.cnt.24:
	.word	0XBFDFFFFF,0XFFFFF967
	.type	.LS0.cnt.24,#object
	.size	.LS0.cnt.24,.-.LS0.cnt.24
	.align	8
.LS0.cnt.25:
	.word	0X3FD55555,0X55554BC9
	.type	.LS0.cnt.25,#object
	.size	.LS0.cnt.25,.-.LS0.cnt.25
	.align	8
.LS0.cnt.26:
	.word	0XBFD00001,0XC6F05DB
	.type	.LS0.cnt.26,#object
	.size	.LS0.cnt.26,.-.LS0.cnt.26
	.align	8
.LS0.cnt.27:
	.word	0X3FC9999B,0XC4D3821A
	.type	.LS0.cnt.27,#object
	.size	.LS0.cnt.27,.-.LS0.cnt.27
	.align	8
.LS0.cnt.28:
	.word	0X3FF00400,0
	.type	.LS0.cnt.28,#object
	.size	.LS0.cnt.28,.-.LS0.cnt.28
	.align	8
.LS0.cnt.29:
	.word	0X3FE62E42,0XFE000000
	.type	.LS0.cnt.29,#object
	.size	.LS0.cnt.29,.-.LS0.cnt.29
	.align	8
.LS0.cnt.30:
	.word	0X3E1F473D,0XE6AF278F
	.type	.LS0.cnt.30,#object
	.size	.LS0.cnt.30,.-.LS0.cnt.30
	.align	8
.LS0.cnt.31:
	.word	0X80400,0
	.type	.LS0.cnt.31,#object
	.size	.LS0.cnt.31,.-.LS0.cnt.31
	.align	8
.LS0.cnt.32:
	.word	0,0X1000
	.type	.LS0.cnt.32,#object
	.size	.LS0.cnt.32,.-.LS0.cnt.32
	.align	8
.LS0.cnt.33:
	.word	0X7FF,0XFFFFFFFF
	.type	.LS0.cnt.33,#object
	.size	.LS0.cnt.33,.-.LS0.cnt.33
	.align	8
.LS0.cnt.34:
	.word	0,0X1
	.type	.LS0.cnt.34,#object
	.size	.LS0.cnt.34,.-.LS0.cnt.34
	.align	8
.LS0.cnt.35:
	.word	0X40900000,0
	.type	.LS0.cnt.35,#object
	.size	.LS0.cnt.35,.-.LS0.cnt.35
	.align	8
.LS0.cnt.36:
	.word	0XC08FF800,0
	.type	.LS0.cnt.36,#object
	.size	.LS0.cnt.36,.-.LS0.cnt.36
	.align	8
.LS0.cnt.37:
	.word	0X3FEFFFFF,0XFFFFFFB1
	.type	.LS0.cnt.37,#object
	.size	.LS0.cnt.37,.-.LS0.cnt.37
	.align	8
.LS0.cnt.38:
	.word	0X3FE00000,0X3B
	.type	.LS0.cnt.38,#object
	.size	.LS0.cnt.38,.-.LS0.cnt.38
	.align	8
.LS0.cnt.39:
	.word	0X3FC55555,0X7E548294
	.type	.LS0.cnt.39,#object
	.size	.LS0.cnt.39,.-.LS0.cnt.39
	.align	8
.LS0.cnt.40:
	.word	0X3FA55555,0X44670F74
	.type	.LS0.cnt.40,#object
	.size	.LS0.cnt.40,.-.LS0.cnt.40
	.align	8
.LS0.cnt.41:
	.word	0X1000000,0
	.type	.LS0.cnt.41,#object
	.size	.LS0.cnt.41,.-.LS0.cnt.41
	.align	8
.LS0.cnt.42:
	.word	0X100000,0
	.type	.LS0.cnt.42,#object
	.size	.LS0.cnt.42,.-.LS0.cnt.42
	.align	8
.LS0.cnt.43:
	.word	0X80000000,0
	.type	.LS0.cnt.43,#object
	.size	.LS0.cnt.43,.-.LS0.cnt.43
	.align	8
.LS0.cnt.44:
	.word	0X3FF71547,0X652B82FE
	.type	.LS0.cnt.44,#object
	.size	.LS0.cnt.44,.-.LS0.cnt.44
	.align	8
.LS0.cnt.45:
	.word	0XBF700000,0
	.type	.LS0.cnt.45,#object
	.size	.LS0.cnt.45,.-.LS0.cnt.45
	.align	8
.LS0.cnt.46:
	.word	0X408FF800,0
	.type	.LS0.cnt.46,#object
	.size	.LS0.cnt.46,.-.LS0.cnt.46
	.align	8
.LS0.cnt.47:
	.word	0XFFFFFFFF,0
	.type	.LS0.cnt.47,#object
	.size	.LS0.cnt.47,.-.LS0.cnt.47
	.align	8
.LS0.cnt.48:
	.word	0X43400000,0
	.type	.LS0.cnt.48,#object
	.size	.LS0.cnt.48,.-.LS0.cnt.48
	.align	8
.LS0.cnt.49:
	.word	0X7FF00000,0
	.type	.LS0.cnt.49,#object
	.size	.LS0.cnt.49,.-.LS0.cnt.49
	.align	8
.LS0.cnt.50:
	.word	0X7FF80000,0
	.type	.LS0.cnt.50,#object
	.size	.LS0.cnt.50,.-.LS0.cnt.50
