	.set noreorder
	.set noat
	#  /usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-495/be::5.421-sw-495

	#-----------------------------------------------------------
	# Compiling master.c (/tmp/ccI#.75SDZI)
	#-----------------------------------------------------------

	#-----------------------------------------------------------
	# Options:
	#-----------------------------------------------------------
	#  Target:SW2, ISA:ISA_1, Endian:little, Pointer Size:64
	#  -O2	(Optimization level)
	#  -g0	(Debug level)
	#  -m2	(Report advisories)
	#-----------------------------------------------------------

	.file	1	"/home/export/online1/riken/sandbox/nitadori_box/ring/nitadori_box/icache/master.c"
	.file	2	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/stdio.h"
	.file	3	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdlib.h"
	.file	4	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/sys/sysmacros.h"
	.file	5	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/mathinline.h"
	.file	6	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigset.h"
	.file	7	"/usr/sw-mpp/swcc/swgcc-binary-sw_64.1206/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/string2.h"

	.comm	sw3_sig_handler, 512, 32
	.comm	sw3_exp_handler, 256, 32
	.comm	START_ADDR, 8, 8

	.section .text, "ax", "progbits"
	.align	4

	.section .rodata, "a", "progbits"
	.align	3
	.comm	sw3_sig_handler, 512, 32
	.comm	sw3_exp_handler, 256, 32
	.comm	START_ADDR, 8, 8
	.section .text

	# Program Unit: main
	.align 4
	.ent	main#
	.globl	main
main:	# 0x0
	# anon32 = 32
	# return_address = 0
	# _temp_gra_spill2 = 40
	.loc	1	14	0
#  10  }
#  11  
#  12  extern void SLAVE_FUN(func)();
#  13  
#  14  int main() {
.BB1_main:
	.prologue
	ldih	$gp,0($27)               	!gpdisp!1	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!1	# [1] 0
	
$ng..main:
	ldi	$sp,-64($sp)              	# [3] 
.LCFI_main_adjustsp:
	stl	$26,0($sp)                	# [4] return_address
.LCFI_main_store26:
	.loc	1	15	0
#  15  	athread_init();
	.globl	athread_init
	ldl	$27,athread_init($gp)     	!literal	# [4] athread_init
	.frame $15,64,$26,0
	.mask 0x4008000,-64
	call	$26,($27),athread_init   	# [7] athread_init
.BB2_main:
	ldih	$gp,0($26)               	!gpdisp!2	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!2	# [1] 0
	rtc $0
	.loc	1	19	0
#  16  	{
#  17  		int j;
#  18  		long t0 = rpcc_();
#  19  		for(j=0; j<1000; j++){
	stw	$31,32($sp)               	# [3] anon32
	.align	4
.Lt_43_2050:
#<loop> Loop body line 19, nesting depth: 1, iterations: 1000
	.loc	1	20	0
#  20  			athread_spawn(func, 0);
	.globl	__real_athread_spawn
	ldl	$27,__real_athread_spawn($gp)	!literal	# [0] __real_athread_spawn
	.globl	slave_func
	ldl	$16,slave_func($gp)       	!literal	# [0] slave_func
	mov	$31,$17                   	# [0] 
	call	$26,($27),__real_athread_spawn 	# [3] __real_athread_spawn
.BB5_main:
#<loop> Part of loop body line 65536, head labeled .Lt_43_2050
	ldih	$gp,0($26)               	!gpdisp!3	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!3	# [1] 0
	.loc	1	21	0
#  21  			athread_join();
	.globl	athread_join
	ldl	$27,athread_join($gp)     	!literal	# [2] athread_join
	call	$26,($27),athread_join   	# [5] athread_join
.BB6_main:
#<loop> Part of loop body line 65536, head labeled .Lt_43_2050
	.loc	1	19	0
	ldw	$0,32($sp)                	# [0] anon32
	ldi	$1,1000($31)              	# [0] 
	.loc	1	21	0
	ldih	$gp,0($26)               	!gpdisp!4	# [0] 0
	ldi	$gp,0+8($gp)              	!gpdisp!4	# [1] 0
	.loc	1	19	0
	addw	$0,1,$0                  	# [3] 
	stw	$0,32($sp)                	# [4] anon32
	cmpeq	$0,$1,$0                	# [4] 
	beq	$0,.Lt_43_2050            	# [5] 
.BB7_main:
.BB16_main:
	rtc $0
.BB8_main:
	.loc	1	26	0
#  22  		}
#  23  		long t1 = rpcc_();
#  24  
#  25  		double dt = (t1-t0)/1.45e9/1000;
#  26  		printf("%e sec per fork-join\n");
	ldl	$16,.rodata($gp)          	!literal	# [0] .rodata
	.globl	__nldbl_printf
	ldl	$27,__nldbl_printf($gp)   	!literal	# [0] __nldbl_printf
	ldi	$16,24($16)               	# [3] 
	call	$26,($27),__nldbl_printf 	# [3] __nldbl_printf
.BB9_main:
	ldih	$gp,0($26)               	!gpdisp!5	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!5	# [1] 0
	.loc	1	28	0
#  27  	}
#  28  	athread_halt();
	.globl	athread_halt
	ldl	$27,athread_halt($gp)     	!literal	# [2] athread_halt
	call	$26,($27),athread_halt   	# [5] athread_halt
.BB10_main:
	ldih	$gp,0($26)               	!gpdisp!6	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!6	# [1] 0
	rtc $0
	.loc	1	8	0
	stl	$0,40($sp)                	# [3] _temp_gra_spill2
.BB11_main:
	.loc	1	31	0
#  29  
#  30  	long t0 = rpcc_();
#  31  	sleep(1);
	mov	1,$16                     	# [3] 
	.globl	sleep
	ldl	$27,sleep($gp)            	!literal	# [3] sleep
	call	$26,($27),sleep          	# [3] sleep
.BB12_main:
	ldih	$gp,0($26)               	!gpdisp!7	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!7	# [1] 0
	rtc $0
.BB13_main:
	.loc	1	34	0
#  32  	long t1 = rpcc_();
#  33  
#  34  	printf("clock %ld\n", t1-t0);
	ldl	$16,.rodata($gp)          	!literal	# [0] .rodata
	ldl	$17,40($sp)               	# [0] _temp_gra_spill2
	ldl	$27,__nldbl_printf($gp)   	!literal	# [1] __nldbl_printf
	ldi	$16,48($16)               	# [3] 
	subl	$0,$17,$17               	# [3] 
	call	$26,($27),__nldbl_printf 	# [4] __nldbl_printf
.BB14_main:
	ldih	$gp,0($26)               	!gpdisp!8	# [0] 0
	.loc	1	36	0
#  35  
#  36  	return 0;
	ldl	$26,0($sp)                	# [0] return_address
	mov	$31,$0                    	# [0] 
	ldi	$sp,64($sp)               	# [0] 
	.loc	1	34	0
	ldi	$gp,0($gp)                	!gpdisp!8	# [1] 0
	.loc	1	36	0
	ret	$31,($26),1               	# [3] 
.L_CC_main:
# PU cycle count: 0.000000
	.end	main

	.section .rodata
	.org 0x0
	.align	0
	# offset 0
	.ascii	"\x00\x00\x00\x00"	# float 0.00000
	.org 0x4
	.align	0
	# offset 4
	.ascii	"\x00\x00\x80\x4b"	# float 1.67772e+07
	.org 0x8
	.align	0
	# offset 8
	.ascii	"\x00\x00\x00\x00\x00\x00\x00\x00"	# double 0.00000
	.org 0x10
	.align	0
	# offset 16
	.ascii	"\x00\x00\x00\x00\x00\x00\x40\x43"	# double 9.00720e+15
	.org 0x18
	.align	0
	# offset 24
	.byte	0x25, 0x65, 0x20, 0x73, 0x65, 0x63, 0x20, 0x70	# %e sec p
	.byte	0x65, 0x72, 0x20, 0x66, 0x6f, 0x72, 0x6b, 0x2d	# er fork-
	.byte	0x6a, 0x6f, 0x69, 0x6e, 0xa, 0x0	# join\n\000
	.org 0x30
	.align	0
	# offset 48
	.byte	0x63, 0x6c, 0x6f, 0x63, 0x6b, 0x20, 0x25, 0x6c	# clock %l
	.byte	0x64, 0xa, 0x0	# d\n\000
	.weak	_tdata_local_start
	.weak	_tdata_local_end
	.weak	_tdata_private_start
	.weak	_tdata_private_end
	.weak	_tdata_local_fix_end
	.section .text
	.align	4
	.section .rodata
	.align	3
#	.gpvalue 0

	.section .debug_info, "", "progbits"
	.align	0
	.byte	0x4d, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_abbrev
	.long	0x616d0108, 0x72657473, 0x6f00632e, 0x436e6570
	.long	0x2e352043, 0x2d313234, 0x342d7773, 0x01003539
	.byte	0x00
	.long	.debug_line
	.long	0x6d0e0102, 0x006e6961, 0x92040101, 0x0200c01e
	.quad	.BB1_main
	.quad	.L_CC_main
	.byte	0x00, 0x00

	.section .debug_frame, "", "progbits"
	.align	0

	.section .debug_aranges, "", "progbits"
	.align	0
	.byte	0x2c, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.byte	0x08, 0x00, 0x00, 0x00, 0x00, 0x00
	.quad	.BB1_main
	.quad	.L_CC_main - .BB1_main
	.long	0x00000000, 0x00000000, 0x00000000, 0x00000000

	.section .debug_pubnames, "", "progbits"
	.align	0
	.byte	0x17, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.long	0x00000051, 0x0000002f, 0x6e69616d, 0x00000000
	.byte	0x00

	.section .eh_frame, "a", "progbits"
	.align	0

.LEHCIE:
	.long	.LEHCIE_end - .LEHCIE_begin
.LEHCIE_begin:
	.long 0x0
	.byte	0x01, 0x00, 0x01, 0x78, 0x1a, 0x0c, 0x1e, 0x00
	.align 3
.LEHCIE_end:

	.section .debug_abbrev, "", "progbits"
	.align	0
	.long	0x03011101, 0x13082508, 0x100b420b, 0x02000006
	.long	0x0b3a002e, 0x08030b3b, 0x408b0c3f, 0x360a400c
	.byte	0x0b, 0x11, 0x01, 0x12, 0x01, 0x00, 0x00, 0x00
	.byte	0x00
	.section	.note.GNU-stack,"",@progbits
	.ident	"#SWCC Version 5.421-sw-495 : master.c compiled with : -O2 -WOPT:warn_uninit=on -msimd "

