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
