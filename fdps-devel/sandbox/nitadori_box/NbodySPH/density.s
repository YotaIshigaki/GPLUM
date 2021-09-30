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
