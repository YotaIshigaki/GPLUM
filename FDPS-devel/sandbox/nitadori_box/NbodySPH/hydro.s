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
