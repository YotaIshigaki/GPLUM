#include "MGPoisson.h"
#include "MGFFT3DImpl.h"
using namespace MultigridModule;

// =====================================================================

// ---------------------------------------------------------------------

MGPoisson::MGPoisson() {

	numLevel = 0;

	size = (double*)0;
	numGrid = (int*)0;

	vectU0 = (double*)0;
	vectF0 = (double*)0;

	numLayout = 0;
	layout = (mgp_layout_struct**)0;

	numRecall = 1;

	numPre = 2;

	numPost = 2;

	fftLevel = -1;

#ifdef USE_MPI

	comm = MPI_COMM_NULL;

	myRank = MPI_PROC_NULL;

	procMx = (int*)0;

#endif
	isInit = 0;

	NUMITERATION = 0;
	RESID_L2NORM = 0.0;
	CONVFACT = 0.0;

	cleanTimer();

}

// ---------------------------------------------------------------------

MGPoisson::~MGPoisson() {

	destroyLayout();

#ifdef USE_MPI
	if( procMx != (int*)0 ) {
		delete [] procMx;
		procMx = (int*)0;
	}
#endif

}

// =====================================================================

// ---------------------------------------------------------------------

int MGPoisson::Initialize( int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0, int a_fftLevel ) {

	int iErr = _ELSE_ERR_;

	if( ( iErr = setNumGrid( a_numGrid ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = setSize( a_size ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = setVectU0( a_vectU0 ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = setVectF0( a_vectF0 ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = setNumLevel( numLevel ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = setFftLevel( a_fftLevel ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = createLayout() ) != _NO_ERR_ ) {
	}
	else {
		isInit = 1;
	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::Initialize( int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0 ) {

	return Initialize( a_numGrid, a_size, a_vectU0, a_vectF0, -1 );

}

// ---------------------------------------------------------------------

int MGPoisson::Initialize( int *a_flagBC, int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0, int a_fftLevel ) {

	int iErr = _ELSE_ERR_;

	if( ( iErr = setFlagBC( a_flagBC ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = Initialize( a_numGrid, a_size, a_vectU0, a_vectF0, a_fftLevel ) ) != _NO_ERR_ ) {
	}
	else {
		isInit = 1;
		iErr = _NO_ERR_;
	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::Initialize( int *a_flagBC, int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0 ) {

	return Initialize( a_flagBC, a_numGrid, a_size, a_vectU0, a_vectF0, -1 );

}

// ---------------------------------------------------------------------

int MGPoisson::InitVCycle( int a_numRecall, int a_numPre, int a_numPost ) {

	int iErr = _ELSE_ERR_;

	int iLayout;

	numRecall = a_numRecall;
	numPre = a_numPre;

	if( layout != (mgp_layout_struct**)0 && numLayout > 0 ) {
		for( iLayout = 0; iLayout < numLayout; iLayout++ ) {
			layout[ iLayout ]->mgNumPre = a_numPre;
			layout[ iLayout ]->mgNumPost = a_numPost;
			if( ( layout[ iLayout ]->mgLevel > 1 ) && ( layout[ iLayout ]->mgLevel > fftLevel + 1 ) ) {
				layout[ iLayout ]->mgNumRecall = a_numRecall;
			}
			else {
				layout[ iLayout ]->mgNumRecall = 1;
			}
		}
		isInit = 1;
		iErr = _NO_ERR_;
	}
	else {
		isInit = 0;
	}

	return iErr;
}

// ---------------------------------------------------------------------

#ifdef USE_MPI
int MGPoisson::Initialize( MPI_Comm a_comm, int *a_RankArr, int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0, int a_fftLevel ) {

	int iErr = _ELSE_ERR_;

	if( ( iErr = setComm( a_comm ) ) != _NO_ERR_ ) {
	}
	else if( ( iErr = Initialize( a_RankArr, a_numGrid, a_size, a_vectU0, a_vectF0, a_fftLevel ) ) != _NO_ERR_ ) {
	}
	else {
		isInit = 1;
		iErr = _NO_ERR_;
	}

	return iErr;

}

int MGPoisson::Initialize( MPI_Comm a_comm, int *a_RankArr, int *a_numGrid, double *a_size, double *a_vectU0, double *a_vectF0 ) {

	return Initialize( a_comm, a_RankArr, a_numGrid, a_size, a_vectU0, a_vectF0, -1 );
}
#endif

// =====================================================================

// ---------------------------------------------------------------------

int MGPoisson::SORNIter( int numIter ) {

	int iErr = _ELSE_ERR_;

	int iter;

	double t0;

	double RESID_L2NORM_old = DBL_MAX;
	double CONVFACT_sum = 0.0;

	if( isInit ) {

		cleanTimer();

		renewVector();

		t0 = MPI_Wtime();

		for( iter = 0; iter < numIter; iter++ ) {
			iErr = SORCore( layout[ 0 ] );
			if( iErr != _NO_ERR_ ) {
				break;
			}

#ifdef USE_MPI
			double t1;
			double local_RESID_L2NORM = RESID_L2NORM;

			t1 = MPI_Wtime();
			/*
			MPI_Allreduce( &local_RESID_L2NORM, &RESID_L2NORM, 1,
					MPI_DOUBLE_PRECISION, MPI_SUM, comm );
			*/
			MPI_Reduce( &local_RESID_L2NORM, &RESID_L2NORM, 1,
					MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm );
			timer[_MPI2_TIME_]+=MPI_Wtime() - t1;

			if( myRank == 0 ) {
#endif
				RESID_L2NORM =
					sqrt( layout[ 0 ]->spaceGrid[_x_]*layout[ 0 ]->spaceGrid[_y_]*layout[ 0 ]->spaceGrid[_z_]*
							RESID_L2NORM );
				if( RESID_L2NORM_old > 0.0 ) {
					CONVFACT = RESID_L2NORM/RESID_L2NORM_old;
				}
				CONVFACT_sum += CONVFACT;
				RESID_L2NORM_old = RESID_L2NORM;

#ifdef _TEST_
				printf( "# %3d %20.12g %20.12f\n", iter, RESID_L2NORM, CONVFACT );
#endif

#ifdef USE_MPI
			}  // end if
#endif
		}  // end for

		timer[_TOTAL_TIME_] += MPI_Wtime() - t0;

		NUMITERATION = iter;
		if( NUMITERATION > 1 ) {
			CONVFACT = CONVFACT_sum/(double)( NUMITERATION - 1 );
		}

		returnVector();

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::VCycleNIter( int numIter ) {

	int iErr = _ELSE_ERR_;

	int iter;

	double t0;

	double RESID_L2NORM_old = DBL_MAX;
	double CONVFACT_sum = 0.0;

	if( isInit ) {

		cleanTimer();

		renewVector();

		t0 = MPI_Wtime();

		for( iter = 0; iter < numIter; iter++ ) {
			iErr = VCycleCore( layout[ 0 ] );
			if( iErr != _NO_ERR_ ) {
				break;
			}

#ifdef USE_MPI
			double t1;
			double local_RESID_L2NORM = RESID_L2NORM;

			t1 = MPI_Wtime();
			MPI_Allreduce( &local_RESID_L2NORM, &RESID_L2NORM, 1,
					MPI_DOUBLE_PRECISION, MPI_SUM, comm );
			timer[_MPI2_TIME_]+=MPI_Wtime() - t1;
#endif

			RESID_L2NORM =
				sqrt( layout[ 0 ]->spaceGrid[_x_]*layout[ 0 ]->spaceGrid[_y_]*layout[ 0 ]->spaceGrid[ 0 ]*
						RESID_L2NORM );

			if( RESID_L2NORM_old > 0.0 ) {
				CONVFACT = RESID_L2NORM/RESID_L2NORM_old;
			}
			CONVFACT_sum += CONVFACT;
			RESID_L2NORM_old = RESID_L2NORM;

#ifdef _TEST_
			timer[_TOTAL_TIME_] += MPI_Wtime() - t0;
#ifdef USE_MPI
			if( myRank == 0 ) {
				printf( "# %3d %20.12g %20.12f\n", iter, RESID_L2NORM, CONVFACT );
			}
#else
			printf( "# %3d %20.12g %20.12f\n", iter, RESID_L2NORM, CONVFACT );
#endif
			t0 = MPI_Wtime();
#endif
		}  // end for

		timer[_TOTAL_TIME_] += MPI_Wtime() - t0;

		NUMITERATION = iter;
		if( NUMITERATION > 1 ) {
			CONVFACT = CONVFACT_sum/(double)( NUMITERATION - 1 );
		}

		returnVector();

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::FMG( int a_numRecall ) {

	int iErr = _ELSE_ERR_;

	double t0;

	if( isInit ) {

		cleanTimer();

		renewVector();

		t0 = MPI_Wtime();

		iErr = FMGCore( a_numRecall, layout[ 0 ] );

		timer[_TOTAL_TIME_] += MPI_Wtime() - t0;

#ifdef USE_MPI
		double t1;
		double local_RESID_L2NORM = RESID_L2NORM;

		t1 = MPI_Wtime();
		MPI_Allreduce( &local_RESID_L2NORM, &RESID_L2NORM, 1,
				MPI_DOUBLE_PRECISION, MPI_SUM, comm );
		timer[_MPI2_TIME_]+=MPI_Wtime() - t1;
#endif

		RESID_L2NORM =
			sqrt( layout[ 0 ]->spaceGrid[_x_]*layout[ 0 ]->spaceGrid[_y_]*layout[ 0 ]->spaceGrid[ 0 ]*
					RESID_L2NORM );

#ifdef _TEST_FMG_
#ifdef USE_MPI
		if( myRank == 0 ) {
			printf( "# Final Residual L2 Norm = %20.12g\n", RESID_L2NORM );
		}
#else
		printf( "# Final Residual L2 Norm = %20.12g\n", RESID_L2NORM );
#endif
#endif

		NUMITERATION = 1;

		returnVector();

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::SORCore( mgp_layout_struct *a_layout ) {

	int idir;
	int ix, iy, iz;
	int izS;
	int icolor, ncolor;
	long index;
	int idX, idY, idZ;
	double hx2Inv, hy2Inv, hz2Inv;
	double coef, coefInv;
	double resid;

	double t0;

	t0 = MPI_Wtime();

	idX = a_layout->deltIndex[_x_];
	idY = a_layout->deltIndex[_y_];
	idZ = a_layout->deltIndex[_z_];

	hx2Inv = a_layout->mxCoef[_x_];
	hy2Inv = a_layout->mxCoef[_y_];
	hz2Inv = a_layout->mxCoef[_z_];
	coef = a_layout->mxCoef[ 3 ];
	coefInv = 1.0/coef;

	timer[_ELSE_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

		for( idir = 0; idir < 2*_3D_; idir++ ) {
			transBC( idir, a_layout, _VECT_U_ );
		}

	timer[_MPI1_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

	RESID_L2NORM = 0.0;
	ncolor = a_layout->ncolor;
	for( icolor = 0; icolor < ncolor; icolor++ ) {
		for( ix = 1; ix < a_layout->numGrid[_x_] - 1; ix++ ) {
			for( iy = 1; iy < a_layout->numGrid[_y_] - 1; iy++ ) {
				izS = 1 + ((ix+iy+icolor)%ncolor);
				for( iz = izS; iz < a_layout->numGrid[_z_] - 1; iz+=ncolor ) {

					index = getIndex( ix, iy, iz, a_layout->numGrid );

					resid =
						a_layout->vectF[ index ]
						+ hx2Inv*( a_layout->vectU[ index - idX ] - a_layout->vectU[ index ]
								+ a_layout->vectU[ index + idX ] - a_layout->vectU[ index ] )
						+ hy2Inv*( a_layout->vectU[ index - idY ] - a_layout->vectU[ index ]
								+ a_layout->vectU[ index + idY ] - a_layout->vectU[ index ] )
						+ hz2Inv*( a_layout->vectU[ index - idZ ] - a_layout->vectU[ index ]
								+ a_layout->vectU[ index + idZ ] - a_layout->vectU[ index ] );

					a_layout->vectU[ index ] -= a_layout->omega*coefInv*resid;
					RESID_L2NORM += resid*resid;

				}
			}
		}
	}

	timer[_CALC_TIME_] += MPI_Wtime() - t0;

	return 0;

}

// ---------------------------------------------------------------------

int MGPoisson::VCycleCore( mgp_layout_struct *a_layout ) {
	int iErr = _ELSE_ERR_;

	int itimes;

	if( a_layout->mgLevel == fftLevel ) {
		iErr = FFTPoisson( a_layout );
	}
	else if( fftLevel < 0 && a_layout->mgLevel == 0 ) {
		iErr = SORNIter( 1, a_layout );
	}
	else {

		iErr = SORNIter( a_layout->mgNumPre, a_layout );

		iErr = CalcResidual( a_layout );

		iErr = Restriction( a_layout );

		for( itimes = 0; itimes < a_layout->mgNumRecall; itimes++ ) {
			iErr = VCycleCore( a_layout->coarse );
			if( iErr != _NO_ERR_ ) {
				break;
			}
		}

		iErr = Correction( a_layout );

		iErr = SORNIter( a_layout->mgNumPost, a_layout );

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::FMGCore( int a_numRecall, mgp_layout_struct *a_layout ) {
	int iErr = _ELSE_ERR_;

	if( a_layout->mgLevel > 0 ) {

		iErr = CalcResidual( a_layout );

		iErr = copyVector( a_layout, 1, 2 );

		iErr = Restriction( a_layout );

		iErr = FMGCore( a_numRecall, a_layout->coarse );

		iErr = Correction( a_layout );

	}

	for( int itimes = 0; itimes < a_numRecall; itimes++ ) {
		iErr = VCycleCore( a_layout );
	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::FFTPoisson( mgp_layout_struct *a_layout ) {

	int iErr = _ELSE_ERR_;

	double coef[_3D_];
	double coefSum;
	double invCoefSum;
	double angCnst[_3D_];

	int idim;
	int ix, iy, iz;
	int nnn;

	double t0;

	t0 = MPI_Wtime();

	iErr = transMGtoFFT( a_layout, _VECT_F_ );

	timer[_FFT_MPI_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();
#ifdef USE_MPI
	if( iErr == _NO_ERR_ && myRank == _ORIGINAL_RANK_ ) {
#else
	if( iErr == _NO_ERR_ ) {
#endif

                //a_layout->fftF->printReal();
                a_layout->fftF->forward();

		nnn = 1;
		for( idim = 0; idim < _3D_; idim++ ) {
			nnn *= a_layout->fftNumGrid[ idim ];
			angCnst[ idim ] = 2.0*M_PI/(double)(a_layout->fftNumGrid[ idim ]);
		}

                FFT3D::Complex3D& FF = a_layout->fftF->getComplex();
                FFT3D::Complex3D& UU = a_layout->fftU->getComplex();
		for( ix = 0; ix < a_layout->fftF->getComplexEnd(_x_); ix++ ) {
			coef[_x_] = 2.0*a_layout->mxCoef[_x_]*( 1.0 - cos( angCnst[_x_]*(double)ix ) );

			for( iy = 0; iy < a_layout->fftF->getComplexEnd(_y_); iy++ ) {
				coef[_y_] = 2.0*a_layout->mxCoef[_y_]*( 1.0 - cos( angCnst[_y_]*(double)iy ) );

				for( iz = 0; iz < a_layout->fftF->getComplexEnd(_z_); iz++ ) {
					coef[_z_] = 2.0*a_layout->mxCoef[_z_]*( 1.0 - cos( angCnst[_z_]*(double)iz ) );

					coefSum = coef[_x_] + coef[_y_] + coef[_z_];
					if( coefSum != 0.0 ) {
						invCoefSum = 1.0/coefSum;
                                                UU[ix][iy][iz] = invCoefSum*FF[ix][iy][iz];
					}
					else {
                                                UU[ix][iy][iz] = FF[ix][iy][iz];
					}

				}
			}
		}

                a_layout->fftU->backward();
                a_layout->fftU->scaling(1.0/(double)nnn);
	}

	timer[_FFT_TIME_] += MPI_Wtime() - t0;


	t0 = MPI_Wtime();

	iErr = transFFTtoMG( a_layout, _VECT_U_ );

	timer[_FFT_MPI_TIME_] += MPI_Wtime() - t0;

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::SORNIter( int a_iter, mgp_layout_struct *a_layout ) {
	int iErr = _ELSE_ERR_;

	int iter;

#ifdef _TEST1_
	printf( "Call SORNIter:\t Level = %d - %d\n", a_layout->mgLevel, a_iter );
#endif

	for( iter = 0; iter < a_iter; iter++ ) {
		iErr = SORCore( a_layout );
		if( iErr != _NO_ERR_ ) {
			break;
		}
	}

	return iErr;
}

// ---------------------------------------------------------------------

int MGPoisson::CalcResidual( mgp_layout_struct *a_layout ) {

	int iErr = _ELSE_ERR_;

	int idir;
	int ix, iy, iz;
	long index;
	int idX, idY, idZ;
	double hx2Inv, hy2Inv, hz2Inv;
	double coef;

	double t0;

#ifdef _TEST1_
	printf( "Call CalcResidual:\t Level = %d\n", a_layout->mgLevel );
#endif

	t0 = MPI_Wtime();

	idX = a_layout->deltIndex[_x_];
	idY = a_layout->deltIndex[_y_];
	idZ = a_layout->deltIndex[_z_];

	hx2Inv = a_layout->mxCoef[_x_];
	hy2Inv = a_layout->mxCoef[_y_];
	hz2Inv = a_layout->mxCoef[_z_];
	coef = a_layout->mxCoef[ 3 ];

	timer[_ELSE_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

	for( idir = 0; idir < 2*_3D_; idir++ ) {
		iErr = transBC( idir, a_layout, _VECT_U_ );
	}

	timer[_MPI1_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

	for( ix = 1; ix < a_layout->numGrid[_x_] - 1; ix++ ) {
		for( iy = 1; iy < a_layout->numGrid[_y_] - 1; iy++ ) {
			for( iz = 1; iz < a_layout->numGrid[_z_] - 1; iz++ ) {

				index = getIndex( ix, iy, iz, a_layout->numGrid );

				a_layout->vectR[ index ] =
					a_layout->vectF[ index ]
					+ hx2Inv*( a_layout->vectU[ index - idX ] - a_layout->vectU[ index ]
							+ a_layout->vectU[ index + idX ] - a_layout->vectU[ index ] )
					+ hy2Inv*( a_layout->vectU[ index - idY ] - a_layout->vectU[ index ]
							+ a_layout->vectU[ index + idY ] - a_layout->vectU[ index ] )
					+ hz2Inv*( a_layout->vectU[ index - idZ ] - a_layout->vectU[ index ]
							+ a_layout->vectU[ index + idZ ] - a_layout->vectU[ index ] );

			}
		}
	}

	timer[_CALC_TIME_]+=MPI_Wtime() - t0;

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::Restriction( mgp_layout_struct *a_layout ) {
	int iErr = _ELSE_ERR_;

	int idir;
	long idxF, idxC;
	int idX, idY, idZ;
	int ixF, iyF, izF;
	int ixC, iyC, izC;

	double fact = 1.0/12.0;

	double t0;

#ifdef _TEST1_
	printf( "Call Restriction:\t Level = %d\n", a_layout->mgLevel );
#endif

	t0 = MPI_Wtime();

	idX = a_layout->deltIndex[_x_];
	idY = a_layout->deltIndex[_y_];
	idZ = a_layout->deltIndex[_z_];

	timer[_ELSE_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();
	for( idir = _3D_; idir < 2*_3D_; idir++ ) {
		iErr = transBC( idir, a_layout, _VECT_R_ );
	}
	timer[_MPI1_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

	for( ixF = 1, ixC = 1; ixF < a_layout->numGrid[_x_] - 1; ixF+=2, ixC++ ) {
		for( iyF = 1, iyC = 1; iyF < a_layout->numGrid[_y_] - 1; iyF+=2, iyC++ ) {
			for( izF = 1, izC = 1; izF < a_layout->numGrid[_z_] - 1; izF+=2, izC++ ) {

				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				idxC = getIndex( ixC, iyC, izC, a_layout->coarse->numGrid );

				a_layout->coarse->vectF[ idxC ] = 0.5*a_layout->vectR[ idxF ]
					+ fact*( a_layout->vectR[ idxF - idX ] + a_layout->vectR[ idxF + idX ]
							+ a_layout->vectR[ idxF - idY ] + a_layout->vectR[ idxF + idY ]
							+ a_layout->vectR[ idxF - idZ ] + a_layout->vectR[ idxF + idZ ] );
				a_layout->coarse->vectU[ idxC ] = 0.0;
			}
		}
	}

	timer[_CALC_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();
	for( idir = 0; idir < _3D_; idir++ ) {
		iErr = transBC( idir, a_layout->coarse, _VECT_F_ );
	}
	timer[_MPI1_TIME_] += MPI_Wtime() - t0;

	iErr = _NO_ERR_;

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::Correction( mgp_layout_struct *a_layout ) {
	int iErr = _ELSE_ERR_;

	int idir;
	long idxF, idxC;
	int ixF, iyF, izF;
	int ixC, iyC, izC;
	int idX, idY, idZ;

	double t0;

#ifdef _TEST1_
	printf( "Call Corection:\t Level = %d\n", a_layout->mgLevel );
#endif

	t0 = MPI_Wtime();

	idX = a_layout->deltIndex[_x_];
	idY = a_layout->deltIndex[_y_];
	idZ = a_layout->deltIndex[_z_];

	timer[_ELSE_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();
	for( idir = 0; idir < 2*_3D_; idir++ ) {
		iErr = transBC( idir, a_layout->coarse, _VECT_U_ );
	}
	timer[_MPI1_TIME_] += MPI_Wtime() - t0;

	t0 = MPI_Wtime();

	for( ixC = 1, ixF = 1; ixF < a_layout->numGrid[_x_]; ixC++, ixF+=2 ) {
		for( iyC = 1, iyF = 1; iyF < a_layout->numGrid[_y_]; iyC++, iyF+=2 ) {
			for( izC = 1, izF = 1; izF < a_layout->numGrid[_z_]; izC++, izF+=2 ) {
				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				idxC = getIndex( ixC, iyC, izC, a_layout->coarse->numGrid );
				a_layout->vectR[ idxF ] = a_layout->coarse->vectU[ idxC ];

			}
		}
	}

	for( ixF = 2; ixF < a_layout->numGrid[_x_] - 1; ixF+=2 ) {
		for( iyF = 1; iyF < a_layout->numGrid[_y_]; iyF+=2 ) {
			for( izF = 1; izF < a_layout->numGrid[_z_]; izF+=2 ) {
				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				a_layout->vectR[ idxF ] =
					0.5*( a_layout->vectR[ idxF + idX ] + a_layout->vectR[ idxF - idX ] );
			}
		}
	}

	for( ixF = 1; ixF < a_layout->numGrid[_x_]; ixF++ ) {
		for( iyF = 2; iyF < a_layout->numGrid[_y_] - 1; iyF+=2 ) {
			for( izF = 1; izF < a_layout->numGrid[_z_]; izF+=2 ) {
				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				a_layout->vectR[ idxF ] =
					0.5*( a_layout->vectR[ idxF + idY ] + a_layout->vectR[ idxF - idY ] );
			}
		}
	}

	for( ixF = 1; ixF < a_layout->numGrid[_x_]; ixF++ ) {
		for( iyF = 1; iyF < a_layout->numGrid[_y_]; iyF++ ) {
			for( izF = 2; izF < a_layout->numGrid[_z_] - 1; izF+=2 ) {
				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				a_layout->vectR[ idxF ] =
					0.5*( a_layout->vectR[ idxF + idZ ] + a_layout->vectR[ idxF - idZ ] );
			}
		}
	}

	for( ixF = 1; ixF < a_layout->numGrid[_x_]; ixF++ ) {
		for( iyF = 1; iyF < a_layout->numGrid[_y_]; iyF++ ) {
			for( izF = 1; izF < a_layout->numGrid[_z_]; izF++ ) {
				idxF = getIndex( ixF, iyF, izF, a_layout->numGrid );
				a_layout->vectU[ idxF ] += a_layout->vectR[ idxF ];
			}
		}
	}

	timer[_CALC_TIME_] += MPI_Wtime() - t0;

	return _NO_ERR_;

}

// ---------------------------------------------------------------------

int MGPoisson::transBC( int a_dir, mgp_layout_struct *a_layout, int a_flgVect ) {

	int iErr = _ELSE_ERR_;
	int ix, iy, iz;
	int idxSend, idxRecv;

	double *transVect;

	if( a_flgVect == _VECT_U_ ) {
		transVect = a_layout->vectU;
	}
	else if( a_flgVect == _VECT_F_ ) {
		transVect = a_layout->vectF;
	}
	else if( a_flgVect == _VECT_R_ ) {
		transVect = a_layout->vectR;
	}
	else {
		return _ILLEGAL_VALUE_;
	}

#ifdef USE_MPI
	int mpi_err_send, mpi_err_recv;
	int sendRank, recvRank;
	MPI_Request sendReq, recvReq;
	MPI_Status mpi_statusSend, mpi_statusRecv;

	ix = 0*(a_dir==_x_) + 2*(a_dir==(_x_+_3D_)) + ((a_dir%_3D_)!=_x_);
	iy = 0*(a_dir==_y_) + 2*(a_dir==(_y_+_3D_)) + ((a_dir%_3D_)!=_y_);
	iz = 0*(a_dir==_z_) + 2*(a_dir==(_z_+_3D_)) + ((a_dir%_3D_)!=_z_);
	idxSend = getIndex( ix, iy, iz, _3D_, _3D_, _3D_ );
	sendRank = a_layout->flagBC[ idxSend ];

	ix = 2*(a_dir==_x_) + 0*(a_dir==(_x_+_3D_)) + ((a_dir%_3D_)!=_x_);
	iy = 2*(a_dir==_y_) + 0*(a_dir==(_y_+_3D_)) + ((a_dir%_3D_)!=_y_);
	iz = 2*(a_dir==_z_) + 0*(a_dir==(_z_+_3D_)) + ((a_dir%_3D_)!=_z_);
	idxRecv = getIndex( ix, iy, iz, _3D_, _3D_, _3D_ );
	recvRank = a_layout->flagBC[ idxRecv ];

	mpi_err_recv = MPI_Irecv( transVect, 1, a_layout->recvType[ a_dir ], recvRank, a_dir, comm, &recvReq );

	mpi_err_send = MPI_Isend( transVect, 1, a_layout->sendType[ a_dir ], sendRank, a_dir, comm, &sendReq );

	if( mpi_err_recv == MPI_SUCCESS && mpi_err_send == MPI_SUCCESS ) {

		mpi_err_recv = MPI_Wait( &recvReq, &mpi_statusRecv );

		mpi_err_send = MPI_Wait( &sendReq, &mpi_statusSend );
	}

	if( mpi_err_recv != MPI_SUCCESS ) {
		iErr = mpi_err_recv;
	}
	else if( mpi_err_send != MPI_SUCCESS ) {
		iErr = mpi_err_send;
	}
	else {
		iErr = _NO_ERR_;
	}

#else
	int ixSend, iySend, izSend;
	int ixRecv, iyRecv, izRecv;

	if( a_dir == 0 || a_dir == 3 ) {
		for( ix = 0; ix < a_layout->transSize[ a_dir ][_x_]; ix++ ) {
			ixSend = a_layout->sendStart[ a_dir ][_x_] + ix;
			ixRecv = a_layout->recvStart[ a_dir ][_x_] + ix;
			for( iy = 0; iy < a_layout->transSize[ a_dir ][_y_]; iy++ ) {
				iySend = a_layout->sendStart[ a_dir ][_y_] + iy;
				iyRecv = a_layout->recvStart[ a_dir ][_y_] + iy;
				for( iz = 0; iz < a_layout->transSize[ a_dir ][_z_]; iz++ ) {
					izSend = a_layout->sendStart[ a_dir ][_z_] + iz;
					izRecv = a_layout->recvStart[ a_dir ][_z_] + iz;
					idxSend = getIndex( ixSend, iySend, izSend, a_layout->numGrid );
					idxRecv = getIndex( ixRecv, iyRecv, izRecv, a_layout->numGrid );

					transVect[ idxRecv ] = transVect[ idxSend ];
				}
			}
		}
		iErr = _NO_ERR_;
	}
	else if( a_dir == 1 || a_dir == 4 ) {
		for( iy = 0; iy < a_layout->transSize[ a_dir ][_y_]; iy++ ) {
			iySend = a_layout->sendStart[ a_dir ][_y_] + iy;
			iyRecv = a_layout->recvStart[ a_dir ][_y_] + iy;
			for( ix = 0; ix < a_layout->transSize[ a_dir ][_x_]; ix++ ) {
				ixSend = a_layout->sendStart[ a_dir ][_x_] + ix;
				ixRecv = a_layout->recvStart[ a_dir ][_x_] + ix;
				for( iz = 0; iz < a_layout->transSize[ a_dir ][_z_]; iz++ ) {
					izSend = a_layout->sendStart[ a_dir ][_z_] + iz;
					izRecv = a_layout->recvStart[ a_dir ][_z_] + iz;
					idxSend = getIndex( ixSend, iySend, izSend, a_layout->numGrid );
					idxRecv = getIndex( ixRecv, iyRecv, izRecv, a_layout->numGrid );

					transVect[ idxRecv ] = transVect[ idxSend ];
				}
			}
		}
		iErr = _NO_ERR_;
	}
	else if( a_dir == 2 || a_dir == 5 ) {
		for( iz = 0; iz < a_layout->transSize[ a_dir ][_z_]; iz++ ) {
			izSend = a_layout->sendStart[ a_dir ][_z_] + iz;
			izRecv = a_layout->recvStart[ a_dir ][_z_] + iz;
			for( ix = 0; ix < a_layout->transSize[ a_dir ][_x_]; ix++ ) {
				ixSend = a_layout->sendStart[ a_dir ][_x_] + ix;
				ixRecv = a_layout->recvStart[ a_dir ][_x_] + ix;
				for( iy = 0; iy < a_layout->transSize[ a_dir ][_y_]; iy++ ) {
					iySend = a_layout->sendStart[ a_dir ][_y_] + iy;
					iyRecv = a_layout->recvStart[ a_dir ][_y_] + iy;
					idxSend = getIndex( ixSend, iySend, izSend, a_layout->numGrid );
					idxRecv = getIndex( ixRecv, iyRecv, izRecv, a_layout->numGrid );

					transVect[ idxRecv ] = transVect[ idxSend ];
				}
			}
		}
		iErr = _NO_ERR_;
	}
	else {
		iErr = _ILLEGAL_VALUE_;
	}
#endif

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::transMGtoFFT( mgp_layout_struct *a_layout, int a_flgVect ) {

	int iErr = _ELSE_ERR_;

	double *sendVect;
	double *recvVect;

	if( a_flgVect == _VECT_U_ ) {
		sendVect = a_layout->vectU;
		recvVect = a_layout->fftBuff;
	}
	else if( a_flgVect == _VECT_F_ ) {
		sendVect = a_layout->vectF;
		recvVect = a_layout->fftBuff;
	}
	else {
		return _ILLEGAL_VALUE_;
	}

#ifdef USE_MPI
	int mpi_err_send, mpi_err_recv;
	int mpi_np = numProc[_x_]*numProc[_y_]*numProc[_z_];
	MPI_Request sendReq, *recvReq;
	MPI_Status sendStatus, recvStatus;

	int idx;

	recvReq = new(std::nothrow) MPI_Request [ mpi_np ];
	if( recvReq == NULL ) {
		return _ALLOC_ERR_;
	}

	// Recv
	mpi_err_recv = MPI_SUCCESS;
	if( myRank == _ORIGINAL_RANK_ ) {
		for( idx = 0; idx < mpi_np; idx++ ) {
			recvReq[ idx ] = 0;
			mpi_err_recv = MPI_Irecv( recvVect, 1, a_layout->globalGridDataType[ idx ], procMx[ idx ],
					1, comm, &recvReq[ idx ] );

			if( mpi_err_recv != MPI_SUCCESS ) {
				break;
			}
		}
	}

	// Send
	mpi_err_send = MPI_Isend( sendVect, 1, a_layout->localGridDataType, _ORIGINAL_RANK_,
			1, comm, &sendReq );

	if( mpi_err_recv == MPI_SUCCESS && mpi_err_send == MPI_SUCCESS ) {

		// Wait Recv
		if( myRank == _ORIGINAL_RANK_ ) {
			for( idx = 0; idx < mpi_np; idx++ ) {
				mpi_err_recv = MPI_Wait( &recvReq[ idx ], &recvStatus );
				if( mpi_err_recv != MPI_SUCCESS ) {
					break;
				}
			}
		}

		// Wait Send
		mpi_err_send = MPI_Wait( &sendReq, &sendStatus );

	}

	if( mpi_err_recv != MPI_SUCCESS ) {
		iErr = mpi_err_recv;
	}
	else if( mpi_err_send != MPI_SUCCESS ) {
		iErr = mpi_err_send;
	}
	else {
		iErr = _NO_ERR_;
	}

	if( recvReq != NULL ) {
		delete [] recvReq;
		recvReq = NULL;
	}

#else

	int ix, iy, iz;
	int sendIdx, recvIdx;

	for( ix = 1; ix < a_layout->numGrid[_x_] - 1; ix++ ) {
		for( iy = 1; iy < a_layout->numGrid[_y_] - 1; iy++ ) {
			for( iz = 1; iz < a_layout->numGrid[_z_] - 1; iz++ ) {
				sendIdx = getIndex( ix, iy, iz, a_layout->numGrid );
				recvIdx = getIndex( ix - 1, iy - 1, iz - 1, a_layout->fftNumGrid );
				recvVect[ recvIdx ] = sendVect[ sendIdx ];
			}
		}
	}

	iErr = _NO_ERR_;

#endif
	if( a_flgVect == _VECT_U_ ) {
		copyToFFT3D(a_layout->fftU, a_layout->fftBuff);
        }
	else if( a_flgVect == _VECT_F_ ) {
		copyToFFT3D(a_layout->fftF, a_layout->fftBuff);
        }

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::transFFTtoMG( mgp_layout_struct *a_layout, int a_flgVect ) {
	if( a_flgVect == _VECT_U_ ) {
		copyFromFFT3D(a_layout->fftU, a_layout->fftBuff);
        }
	else if( a_flgVect == _VECT_F_ ) {
		copyFromFFT3D(a_layout->fftF, a_layout->fftBuff);
        }

	int iErr = _ELSE_ERR_;

	double *sendVect;
	double *recvVect;

	if( a_flgVect == _VECT_U_ ) {
		sendVect = a_layout->fftBuff;
		recvVect = a_layout->vectU;
	}
	else if( a_flgVect == _VECT_F_ ) {
		sendVect = a_layout->fftBuff;
		recvVect = a_layout->vectF;
	}
	else {
		return _ILLEGAL_VALUE_;
	}

#ifdef USE_MPI
	int mpi_err_send, mpi_err_recv;
	int mpi_np = numProc[_x_]*numProc[_y_]*numProc[_z_];
	MPI_Request *sendReq, recvReq;
	MPI_Status sendStatus, recvStatus;

	int idx;

	sendReq = new(std::nothrow) MPI_Request [ mpi_np ];
	if( sendReq == NULL ) {
		return _ALLOC_ERR_;
	}

	// Recv
	mpi_err_recv = MPI_Irecv( recvVect, 1, a_layout->localGridDataType, _ORIGINAL_RANK_,
			1, comm, &recvReq );

	// Send
	mpi_err_send = MPI_SUCCESS;
	if( myRank == _ORIGINAL_RANK_ ) {
		for( idx = 0; idx < mpi_np; idx++ ) {
			sendReq[ idx ] = 0;
			mpi_err_send = MPI_Isend( sendVect, 1, a_layout->globalGridDataType[ idx ], procMx[ idx ],
					1, comm, &sendReq[ idx ] );
			if( mpi_err_send != MPI_SUCCESS ) {
				break;
			}
		}
	}

	if( mpi_err_recv == MPI_SUCCESS && mpi_err_send == MPI_SUCCESS ) {

		// Wait Recv
		mpi_err_recv = MPI_Wait( &recvReq, &recvStatus );

		// Wait Send
		if( myRank == _ORIGINAL_RANK_ ) {
			for( idx = 0; idx < mpi_np; idx++ ) {
				mpi_err_send = MPI_Wait( &sendReq[ idx ], &sendStatus );
				if( mpi_err_send != MPI_SUCCESS ) {
					break;
				}
			}
		}

	}

	if( mpi_err_recv != MPI_SUCCESS ) {
		iErr = mpi_err_recv;
	}
	else if( mpi_err_send != MPI_SUCCESS ) {
		iErr = mpi_err_send;
	}
	else {
		iErr = _NO_ERR_;
	}

	if( sendReq != NULL ) {
		delete [] sendReq;
		sendReq = NULL;
	}

#else

	int ix, iy, iz;
	int sendIdx, recvIdx;

	for( ix = 1; ix < a_layout->numGrid[_x_] - 1; ix++ ) {
		for( iy = 1; iy < a_layout->numGrid[_y_] - 1; iy++ ) {
			for( iz = 1; iz < a_layout->numGrid[_z_] - 1; iz++ ) {
				sendIdx = getIndex( ix - 1, iy - 1, iz - 1, a_layout->fftNumGrid );
				recvIdx = getIndex( ix, iy, iz, a_layout->numGrid );
				recvVect[ recvIdx ] = sendVect[ sendIdx ];
			}
		}
	}

	iErr = _NO_ERR_;

#endif

	return iErr;

}


// =====================================================================

// ---------------------------------------------------------------------

int MGPoisson::createLayout() {

	int iErr = _ELSE_ERR_;
	int tmp_numGrid[_3D_];
	double tmp_size[_3D_];

	int iLevel;
	int iLayout;
	int idim;
	int ndiv;

	if( numLayout < 1 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		layout = new(std::nothrow) mgp_layout_struct * [ numLayout ];
		if( layout == (mgp_layout_struct**)0 ) {
			iErr = _ALLOC_ERR_;
		}
		else {

			iLevel = numLevel - 1;
			ndiv = 1;

			for( iLayout = 0; iLayout < numLayout; iLayout++ ) {

				for( idim = 0; idim < _3D_; idim++ ) {
					tmp_numGrid[ idim ] = numGrid[ idim ]/ndiv;
					tmp_size[ idim ] = size[ idim ];
				}

				layout[ iLayout ] = new(std::nothrow) mgp_layout_struct;
				if( layout[ iLayout ] == (mgp_layout_struct*)0 ) {
					iErr = _ALLOC_ERR_;
					break;
				}
				else if( ( iErr = allocLayout( layout[ iLayout ],
								iLevel, tmp_numGrid, tmp_size ) ) != _NO_ERR_ ) {
					break;
				}

				iLevel--;
				ndiv *= 2;

			}  // end for iLayout

			if( iErr == _NO_ERR_ ) {
				for( iLayout = 0; iLayout < numLayout; iLayout++ ) {
					if( iLayout == 0 ) {
						layout[ iLayout ]->fine = (mgp_layout_struct*)0;
						layout[ iLayout ]->coarse = layout[ iLayout + 1 ];
					}
					else if( iLayout == numLayout - 1 ) {
						layout[ iLayout ]->fine = layout[ iLayout - 1 ];
						layout[ iLayout ]->coarse = (mgp_layout_struct*)0;
					}
					else {
						layout[ iLayout ]->fine = layout[ iLayout - 1 ];
						layout[ iLayout ]->coarse = layout[ iLayout + 1 ];
					}
				}
			}  // end if iErr

		}  // end if layout

	}

	return iErr;

}

// ---------------------------------------------------------------------

void MGPoisson::destroyLayout() {

	int iLayout;

	if( layout != (mgp_layout_struct**)0 ) {

		for( iLayout = 0; iLayout < numLayout; iLayout++ ) {

			if( layout[ iLayout ] != (mgp_layout_struct*)0 ) {
				freeLayout( layout[ iLayout ] );
				delete layout[ iLayout ];
				layout[ iLayout ] = (mgp_layout_struct*)0;
			}

		}

		delete [] layout;
		layout = (mgp_layout_struct**)0;

	}

}

// ---------------------------------------------------------------------

void MGPoisson::renewVector() {

	long idx0, idx1;
	int ix, iy, iz;

	for( ix = 0; ix < numGrid[_x_]; ix++ ) {
		for( iy = 0; iy < numGrid[_y_]; iy++ ) {
			for( iz = 0; iz < numGrid[_z_]; iz++ ) {
				idx0 = getIndex( ix, iy, iz, numGrid );
				idx1 = getIndex( ix+1, iy+1, iz+1, layout[ 0 ]->numGrid );
				layout[ 0 ]->vectU[ idx1 ] = vectU0[ idx0 ];
				layout[ 0 ]->vectF[ idx1 ] = vectF0[ idx0 ];
			}
		}
	}

}

// ---------------------------------------------------------------------

void MGPoisson::returnVector() {

	long idx0, idx1;
	int ix, iy, iz;

	for( ix = 0; ix < numGrid[_x_]; ix++ ) {
		for( iy = 0; iy < numGrid[_y_]; iy++ ) {
			for( iz = 0; iz < numGrid[_z_]; iz++ ) {
				idx0 = getIndex( ix, iy, iz, numGrid );
				idx1 = getIndex( ix+1, iy+1, iz+1, layout[ 0 ]->numGrid );
				vectU0[ idx0 ] = layout[ 0 ]->vectU[ idx1 ];
			}
		}
	}

}

// ---------------------------------------------------------------------
int MGPoisson::copyVector( mgp_layout_struct *a_layout, int a_flgVectFrom, int a_flgVectTo ) {

	int iErr = _ELSE_ERR_;

	long iElem;
	int ix, iy, iz;

	double *vectFrom, *vectTo;

	if( a_flgVectFrom != a_flgVectTo && a_flgVectFrom < 3 && a_flgVectTo < 3) {

		if( a_flgVectFrom == _VECT_U_ ) {
			vectFrom = a_layout->vectU;
		}
		else if( a_flgVectFrom == _VECT_F_ ) {
			vectFrom = a_layout->vectF;
		}
		else if( a_flgVectFrom == _VECT_R_ ) {
			vectFrom = a_layout->vectR;
		}
		else {
			return _ILLEGAL_VALUE_;
		}

		if( a_flgVectTo == _VECT_U_ ) {
			vectTo = a_layout->vectU;
		}
		else if( a_flgVectTo == _VECT_F_ ) {
			vectTo = a_layout->vectF;
		}
		else if( a_flgVectTo == _VECT_R_ ) {
			vectTo = a_layout->vectR;
		}
		else {
			return _ILLEGAL_VALUE_;
		}

		for( ix = 0; ix < a_layout->numGrid[_x_]; ix++ ) {
			for( iy = 0; iy < a_layout->numGrid[_y_]; iy++ ) {
				for( iz = 0; iz < a_layout->numGrid[_z_]; iz++ ) {
					iElem = getIndex( ix, iy, iz, a_layout->numGrid );
					vectTo[ iElem ] = vectFrom[ iElem ];
				}
			}
		}
		/*
		for( iElem = 0; iElem < a_layout->numElem; iElem++ ) {
			vectTo[ iElem ] = vectFrom[ iElem ];
		}
		*/

		iErr = _NO_ERR_;

	}
	else {
		iErr = _ILLEGAL_VALUE_;
	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::allocLayout( mgp_layout_struct *a_layout,
		int a_level, int *a_numGrid, double *a_size ) {

	int iErr = _ELSE_ERR_;
	int deltaIndex;

	if( a_layout == (mgp_layout_struct*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( a_level < 0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( a_numGrid == (int*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( a_size == (double*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		a_layout->numElem = 1;
		for( int idim = 0; idim < _3D_; idim++ ) {
			a_layout->numGrid[ idim ] = a_numGrid[ idim ] + 2;
			a_layout->size[ idim ] = a_size[ idim ];
			a_layout->spaceGrid[ idim ] = a_layout->size[ idim ]/(double)a_numGrid[ idim ];
			a_layout->numElem *= a_layout->numGrid[ idim ];
		}

		for( int ibc = 0; ibc < _numBC_; ibc++ ) {
			a_layout->flagBC[ ibc ] = flagBC[ ibc ];
		}

		deltaIndex = a_layout->numElem;
		a_layout->mxCoef[_3D_] = 0.0;
		for( int idim = 0; idim < _3D_; idim++ ) {
			deltaIndex /= a_layout->numGrid[ idim ];

			a_layout->deltIndex[ idim ] = -deltaIndex;
			a_layout->deltIndex[ 2*_3D_ - idim ] = deltaIndex;

			a_layout->mxCoef[ idim ] = 1.0/( a_layout->spaceGrid[ idim ]*a_layout->spaceGrid[ idim ] );
			a_layout->mxCoef[ 2*_3D_ - idim ] = a_layout->mxCoef[ idim ];
			a_layout->mxCoef[_3D_] += a_layout->mxCoef[ idim ];
		}
		a_layout->deltIndex[_3D_] = 0;
		a_layout->mxCoef[_3D_] *= -2.0;

		if( a_level != 0 ) {
			a_layout->ncolor = 2;
		}
		else {
			a_layout->ncolor = 1;
		}
		a_layout->omega = 1.0;
		a_layout->mgLevel = a_level;
		a_layout->mgNumPre = numPre;
		a_layout->mgNumPost = numPost;

		if( a_layout->mgLevel > 1 && a_layout->mgLevel > fftLevel + 1 ) {
			a_layout->mgNumRecall = numRecall;
		}
		else {
			a_layout->mgNumRecall = 1;
		}

		a_layout->vectF = new(std::nothrow) double [ a_layout->numElem ];
		a_layout->vectU = new(std::nothrow) double [ a_layout->numElem ];
		a_layout->vectR = new(std::nothrow) double [ a_layout->numElem ];
		if( a_layout->vectF == (double*)0 || a_layout->vectU == (double*)0 ||
				a_layout->vectR == (double*)0 ) {
			return  _ALLOC_ERR_;
		}
		else {
			for( long ielem = 0; ielem < a_layout->numElem; ielem++ ) {
				a_layout->vectF[ ielem ] = 0.0;
				a_layout->vectU[ ielem ] = 0.0;
				a_layout->vectR[ ielem ] = 0.0;
			}
			iErr = _NO_ERR_;
		}

		for( int idir = 0; idir < 2*_3D_; idir++ ) {
			for( int idim = 0; idim < _3D_; idim++ ) {

				a_layout->transSize[ idir ][ idim ] =
					((idir%_3D_)==idim) + ((idir%_3D_)!=idim)*( a_layout->numGrid[ idim ] - 1 );

				a_layout->sendStart[ idir ][ idim ] =
					(idir==idim) + (idir==(idim+_3D_))*( a_layout->numGrid[ idim ] - 2 ) + ((idir%_3D_)!=idim);

				a_layout->recvStart[ idir ][ idim ] =
					(idir==idim)*( a_layout->numGrid[ idim ] - 1 ) + (idir==(idim+_3D_))*0 + ((idir%_3D_)!=idim);

			}
		}

#ifdef USE_MPI

		int mpi_err;

		for( int idir = 0; idir < 2*_3D_; idir++ ) {

			if( ( mpi_err = MPI_Type_create_subarray( _3D_, a_layout->numGrid,
					a_layout->transSize[ idir ], a_layout->sendStart[ idir ],
					MPI_ORDER_C, MPI_DOUBLE_PRECISION, &(a_layout->sendType[ idir ]) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}
			else if( ( mpi_err = MPI_Type_commit( &(a_layout->sendType[ idir ] ) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}

			else if( ( mpi_err = MPI_Type_create_subarray( _3D_, a_layout->numGrid,
					a_layout->transSize[ idir ], a_layout->recvStart[ idir ],
					MPI_ORDER_C, MPI_DOUBLE_PRECISION, &(a_layout->recvType[ idir ]) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}
			else if( ( mpi_err = MPI_Type_commit( &(a_layout->recvType[ idir ] ) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}
			else {
				iErr = _NO_ERR_;
			}

		}  // end for idir

#endif

		if( iErr == _NO_ERR_ && a_layout->mgLevel == fftLevel ) {
			iErr = allocFFTLayout( a_layout, a_numGrid );
		}

	}  // end if INPUT CHECK

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::allocFFTLayout( mgp_layout_struct *a_layout, int *a_numGrid ) {

	int iErr = _ELSE_ERR_;

//	int fftNumGrid[_3D_];
	int fftNumElem;
	int rfftNumElem;
	int idim;
	int idx;

#ifdef USE_MPI
	int mpi_err;
	int totNumProc;
	int fftGridDisp[_3D_];
	int irank;

#else

#endif


#ifdef USE_MPI

	fftNumElem = 1;
	totNumProc = 1;
	for( idim = 0; idim < _3D_; idim++ ) {
		a_layout->fftNumGrid[ idim ] = a_numGrid[ idim ]*numProc[ idim ];
		a_layout->globalGridDisp[ idim ] = a_numGrid[ idim ]*idxProc[ idim ];
		fftNumElem *= a_layout->fftNumGrid[ idim ];
		totNumProc *= numProc[ idim ];
	}
	rfftNumElem = (a_layout->fftNumGrid[_z_]/2 + 1 )*fftNumElem/a_layout->fftNumGrid[_z_];

#else

	fftNumElem = 1;
	for( idim = 0; idim < _3D_; idim++ ) {
		a_layout->fftNumGrid[ idim ] = a_numGrid[ idim ];
		fftNumElem *= a_layout->fftNumGrid[ idim ];
	}
	rfftNumElem = (a_layout->fftNumGrid[_z_]/2 + 1 )*fftNumElem/a_layout->fftNumGrid[_z_];

#endif

        a_layout->fftU = FFT3D::createFFT3D(a_layout->fftNumGrid[_x_],
                                            a_layout->fftNumGrid[_y_],
                                            a_layout->fftNumGrid[_z_], true);
        a_layout->fftF = FFT3D::createFFT3D(a_layout->fftNumGrid[_x_],
                                            a_layout->fftNumGrid[_y_],
                                            a_layout->fftNumGrid[_z_], true);
        a_layout->fftBuff = new double[ fftNumElem ];
	if( a_layout->fftU == 0 || a_layout->fftF == 0 || a_layout->fftBuff == 0) {
		return  _ALLOC_ERR_;
	}
	else {
		iErr = _NO_ERR_;
	}

#ifdef USE_MPI
	a_layout->globalGridDataType = new(std::nothrow) MPI_Datatype [ totNumProc ];
	if( a_layout->globalGridDataType == (MPI_Datatype*)0 ) {
		return _ALLOC_ERR_;
	}
	else {
		iErr = _NO_ERR_;
	}

	fftGridDisp[_x_] = 1;
	fftGridDisp[_y_] = 1;
	fftGridDisp[_z_] = 1;

	if( ( mpi_err = MPI_Type_create_subarray( _3D_, a_layout->numGrid, a_numGrid, fftGridDisp,
			MPI_ORDER_C, MPI_DOUBLE_PRECISION, &(a_layout->localGridDataType) ) ) != MPI_SUCCESS ) {
		iErr = mpi_err;
	}
	else if( ( mpi_err = MPI_Type_commit( &(a_layout->localGridDataType) ) ) != MPI_SUCCESS ) {
		iErr = mpi_err;
	}
	else {

		for( idx = 0; idx < totNumProc; idx++ ) {

			irank = procMx[ idx ];
			for( idim = _3D_ - 1; idim >= 0; idim-- ) {
				fftGridDisp[ idim ] = ( irank%numProc[ idim ] )*a_numGrid[ idim ];
				irank /= numProc[ idim ];
			}

			irank = procMx[ idx ];
			if( ( mpi_err = MPI_Type_create_subarray( _3D_, a_layout->fftNumGrid, a_numGrid, fftGridDisp,
					MPI_ORDER_C, MPI_DOUBLE_PRECISION, &(a_layout->globalGridDataType[ irank ]) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}
			else if( ( mpi_err = MPI_Type_commit( &(a_layout->globalGridDataType[ irank ]) ) )
					!= MPI_SUCCESS ) {
				iErr = mpi_err;
				break;
			}
			else {
				iErr = _NO_ERR_;
			}

		}  // end for idx

	}

#endif

	return iErr;

}

// ---------------------------------------------------------------------

void MGPoisson::freeLayout( mgp_layout_struct *a_layout ) {

	a_layout->fine = a_layout;
	a_layout->coarse = a_layout;

	if( a_layout->vectR != (double*)0 ) {
		delete [] a_layout->vectR;
		a_layout->vectR = (double*)0;
	}
	if( a_layout->vectF != (double*)0 ) {
		delete [] a_layout->vectF;
		a_layout->vectF = (double*)0;
	}
	if( a_layout->vectU != (double*)0 ) {
		delete [] a_layout->vectU;
		a_layout->vectU = (double*)0;
	}

	if( a_layout->mgLevel == fftLevel ) {

                if ( a_layout->fftU ) {
                    delete(a_layout->fftU);
                }
                if ( a_layout->fftF ) {
                    delete(a_layout->fftF);
                }
                if ( a_layout->fftBuff ) {
                    delete(a_layout->fftBuff);
                }

#ifdef USE_MPI
		if( a_layout->globalGridDataType != (MPI_Datatype*)0 ) {
			delete [] a_layout->globalGridDataType;
			a_layout->globalGridDataType = (MPI_Datatype*)0;
		}
#endif

	}  // end if mgLevel

}

// =====================================================================

// ---------------------------------------------------------------------

int MGPoisson::setNumGrid( int *a_numGrid ) {

	int iErr = _ELSE_ERR_;
	int numGridX, tmp_numElem;
	int idim;

	if( a_numGrid == (int*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		numGridX = a_numGrid[_x_];
		tmp_numElem = 1;
		for( idim = 0; idim < _3D_; idim++ ) {
			if( a_numGrid[ idim ] < 2 ) {
				iErr = _ILLEGAL_VALUE_;
				break;
			}
			else if( numGridX != a_numGrid[ idim ] ) {
				iErr = _ILLEGAL_VALUE_;
				break;
			}
			else {
				tmp_numElem *= a_numGrid[ idim ];
			}
		}

		if( idim == _3D_ ) {
			numGrid = a_numGrid;
			numElem = tmp_numElem;
			iErr = _NO_ERR_;
		}

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::setSize( double *a_size ) {

	int iErr = _ELSE_ERR_;
	int idim;

	if( a_size == (double*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		for( idim = 0; idim < _3D_; idim++ ) {
			if( a_size[ idim ] < 0.0 ) {
				iErr = _ILLEGAL_VALUE_;
				break;
			}
		}

		if( idim == _3D_ ) {
			size = a_size;
			iErr = _NO_ERR_;
		}

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::setNumLevel( int a_numLevel ) {

	int iErr = _ELSE_ERR_;

	int tmp_numGrid;
	int tmp_numLevel;

	if( a_numLevel < 0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( numGrid == (int*)0 ) {
		iErr = _NOT_INIT_;
	}
	else {

		tmp_numGrid = numGrid[_x_];
		tmp_numLevel = 1;
		while( tmp_numGrid%2 == 0 ) {
			tmp_numLevel++;
			tmp_numGrid /= 2;
		}

		if( tmp_numGrid != 1 ) {
			iErr = _ILLEGAL_VALUE_;
		}
		else if( tmp_numLevel < a_numLevel ) {
			iErr = _ILLEGAL_VALUE_;
		}
		else if( a_numLevel == 0 ) {
			numLevel = tmp_numLevel;
			numLayout = tmp_numLevel;
			iErr = _NO_ERR_;
		}
		else {
			numLevel = a_numLevel;
			numLayout = a_numLevel;
			iErr = _NO_ERR_;
		}

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::setFftLevel( int a_fftLevel ) {

	int iErr = _ELSE_ERR_;

	if( a_fftLevel >= numLevel ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {
		fftLevel = a_fftLevel;
		iErr = _NO_ERR_;
	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::setFlagBC( int *a_flagBC ) {

	int iErr = _ELSE_ERR_;

	int ibc;

#ifdef USE_MPI
	MPI_Status mpi_status;

	int mpi_np;
	int mpi_err;

	int MRank, PRank;
	int tmpMRank, tmpPRank;
	int Midx, Pidx;

	int idx;
	int tmpIdx;
	int count;
	int idim;
#endif

	if( a_flagBC == (int*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {
		for( ibc = 0; ibc < _numBC_; ibc++ ) {
			flagBC[ ibc ] = a_flagBC[ ibc ];
		}

		iErr = _NO_ERR_;

#ifdef USE_MPI

		if( comm == MPI_COMM_NULL ) {
			iErr = _NOT_INIT_;
		}
		else if( myRank == MPI_PROC_NULL ) {
			iErr = _NOT_INIT_;
		}
		else if( ( mpi_err = MPI_Comm_size( comm, &mpi_np ) ) != MPI_SUCCESS ) {
			iErr = mpi_err;
		}
		else {

			procMx = new(std::nothrow) int [ mpi_np ];
			if( procMx == (int*)0 ) {
				iErr = _ALLOC_ERR_;
			}
			else {

				procMx[ 0 ] = flagBC[_centerBC_];

				idx = 1;
				count = 1;
				for( idim = _3D_ - 1; idim >= 0; idim-- ) {

					Midx = ( (idim != _x_)*_numBCy_ + (idim != _y_) )*_numBCz_ + (idim != _z_);
					Pidx = ( ((idim == _x_)+1)*_numBCy_ + ((idim == _y_)+1) )*_numBCz_ + ((idim == _z_)+1);
					MRank = flagBC[ Midx ];
					PRank = flagBC[ Pidx ];

					numProc[ idim ] = 1;
					while( PRank != myRank ) {

						mpi_err = MPI_Sendrecv( &procMx[ 0 ], count, MPI_INTEGER, MRank, 1,
								&procMx[ idx ], count, MPI_INTEGER, PRank, 1,
								comm, &mpi_status );
						if( mpi_err != MPI_SUCCESS ) {
							break;
						}

						mpi_err = MPI_Sendrecv( &flagBC[ Pidx ], 1, MPI_INTEGER, MRank, 2,
								&tmpPRank, 1, MPI_INTEGER, PRank, 2,
								comm, &mpi_status );
						if( mpi_err != MPI_SUCCESS ) {
							break;
						}

						mpi_err = MPI_Sendrecv( &flagBC[ Midx ], 1, MPI_INTEGER, PRank, 3,
								&tmpMRank, 1, MPI_INTEGER, MRank, 3,
								comm, &mpi_status );
						if( mpi_err != MPI_SUCCESS ) {
							break;
						}

						MRank = tmpMRank;
						PRank = tmpPRank;

						numProc[ idim ]++;
						idx += count;

					}  // end while PRank

					if( mpi_err != MPI_SUCCESS ) {
						break;
					}

					count *= numProc[ idim ];

				} // end for idim

				if( mpi_err != MPI_SUCCESS ) {
					iErr = mpi_err;
				}
				else {
					for( idx = 0; idx < mpi_np; idx++ ) {
						if( procMx[ idx ] == _ORIGINAL_RANK_ ) {
							tmpIdx = idx;
							for( idim = _3D_ - 1; idim >= 0; idim-- ) {
								idxProc[ idim ] = ( numProc[ idim ] - tmpIdx%numProc[ idim ] )%numProc[ idim ];
								tmpIdx /= numProc[ idim ];
							}
							break;
						}
					}
					iErr = _NO_ERR_;
				}

			}  // end if procMx

		}  // end if

#endif

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::setVectU0( double *a_vectU0 ) {

	int iErr = _ELSE_ERR_;

	if( a_vectU0 == (double*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {
		vectU0 = a_vectU0;
		iErr = _NO_ERR_;
	}

	return iErr;
}

// ---------------------------------------------------------------------

int MGPoisson::setVectF0( double *a_vectF0 ) {

	int iErr = _ELSE_ERR_;

	if( a_vectF0 == (double*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {
		vectF0 = a_vectF0;
		iErr = _NO_ERR_;
	}

	return iErr;

}

#ifdef USE_MPI
// ---------------------------------------------------------------------

int MGPoisson::setComm( MPI_Comm a_comm ) {

	int iErr = _ELSE_ERR_;
	int mpi_err;

	if( a_comm == MPI_COMM_NULL ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {
		comm = a_comm;
		mpi_err = MPI_Comm_rank( comm, &myRank );
		if( mpi_err != MPI_SUCCESS ) {
			iErr = mpi_err;
		}
		else {
			iErr = _NO_ERR_;
		}
	}

	return iErr;

}

#endif

// ---------------------------------------------------------------------

int MGPoisson::getNumLevel() {
	return numLevel;
}

// ---------------------------------------------------------------------

double MGPoisson::getResidL2Norm() {
	return RESID_L2NORM;
}

// ---------------------------------------------------------------------

double MGPoisson::getConvFact() {
	return CONVFACT;
}

// ---------------------------------------------------------------------

int MGPoisson::getNumIteration() {
	return NUMITERATION;
}

// ---------------------------------------------------------------------

double * MGPoisson::getTimer() {
	return timer;
}

// =====================================================================

// ---------------------------------------------------------------------

void MGPoisson::cleanTimer() {

	for( int idx = 0; idx < _MAX_INDEX_TIMER_; idx++ ) {
		timer[ idx ] = 0.0;
	}

}

// ---------------------------------------------------------------------

void MGPoisson::printTimer( FILE *fp ) {

	fprintf( fp, "# TOTAL CALC MPI-1 MPI-2 FFT FFT_MPI ELSE\n" );
	fprintf( fp, "#" );
	for( int idx = 0; idx < _MAX_INDEX_TIMER_; idx++ ) {
		fprintf( fp, " %f", timer[ idx ] );
	}
	fprintf( fp, "\n" );

}

// ---------------------------------------------------------------------

void MGPoisson::printTimerAvg( FILE *fp ) {

	fprintf( fp, "# TOTAL CALC MPI-1 MPI-2 FFT FFT_MPI ELSE (/ITERATION)\n" );
	fprintf( fp, "#" );
	for( int idx = 0; idx < _MAX_INDEX_TIMER_; idx++ ) {
		fprintf( fp, " %f", timer[ idx ]/NUMITERATION );
	}
	fprintf( fp, "\n" );

}

// ---------------------------------------------------------------------

int MGPoisson::printLayout( FILE *fp, int a_iLayout ) {

	int iErr = _ELSE_ERR_;

	int idim;
	int iband;
	mgp_layout_struct *tmp_layout;

	if( fp == (FILE*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( a_iLayout < 0 || a_iLayout >= numLayout ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( layout[ a_iLayout ] == (mgp_layout_struct*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		tmp_layout = layout[ a_iLayout ];

		fprintf( fp, "# ----- %d Gridlayout menber -----\n", a_iLayout );

		fprintf( fp, "# layout pointer\n" );
		fprintf( fp, "#\tthis = %p\n", tmp_layout );
		fprintf( fp, "#\tfine = %p\n", tmp_layout->fine );
		fprintf( fp, "#\tcoarse = %p\n", tmp_layout->coarse );

		fprintf( fp, "# size\n" );
		for( idim = 0; idim < _3D_; idim++ ) {
			fprintf( fp, "#\t[%d] : %f\n", idim, tmp_layout->size[ idim ] );
		}

		fprintf( fp, "# numGrid\n" );
		for( idim = 0; idim < _3D_; idim++ ) {
			fprintf( fp, "#\t[%d] : %d\n", idim, tmp_layout->numGrid[ idim ] );
		}

		fprintf( fp, "# flagBC\n" );
		for( idim = 0; idim < _3D_; idim++ ) {
			int idxD = (int)pow(3,_3D_-(idim+1));
			int idx = 13 - idxD;
			fprintf( fp, "#\t" );
			for( int idir = 0; idir < 3; idir++ ) {
				fprintf( fp, "(%d,%d,%d)=%d ", idim, idir, idx, tmp_layout->flagBC[ idx ] );
				idx += idxD;
			}
			fprintf( fp, "\n" );
		}

		fprintf( fp, "# deltaIndex\n" );
		for( iband = 0; iband < 2*_3D_ + 1; iband++ ) {
			fprintf( fp, "#\t[%d] : %d\n", iband, tmp_layout->deltIndex[ iband ] );
		}

		fprintf( fp, "# spaceGrid\n" );
		for( idim = 0; idim < _3D_; idim++ ) {
			fprintf( fp, "#\t[%d] : %f\n", idim, tmp_layout->spaceGrid[ idim ] );
		}

		fprintf( fp, "# mxCoef\n" );
		for( iband = 0; iband < 2*_3D_ + 1; iband++ ) {
			fprintf( fp, "#\t[%d] : %f\n", iband, tmp_layout->mxCoef[ iband ] );
		}

		fprintf( fp, "# gmLevel = %d\n", tmp_layout->mgLevel );
		fprintf( fp, "# gmNumRecall = %d\n", tmp_layout->mgNumRecall );
		fprintf( fp, "# gmNumPre = %d\n", tmp_layout->mgNumPre );
		fprintf( fp, "# gmNumPost = %d\n", tmp_layout->mgNumRecall );

		fprintf( fp, "# transSize, sendStart -> recvStart\n" );
		for( int idir = 0; idir < 2*_3D_; idir++ ) {
			fprintf( fp, "#\t[%d]:(", idir );
			for( idim = 0; idim < _3D_; idim++ ) {
				fprintf( fp, " %d", tmp_layout->transSize[ idir ][ idim ] );
			}
			fprintf( fp, " ), (" );
			for( idim = 0; idim < _3D_; idim++ ) {
				fprintf( fp, " %d", tmp_layout->sendStart[ idir ][ idim ] );
			}
			fprintf( fp, " ) -> (" );
			for( idim = 0; idim < _3D_; idim++ ) {
				fprintf( fp, " %d", tmp_layout->recvStart[ idir ][ idim ] );
			}
			fprintf( fp, " )\n" );
		}
		fprintf( fp, "\n" );

		fprintf( fp, "# --------------------------------\n" );

		iErr = _NO_ERR_;

	}

	return iErr;

}

// ---------------------------------------------------------------------

int MGPoisson::printAllLayout( FILE *fp ) {

	int iErr = _ELSE_ERR_;

	int iLayout;

	if( fp == (FILE*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else {

		for( iLayout = 0; iLayout < numLayout; iLayout++ ) {
			if( ( iErr = printLayout( fp, iLayout ) ) != _NO_ERR_ ) {
				break;
			}
		}

	}

	return iErr;

}

#ifdef USE_MPI
// ---------------------------------------------------------------------

int MGPoisson::printProcMx( FILE *fp ) {

	int iErr = _ELSE_ERR_;

	int ix, iy, iz;
	long idx;

	if( fp == (FILE*)0 ) {
		iErr = _ILLEGAL_VALUE_;
	}
	else if( procMx == (int*)0 ) {
		iErr = _NOT_INIT_;
	}
	else {

		fprintf( fp, "---------- %3d ----------\n", myRank );
		fprintf( fp, "  X->\n" );
		fprintf( fp, "  Z | Y->\n" );
		for( iz = 0; iz < numProc[_z_]; iz++ ) {
			fprintf( fp, " %2d |", iz );
			for( ix = 0; ix < numProc[_x_]; ix++ ) {
				for( iy = 0; iy < numProc[_y_]; iy++ ) {
					idx = getIndex( ix, iy, iz, numProc );
					fprintf( fp, " %2d", procMx[ idx ] );
				}
				fprintf( fp, " |" );
			}
			fprintf( fp, "\n" );
		}

		fprintf( fp, "-------------------------\n" );

		iErr = _NO_ERR_;

	}

	return iErr;

}
#endif
