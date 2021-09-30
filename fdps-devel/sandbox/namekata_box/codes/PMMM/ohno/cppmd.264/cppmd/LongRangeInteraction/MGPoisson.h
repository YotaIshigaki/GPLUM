#ifndef _MGPOISSON_H_
#define _MGPOISSON_H_

#include "UseMPI.h"
#ifndef USE_MPI
#include <time.h>
#define MPI_Wtime() ( (double)clock()/(double)CLOCKS_PER_SEC )
#endif

#define _2D_ 2		/**< 2次元を表す定数 */
#define _3D_ 3		/**< 3次元を表す定数 */

#define _numBC_ 27	/**< 隣接する境界面/境界条件を指定する配列サイズの定数 */
#define _numBCx_ 3
#define _numBCy_ 3
#define _numBCz_ 3
#define _centerBC_ 13	/**< 隣接する境界面/境界条件を指定する配列の中心インデックス */

#define _x_ 0		/**< x次元を指定するインデックス定数 */
#define _y_ 1		/**< y次元を指定するインデックス定数 */
#define _z_ 2		/**< z次元を指定するインデックス定数 */

#define _JACOB_ 0		/**< JACOB法を表す定数 */
#define _SOR_ 1		/**< SOR法を表す定数 */
#define _VCYC_ 2		/**< V-Cycle法を表す定数 */
#define _FMG_ 3		/**< FMG法を現す定数 */

#define _DIRECHLET_ -1		/**< 固定境界条件を表す定数 */
#define _PERIODIC_ -2		/**< 周期境界条件を表す定数 */

#define _ORIGINAL_RANK_ 0	/**< 並列分割の原点 */

#define _VECT_U_ 0			/**< 解ベクトルを指定するインデックス */
#define _VECT_F_ 1			/**< 定数ベクトルを指定するインデックス */
#define _VECT_R_ 2			/**< 残差ベクトルを指定するインデックス */

#define _NO_ERR_ 0			/**< エラー番号: 正常終了 */
#define _ILLEGAL_VALUE_ -1	/**< エラー番号：不正な値の使用 */
#define _NOT_INIT_ -2		/**< エラー番号：初期化していない */
#define _ALLOC_ERR_ -3		/**< エラー番号：メモリ確保に失敗 */
#define _ELSE_ERR_ -9999	/**< エラー番号：未定義のエラー */

#define _MAX_INDEX_TIMER_ 7	/**< タイマー配列のサイズ */
#define _TOTAL_TIME_ 0			/**< タイマー配列の全実行時間を指定するインデックス */
#define _CALC_TIME_  1			/**< タイマー配列の計算部時間を指定するインデックス */
#define _MPI1_TIME_  2			/**< タイマー配列のSendRecv時間を指定するインデックス */
#define _MPI2_TIME_  3			/**< タイマー配列のAllReduce時間を指定するインデックス */
#define _ELSE_TIME_  4			/**< タイマー配列のその他の時間を指定するインデックス */

#define _FFT_TIME_ 5				/**< タイマー配列のFFT計算時間を指定するインデックス */
#define _FFT_MPI_TIME_ 6		/**< タイマー配列のFFTでのSendRecv時間を指定するインデックス */


#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <iostream>

#include "mgp_layout_struct.h"

/**
 * マルチグリッドポアッソンソルバークラス.
 */

namespace MultigridModule {

class MGPoisson {

	public:

		/**
		 * デフォルトコンストラクタ
		 */
		MGPoisson();

		/**
		 * デフォルトデストラクタ
		 */
		~MGPoisson();

		//
		// Interface
		//

		/**
		 * 初期化（シリアル版）.
		 *
		 * @param[in] *a_numGrid - 計算領域のグリッド数を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 * @param[in] *a_size - 計算領域の寸法を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 * @param[in,out] *a_vectU0 - 解ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 * @param[in] *a_vectF0 - 定数ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 *
		 * @return エラー番号
		 */
		int Initialize( int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0, int a_fftLevel );

		int Initialize( int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0 );

		int Initialize( int *a_flagBC,
				int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0, int a_fftLevel );

		int Initialize( int *a_flagBC,
				int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0 );

		/**
		 * マルチグリッドパラメータの初期化.
		 *
		 * @param[in] a_numRecall - 再帰呼び出し回数
		 * @param[in] a_numPre - Pre-Smoothing回数
		 * @param[in] a_numPost - Post-Smoothing回数
		 *
		 * @return エラー番号
		 */
		int InitVCycle( int a_numRecall, int a_numPre, int a_numPost );

#ifdef USE_MPI
		/**
		 * 初期化（MPI並列版）.
		 *
		 * @param[in] a_comm - MPIコミュニケータ
		 * @param[in] *a_RankArr - 隣接領域の計算を担当しているプロセスの番号を与える配列へのポインタ (サイズ = 27, [4]=(-1,0,0)面, [22]=(1,0,0)面, [10]=(0,-1,0)面, [16]=(0,1,0)面, [12]=(0,0,-1)面, [14]=(0,0,1)面)
		 * @param[in] *a_numGrid - このプロセスが担当する計算領域のグリッド数を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 * @param[in] *a_size - このプロセスが担当する計算領域の寸法を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 * @param[in,out] *a_vectU0 -このプロセスが担当する計算領域の解ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 * @param[in] *a_vectF0 - このプロセスが担当する計算領域の定数ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 */
		int Initialize( MPI_Comm a_comm, int *a_RankArr,
				int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0, int a_fftLevel );

		int Initialize( MPI_Comm a_comm, int *a_RankArr,
				int *a_numGrid, double *a_size,
				double *a_vectU0, double *a_vectF0 );
#endif

		//
		// Solver
		//

		/**
		 * 指定した繰り返し回数（numIter)だけSORを実行.
		 * 
		 * @param[in] numIter - 繰り返し回数
		 *
		 * @return エラー番号
		 */
		int SORNIter( int numIter );

		/**
		 * 指定した繰り返し回数（numIter）だけVCycleを実行.
		 *
		 * @param[in] numIter - 繰り返し回数
		 *
		 * @return エラー番号
		 */
		int VCycleNIter( int numIter );

		/**
		 * 指定した再帰呼び出し回数でFMGを実行.
		 *
		 * @param [in] numRecall - 再帰呼び出し回数
		 *
		 * @return エラー番号
		 */
		int FMG( int numRecall );

		//
		// Seter and Getter
		//

		/**
		 * 計算領域のグリッド数を与える配列へのポインタの設定.
		 *
		 * @param[in] *a_numGrid - 計算領域のグリッド数を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 *
		 * @return エラー番号
		 */
		int setNumGrid( int *a_numGrid );

		/**
		 * 計算領域の寸法を与える配列へのポインタの設定.
		 *
		 * @param[in] *a_size - 計算領域の寸法を与える配列へのポインタ (サイズ=3, [0]=x, [1]=y, [2]=z)
		 *
		 * @return エラー番号
		 */
		int setSize( double *a_size );

		/**
		 * 解ベクトル配列へのポインタの設定.
		 *
		 * @param[in,out] *a_vectU0 - 解ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 *
		 * @return エラー番号
		 */
		int setVectU0( double *a_vectU0 );

		/**
		 * 定数ベクトル配列へのポインタの設定.
		 *
		 * @param[in] *a_vectF0 - 定数ベクトル配列へのポインタ (サイズ= a_numGrid[x] x a_numGrid[y] x a_numGrid[z] の1次元 row-major order 配列)
		 *
		 * @return エラー番号
		 */
		int setVectF0( double *a_vectF0 );

		/**
		 * マルチグリッドレベル数の設定.
		 *
		 * @param[in] a_numLevel - マルチグリッドレベル数
		 *
		 * @return エラー番号
		 */
		int setNumLevel( int a_numLevel );

		/**
		 * 
		 *
		 */
		int setFlagBC( int *a_flagBC );

		/**
		 * FFTを適用するマルチグリッドレベルの設定.
		 *
		 * @param[in] a_fftLevel - FFT適用マルチグリッドレベル
		 *
		 * @return エラー番号
		 */
		int setFftLevel( int a_fftLevel );

		/**
		 * マルチグリッドレベル数の取得.
		 *
		 * @return マルチグリッドレベル数
		 */
		int getNumLevel();

		/**
		 * 残差L2ノルムの取得.
		 *
		 * @return 残差L2ノルム
		 */
		double getResidL2Norm();

		/**
		 * 収束係数の取得.
		 *
		 * @return 収束係数
		 */
		double getConvFact();

		/**
		 * 反復回数の取得.
		 *
		 * @return 反復回数
		 */
		int getNumIteration();

		//
		// UTILITY
		//

		/**
		 * タイマー配列のクリアー.
		 */
		void cleanTimer();

		/**
		 * タイマー配列へのポインタの取得.
		 *
		 * @return タイマー配列へのポインタ ( サイズ=_MAX_INDEX_TIMER_ )
		 */
		double *getTimer();

		/**
		 * タイマー配列のファイル出力.
		 *
		 * @param[in] *fp - 出力ファイルポインタ
		 */
		void printTimer( FILE *fp );

		/**
		 * 反復回数で平均したタイマー配列のファイル出力.
		 *
		 * @param[in] *fp - 出力ファイルポインタ
		 */
		void printTimerAvg( FILE *fp );

		/**
		 * 指定したマルチグリッドレイアウト構造体のファイル出力.
		 *
		 * @param[in] *fp - 出力ファイルポインタ
		 * @param[in] a_iLayout - マルチグリッドレイアウト構造体配列のインデックス
		 *
		 * @return エラー番号
		 */
		int printLayout( FILE *fp, int a_iLayout );

		/**
		 * すべてのマルチグリッドレイアウト構造体のファイル出力.
		 *
		 * @param[in] *fp - 出力ファイルポインタ
		 *
		 * @return エラー番号
		 */
		int printAllLayout( FILE *fp );

#ifdef USE_MPI
		/**
		 * 3次元分割されたCPUランク番号の配列のファイル出力.
		 *
		 * @param[in] *fp - 出力ファイルポインタ
		 *
		 * @return エラー番号
		 */
		int printProcMx( FILE *fp );

#endif

//	private:

		double *size;		/**< 計算領域の寸法を与える配列へのポインタ ( サイズ = 3, [0]=x, [1]=y, [2]=z ) */

		int *numGrid;		/**< 計算領域のグリッド数を与える配列へのポインタ ( サイズ = 3, [0]=x, [1]=y, [2]=z ) */

		int flagBC[_numBC_];		/**< 計算領域に隣接する境界面/境界条件を指定する配列へのポインタ (サイズ27, [4]=(-1,0,0)面, [22]=(1,0,0)面, [10]=(0,-1,0)面, [16]=(0,1,0)面, [12]=(0,0,-1)面, [14]=(0,0,1)面) */

		double *vectF0;	/**< 解ベクトル配列へのポインタ (サイズ = numGrid[0]xnumGrid[1]xnumGrid[2] 1次元 row-majer order 配列 ) */

		double *vectU0;	/**< 定数ベクトル配列へのポインタ (サイズ = numGrid[0]xnumGrid[1]xnumGrid[2] 1次元 row-majer order 配列 ) */

		mgp_layout_struct **layout;		/**< マルチグリッドレイアウト構造体配列へのポインタの配列 (*layout[numLayout]) */

		long numElem;		/**< 計算領域の全グリッド数 */

		int numLayout;		/**< マルチグリッドレイアウトの数 */

		int numLevel;		/**< マルチグリッドのレベル数 */

		int numRecall;		/**< VCycleの再帰呼び出し回数 */

		int numPre;			/**< VCycleのPre-smoothing回数 */

		int numPost;		/**< VCycleのPost-smoothing回数 */

		int fmgNumRecall;	/**< FMGの再帰呼び出し回数 */

		int fftLevel;		/**< FFT解法を適用するマルチグリッドレベル */

#ifdef USE_MPI
		MPI_Comm comm;		/**< MPIコミュニケータ */
		int myRank;			/**< プロセス(ランク)番号 */

		int numProc[_3D_];		/**< 各次元の並列数配列 ( サイズ = 3, [0]=x, [1]=y, [2]=z ) */
		int idxProc[_3D_];		/**< ランク０を基準とした自CPUのインデックス配列 ( サイズ=3, [0]=x, [1]=y, [2]=z ) */
		int *procMx;				/**< 3次元分割されたCPUランク番号の配列 ( サイズ numProc[0] x numProc[1] x numProc[3] 1次元 row-majer order 配列 ) */

		/**
		 *
		 */
		int setComm( MPI_Comm a_comm );

#endif

		int isInit;			/**< 初期化状態 */

		int NUMITERATION;	/**< 反復回数 */

		double RESID_L2NORM;	/**< 残差L2ノルム */

		double CONVFACT;		/**< 収束係数 */

		double timer[ _MAX_INDEX_TIMER_ ];	/**< タイマー配列 */

		/**
		 * マルチグリッドレイアウト構造体配列の生成.
		 *
		 * @return エラー番号
		 */
		int createLayout();

		/**
		 * マルチグリッドレイアウト構造体配列の破棄.
		 */
		void destroyLayout();

		/**
		 * マルチグリッドレイアウト構造体の初期化.
		 *
		 * @param[in,out] *a_layout - 初期化するマルチグリッドレイアウト構造体へのポインタ
		 * @param[in] a_level - 初期化するマルチグリッドレイアウト構造体のマルチグリッドレベル
		 * @param[in] *a_numGrid - 初期化するマルチグリッドレイアウト構造体にグリッド数を与える配列へのポインタ
		 * @param[in] *a_size - 初期化するマルチグリッドレイアウト構造体に寸法を与える配列へのポインタ
		 *
		 * @return エラー番号
		 */
		int allocLayout( mgp_layout_struct *a_layout,
				int a_level, int *a_numGrid, double *a_size );

		int allocFFTLayout( mgp_layout_struct *a_layout, int *a_numGrid );

		/**
		 * マルチグリッドレイアウト構造体の破棄.
		 *
		 * @param[in] *a_layout - 破棄するマルチグリッドレイアウト構造体
		 */
		void freeLayout( mgp_layout_struct *a_layout );

		/**
		 * 解ベクトル配列ポインタと定数ベクトル配列ポインタの指す配列の値で、マルチグリッドレイアウト構造体の解ベクトル配列と定数ベクトル配列を更新.
		 */
		void renewVector();

		/**
		 * マルチグリッドレイアウト構造体の解ベクトル配列の値を、解ベクトル配列ポインタの指す配列に返す.
		 */
		void returnVector();

		int copyVector( mgp_layout_struct *a_layout, int a_flgVectFrom, int a_flgVectTo );

		/**
		 * 指定したマルチグリッドレイアウト構造体でSORを実行.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int SORCore( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体でVCycleを実行.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int VCycleCore( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体でFMGを実行.
		 *
		 * @param[in] a_numReacll - 再帰呼び出し回数
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int FMGCore( int a_numRecall, mgp_layout_struct *a_layout );

		/**
		 *
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * preturn エラー番号
		 */
		int FFTPoisson( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体で指定した繰り返し回数だけSORを実行.
		 *
		 * @param[in] a_iter - 繰り返し回数
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int SORNIter( int a_iter, mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体で残差を計算.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int CalcResidual( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体に対して、細かいグリッドから粗いグリッドへの残差割付.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int Restriction( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体に対して、粗いグリッドから細かいグリッドへの解ベクトルの修正.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体へのポインタ
		 *
		 * @return エラー番号
		 */
		int Correction( mgp_layout_struct *a_layout );

		/**
		 * 指定したマルチグリッドレイアウト構造体の指定したベクトル配列の仮想境界面を更新.
		 *
		 * @param[in] a_dir - 更新する仮想境界面の指定 ( 0=(-1,0,0)面, 1=(0,-1,0)面, 2=(0,0,-1)面, 3=(+1,0,0)面, 4=(0,+1,0)面, 5=(0,0,+1)面 )
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体
		 * @param[in] a_flgVect - 更新するベクトル配列の指定 ( 0=解ベクトル、1=定数ベクトル、2=残差ベクトル )
		 *
		 * @return エラー番号
		 */
		int transBC( int a_dir, mgp_layout_struct *a_layout, int a_flgVect );


		/**
		 * 指定したマルチグリッドレイアウト構造体の指定したベクトル配列をFFT用の配列に転送.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体
		 * @param[in] a_flgVect - FFT配列に転送するベクトル配列の指定 ( 0=解ベクトル、1=定数ベクトル )
		 *
		 * @return エラー番号
		 */
		int transMGtoFFT( mgp_layout_struct *a_layout, int a_flgVect );

		/**
		 * 指定したマルチグリッドレイアウト構造体の指定したFFT用の配列をベクトル配列に転送.
		 *
		 * @param[in,out] *a_layout - マルチグリッドレイアウト構造体
		 * @param[in] a_flgVect - FFT配列を転送するベクトル配列の指定 ( 0=解ベクトル、1=定数ベクトル )
		 *
		 * @return エラー番号
		 */
		int transFFTtoMG( mgp_layout_struct *a_layout, int a_flgVect );

		/**
		 * 3次元インデックスから1次元row-major order配列のインデックスを取得.
		 *
		 * @param[in] ix - x次元のインデックス
		 * @param[in] iy - y次元のインデックス
		 * @param[in] iz - z次元のインデックス
		 * @param[in] *n - 3次元表現した場合の配列の大きさ (サイズ=3, [0]=x, [1]=y, [2]=z )
		 *
		 * @return 1次元row-major order配列のインデックス
		 */
		long getIndex( int ix, int iy, int iz, int *n ) {
			return iz + n[_z_]*( iy + n[_y_]*ix );
		}

		/**
		 * 3次元インデックスから1次元row-major order配列のインデックスを取得.
		 *
		 * @param[in] ix - x次元のインデックス
		 * @param[in] iy - y次元のインデックス
		 * @param[in] iz - z次元のインデックス
		 * @param[in] nx - x次元の大きさ
		 * @param[in] ny - y次元の大きさ
		 * @param[in] nz - z次元の大きさ
		 *
		 * @return 1次元row-major order配列のインデックス
		 */
		long getIndex( int ix, int iy, int iz, int nx, int ny, int nz ) {
			return iz + nz*( iy + ny*ix );
		}

		/**
		 *
		 */
		void getIndex3D( int idx1d, int *n, int *idx3d ) {
			int idx;
			int idim;

			idx = idx1d;
			for( idim = _3D_ - 1; idim >= 0; idim-- ) {
				idx3d[ idim ] = idx%n[ idim ];
				idx /= n[ idim ];
			}

		}
};
}
#endif // _MGPOISSON_H_
