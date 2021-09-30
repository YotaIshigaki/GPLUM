#ifndef _MGP_LAYOUT_STRUCT_H_
#define _MGP_LAYOUT_STRUCT_H_

#include "UseMPI.h"
#include "MGFFT3D.h"

/**
 * マルチグリッドレイアウト構造体
 */
namespace MultigridModule {

struct mgp_layout_struct {

	mgp_layout_struct *fine;	/**< 細かいグリッドのマルチグリッドレイアウト構造体へのポインタ */

	mgp_layout_struct *coarse;	/**< 粗いグリッドのマルチグリッドレイアウト構造体へのポインタ */

	double size[_3D_];		/**< 計算領域の寸法を与える配列 (サイズ=3, [0]=x, [1]=y, [2]=z ) */

	int numGrid[_3D_];		/**< 計算領域のグリッド数を与える配列 (サイズ=3, [0]=x, [1]=y, [2]=z ) */

	int flagBC[_numBC_];		/**< 計算領域に隣接する境界面/境界条件を指定する配列へのポインタ (サイズ27, [4]=(-1,0,0)面, [22]=(1,0,0)面, [10]=(0,-1,0)面, [16]=(0,1,0)面, [12]=(0,0,-1)面, [14]=(0,0,1)面) */

	double *vectF;		/**< 定数ベクトル配列 (サイズ = numGrid[0]xnumGrid[1]xnumGrid[2] 1次元 row-majer order 配列 ) */

	double *vectU;		/**< 解ベクトル配列 ( サイズ = numGrid[0]xnumGrid[1]xnumGrid[2] 1次元 row-majer order 配列 ) */

	double *vectR;		/**< 残差ベクトル配列 ( サイズ = numGrid[0]xnumGrid[1]xnumGrid[2] 1次元 row-majer order 配列 ) */

	double spaceGrid[_3D_];	/**< グリッド間隔を与える配列 ( サイズ = 3, [0]=x, [1]=y, [2]=z ) */

	double mxCoef[2*_3D_+1];	/**< 係数行列 ( サイズ = 7 ) */

	int deltIndex[2*_3D_+1];	/**< ベクトル配列のインデックス間隔配列 ( サイズ = 7 ) */

	long numElem;		/**< ベクトル配列の要素サイズ ( numGrid[0] x numGrid[1] x numGrid[2] ) */

	double omega;		/**< SORの緩和係数 ( デフォルト = 1.0 ) */

	int ncolor;			/**< SORのカラー数 ( デフォルト = 2 ) # 2まで対応 */

	int mgLevel;		/**< マルチグリッドレベル番号 */

	int mgNumRecall;	/**< VCycleの再帰呼び出し回数 ( デフォルト = 1 ) */

	int mgNumPre;		/**< VCycleのPre-smoothing回数 ( デフォルト = 2 ) */

	int mgNumPost;		/**< VCycleのPost-smoothing回数 ( デフォルト = 2 ) */

	int transSize[2*_3D_][_3D_];	/**< 仮想境界転送サイズを与える配列 ( サイズ = [6][3], 第1インデックスは仮想境界面( [0]=(-1,0,0)面, [1]=(0,-1,0)面, [2]=(0,0,-1)面, [3]=(+1,0,0)面, [4]=(0,+1,0)面, [5]=(0,0,+1)面 )を、第2インデックスは次元( [0]=x, [1]=y, [2]=z )を表す ) */

	int sendStart[2*_3D_][_3D_];	/**< ベクトル配列での仮想境界転送の送信開始位置を与える配列 ( サイズ = [6][3], 第1インデックスは仮想境界面( [0]=(-1,0,0)面, [1]=(0,-1,0)面, [2]=(0,0,-1)面, [3]=(+1,0,0)面, [4]=(0,+1,0)面, [5]=(0,0,+1)面 )を、第2インデックスは次元( [0]=x, [1]=y, [2]=z )を表す ) */

	int recvStart[2*_3D_][_3D_];	/**< ベクトル配列での仮想境界転送の受信開始位置を与える配列 ( サイズ = [6][3], 第1インデックスは仮想境界面( [0]=(-1,0,0)面, [1]=(0,-1,0)面, [2]=(0,0,-1)面, [3]=(+1,0,0)面, [4]=(0,+1,0)面, [5]=(0,0,+1)面 )を、第2インデックスは次元( [0]=x, [1]=y, [2]=z )を表す ) */

   // FFTW2
	int fftNumGrid[_3D_];
        FFT3D::Ptr fftF;
        FFT3D::Ptr fftU;
        double *fftBuff;

#ifdef USE_MPI
	MPI_Datatype sendType[2*_3D_];		/**< 仮想境界転送の送信データタイプを与える配列 ( サイズ = 6, [0]=(-1,0,0)面, [1]=(0,-1,0)面, [2]=(0,0,-1)面, [3]=(+1,0,0)面, [4]=(0,+1,0)面, [5]=(0,0,+1)面 ) */
	MPI_Datatype recvType[2*_3D_];		/**< 仮想境界転送の受信データタイプを与える配列 ( サイズ = 6, [0]=(-1,0,0)面, [1]=(0,-1,0)面, [2]=(0,0,-1)面, [3]=(+1,0,0)面, [4]=(0,+1,0)面, [5]=(0,0,+1)面 ) */

	int globalGridDisp[_3D_];		/**< 全体グリッドデータ内に含まれる部分グリッドデータの位置を示す配列（ サイズ = 3, [0]=x, [1]=y, [2]=z ) */

	MPI_Datatype localGridDataType;		/**< 部分グリッドデータ配列（仮想境界を除いた）の送受信データタイプ */
	MPI_Datatype *globalGridDataType;	/**< 全体グリッドデータ内に含まれる各CPUの部分グリッドデータの送受信タイプ */

#endif

};
}
#endif // _MGP_LAYOUT_STRUCT_H_
