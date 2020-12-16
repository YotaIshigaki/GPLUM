# 開発環境構築メモ

### この文書の内容
- [md_fdpsの実行に必要な環境について](#md_fdps_env)
- [インストールの準備](#インストールの準備)
- [GCC のコンパイルとインストール](#GCCのコンパイルとインストール)
- [OpenMPI のコンパイルとインストール](#OpenMPIのコンパイルとインストール)
- [FFTw のコンパイルとインストール](#FFTwのコンパイルとインストール)
- [おまけ1：　FDPSのParticleMesh機能拡張のコンパイル](#FDPSのparticle_mesh拡張機能のコンパイル)
- [おまけ2：　複数のGCCバージョンの切り替え](#複数のgccバージョンの切り替え)
- [おまけ3：　Google C++ test framework (google test)のインストール](#gtestのインストール)

<a id="md_fdps_env"></a>
<a href="#md_fdps_env"></a>
## md_fdpsの実行に必要な環境について
FDPS (Framework for Developing Particle Simulator) を使用して分子動力学プログラムを作成するため，FDPSを利用したプログラムをコンパイルするために必要な環境を構築する．

FDPS ver4.0b 向け  
参考：
- https://github.com/FDPS/FDPS
- http://www.jmlab.jp/?p=530

#### FDPS開発チームの動作確認環境 (FDPSのチュートリアル参照)
C++ コンパイラ (gcc 4.4.5以降，あるいはKコンパイラ 1.2.0)  
MPI version 1.3環境　(OpenMPI 1.8.1で動作確認)  
FFTw 3.3以降  

クラスタ計算機システムの環境の例  
- シェル bash
- gcc 4.3.4  
MPI環境:  
  - OpenMPI 1.5.3  
  - Intel MPI 4.1.1.036  

GCCが古すぎるのもあるが，C++11以降の規格が使用できると便利なこともあり，開発環境として以下の環境をユーザのローカルフォルダにインストールした．
- gcc 6.3.0
- OpenMPI 2.1.1
- FFTW 3.3.6


<a id="インストールの準備"></a>
<a href="#インストールの準備"></a>
## インストールの準備
適当な作業フォルダ(ここでは `$HOME/install/` )に上記各ソフトウェアの圧縮ソースファイルをダウンロードし，展開しておく．  
下記で例示しているディレクトリ( `/home/hogehoge` )は作業者のホームディレクトリである．各自 `$ pwd` コマンド等で確認すること．  

```bash
$ pwd
/home/hogehoge/install

$ ls
fftw-3.3.6-pl1.tar.gz
gcc-6.3.0.tar.bz2
openmpi-2.1.1.tar.bz2

$ tar -xvf gcc-6.3.0.tar.bz2
$ tar -xvf openmpi-2.1.1.tar.bz2
$ tar -zxvf fftw-3.3.6-pl1.tar.gz
```

既存環境(ここではgcc-4.3)との切り替えのため，コマンドをシンボリックリンクとして作成する．  
パスを登録するディレクトリ(ここでは `$HOME/local_path` )を作成

```bash
$ cd
$ mkdir local_path
$ cd local_path
```

既存環境のコンパイラのシンボリックリンクの作成

```bash
$ pwd
/home/hogehoge/local_path

$ ln -s /usr/bin/gcc gcc-4.3
$ ln -s /usr/bin/g++ g++-4.3
$ ln -s /usr/bin/gfortran gfortran-4.3
```

シンボリックリンクのパスの登録  
.bashrc の最後に以下を追加．　あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```bash
export PATH=$HOME/local_path:$PATH
```
.bashrcを変更後，__ログインしなおして__ 登録したシンボリックリンク(例えば `gcc-4.3` )が有効であることを確認する．
```bash
$ gcc-4.3 -v
Using built-in specs.
Target: x86_64-suse-linux
Configured with: ../configure --prefix=/usr --infodir=/usr/share/info --mandir=/usr/share/man --libdir=/usr/lib64 --libexecdir=/usr/lib64 --enable-languages=c,c++,objc,fortran,obj-c++,java,ada --enable-checking=release --with-gxx-include-dir=/usr/include/c++/4.3 --enable-ssp --disable-libssp --with-bugurl=http://bugs.opensuse.org/ --with-pkgversion='SUSE Linux' --disable-libgcj --disable-libmudflap --with-slibdir=/lib64 --with-system-zlib --enable-__cxa_atexit --enable-libstdcxx-allocator=new --disable-libstdcxx-pch --enable-version-specific-runtime-libs --program-suffix=-4.3 --enable-linux-futex --without-system-libunwind --with-cpu=generic --build=x86_64-suse-linux
Thread model: posix
```

コンパイル作業用のディレクトリを作成  

```bash
$ cd
$ cd install

$ mkdir build
$ cd build

$ pwd
/home/hogehoge/install/build
```

これでインストール作業の準備は完了である．

<a id="GCCのコンパイルとインストール"></a>
<a href="#GCCのコンパイルとインストール"></a>
## GCC のコンパイルとインストール
GCC のコンパイルに必要な依存関係ファイルをダウンロードする．
```bash
$ cd ../gcc-6.3.0/
$ ./contrib/download_prerequisites
```
makeファイルを作成する．  
下記は見やすいように改行しているが， `./gcc-6.3.0/configure` 以降のすべてのオプションを改行せずに続けて入力すること．  
(ここでは `$HOME/local/gcc-6.3.0` にインストールする設定を与えている)
```bash
$ cd ../build
$ ../gcc-6.3.0/configure --prefix=$HOME/local/gcc-6.3.0  
--with-local-prefix=$HOME/local/libgcc63  
--enable-checking=release --disable-multilib  
--enable-languages=c,c++,fortran  
```

無事にmakefileが生成されたら，コンパイルおよびインストールを行う．
```bash
$ make
$ make install
```

`$ make` は早くとも30分～1時間程度はかかる．  
無事に `$ make install` が完了したら，`.bashrc` の末尾に以下の設定を追記する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．

```bash
target=$HOME"/local/gcc-6.3.0"  # prefixで指定したGCCのフォルダ
export PATH=${target}/bin:$PATH
export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${target}/lib64:$LD_LIBRARY_PATH
export LD_RUN_PATH=${target}/lib:$LD_RUN_PATH
export LD_RUN_PATH=${target}/lib64:$LD_RUN_PATH
```

.bashrcを変更後，__ログインしなおして__ 新しくインストールしたGCCとパスが有効になっているかどうかを確認する  

```bash
$ gcc -v
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/home/hogehoge/local/gcc-6.3.0/libexec/gcc/x86_64-pc-linux-gnu/6.3.0/lto-wrapper
Target: x86_64-pc-linux-gnu
Configured with: ../gcc-6.3.0/configure --prefix=/home/hogehoge/local/gcc-6.3.0 --with-local-prefix=/home/hogehoge/local/libgcc63 --enable-checking=release --disable-multilib --enable-languages=c,c++,fortran
Thread model: posix
gcc version 6.3.0 (GCC)

$ echo $PATH             # 各々追加したパスが先頭にあるか確認する
$ echo $LD_LIBRARY_PATH
$ echo $LD_RUN_PATH
```
これでGCCのインストールは完了である．
また，`.bashrc` に追記した部分を削除すればデフォルトの環境に戻すこともできる．

<a id="OpenMPIのコンパイルとインストール"></a>
<a href="#OpenMPIのコンパイルとインストール"></a>
## OpenMPI のコンパイルとインストール
既存のbuildディレクトリを使いまわすので gcc 関係のファイルの消去
```bash
$ rm * -r
```
makeファイルの作成  
Fortranを使用しない場合は `F77=***` と `FC=***` は書かない．  
(ここでは `$HOME/local/openmpi-2.1.1` にインストールする)
```bash
$ ../openmpi-2.1.1/configure --prefix=$HOME/local/openmpi-2.1.1 --enable-mpi-cxx CC=gcc CXX=g++ F77=gfortran FC=gfortran
```
コンパイルおよびインストール
```bash
$ make
$ make install
```

GCCと同様に， `.bashrc` に下記のように追記し，パスを登録する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```bash
target=$HOME"/local/openmpi-2.1.1"      # prefixで指定したOpenMPIのフォルダ
$ export PATH=${target}/bin:$PATH
$ export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
```
ログインしなおして以下のコマンドの出力を確認する．  
MPIプログラムのコンパイルコマンドが有効かどうか確認
```bash
$ mpicxx -v
```
出力例： `CXX=[c++ compiler]` で指定したコンパイラのバージョン情報が表示される．  

OpenMPIライブラリのバージョンを確認する
```bash
$ ompi_info
```
たくさん情報が出てくるが，2行目の Open MPI: (version) や，
後ろのほうの各対応APIのバージョンの右の ", Component v (version)" の数字が
インストールしたものとあっているか確認する．


<a id="FFTwのコンパイルとインストール"></a>
<a href="#FFTwのコンパイルとインストール"></a>
## FFTw のコンパイルとインストール
既存のbuildディレクトリを使いまわすのでOpenMPI関係のファイルの消去
```bash
$ rm * -r
```

makeファイルの作成  
参考：http://www.fftw.org/doc/Installation-on-Unix.html  
Fortranを使用しない場合は `--with-g77-wrappers` は書かない．  
GCCの `configure` と同様すべてのオプションを改行せずに続けて入力すること．  
(ここでは `$HOME/local/fftw-3.3.6` にインストールする)  
```bash
$ ../fftw-3.3.6/configure --prefix=$HOME/local/fftw-3.3.6  
--enable-mpi --enable-threads --enable-openmp  
--enable-static --enable-shared  
--enable-sse2 --enable-avx --with-g77-wrappers  
--enable-float --enable-sse  
```
コンパイルおよびインストール
```bash
$ make
$ make install
```
このFFTw (32bit単精度版) を使用する場合，コンパイル時に
```bash
-I$HOME/local/fftw-3.3.6/include
-L$HOME/local/fftw-3.3.6/lib
-lfftw3f_mpi
-lfftw3f
-lm
```
をコンパイルオプションに追加する．  
ここで，`include` フォルダと `lib` フォルダは path を明示的に指定する必要があるが，コンパイル時の `--prefix=***` に指定したディレクトリ以下の構成は通常変わらないので，下記のように環境変数を追加しておくと便利である．
```bash
export FFTW_ROOT=$HOME/local/fftw-3.3.6
```
この環境変数 `FFTW_ROOT` を利用する場合，コンパイルオプションは以下のように書ける
```bash
-I$FFTW_ROOT/include
-L$FFTW_ROOT/lib
-lfftw3f_mpi
-lfftw3f
-lm
```

FDPSのParticleMesh機能拡張は単精度版のfftwを使用するが，ほかの精度のバージョンも使いたければ
configure 設定の最後につけた `--enable-float --enable-sse` を以下のように書き換える．
```
(消去，このオプションを付けない)  : 標準の倍精度(64bit)浮動小数点
--enable-long-double           : 標準の４倍精度(128bit)浮動小数点
--enable-quad-precision        : 非標準 __float128 四倍精度(128bit)浮動小数点
```
これらの詳細は参考HPを参照のこと．
また，異なる精度のFFTwの指定は上記コンパイルオプションの `-lfftw3*` を書き換えることで行う．上記は32bit単精度版の場合で，64bit倍精度版なら `-lfftw3_mpi -lfftw3` となる．


GCC, OpenMPIと同様に `.bashrc` に下記のように追記し，パスを登録する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```bash
target=$HOME"/local/fftw-3.3.6"      # prefixで指定したFFTwのフォルダ名
$ export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
```
FFTwを使用するプログラムをコンパイルしてできた実行ファイル(例えば a.out )にlddコマンドを使ってFFTwが存在しているかどうかを確認できる．
```bash
$ ldd a.out
```


<a id="FDPSのparticle_mesh拡張機能のコンパイル"></a>
<a href="#FDPSのparticle_mesh拡張機能のコンパイル"></a>
## おまけ1：FDPSのparticle_mesh拡張機能のコンパイル
ダウンロードしたFDPSのファイルを適当なフォルダ(ここでは `$HOME/local/FDPS-3.0` )に展開しておく．

Makefileの編集(configureは付属していない)  
上記のように環境を構築しているのなら,
```bash
$ cd $HOME/local/FDPS-3.0/src/particle_mesh
$ emacs Makefile
```
このMakefileのコンパイルオプションとFFTwライブラリの参照先を編集する．

まずデバッグ用の派生版を作る

```bash
CC = mpicxx
CFLAGS = -O0 -DMPICH_IGNORE_CXX_SEEK -g3 -Wall
INCLUDE_FFTW = -I$HOME/local/fftw-3.3.6/include
```

デバッグ版のコンパイル

```bash
$ make
```

`libpm.a` (デバッグ版)が作成される．これを別名で (例えば `libpm_debug.a` 等) 保存する．

```bash
$ mv libpm.a libpm_debug.a
```

Makefile を正式版仕様に再編集

```bash
$ emacs Makefile
CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
```

デバッグ用バージョンの中間生成物を除去

```bash
$ rm -r *.o
```

コンパイル

```bash
$ make
```
`libpm.a` (正式版)が生成される．


<a id="複数のgccバージョンの切り替え"></a>
<a href="#複数のgccバージョンの切り替え"></a>
## おまけ2: 複数のgccバージョンの切り替え
複数の gcc のバージョンを切り替えたい場合，下記のようなスクリプトを作成しておき 自分の `.bashrc` に読み込むと便利である．

```bash
#!/bin/sh
# additional environment setting for FDPS

#--- select GCC & library (write directory name)
#use_gcc="default"
use_gcc="gcc-6.3.0"

use_ompi="openmpi-2.1.1"
use_fftw="fftw-3.3.6"

#--- user install path
usr_env=$HOME"/local"

#------ symbolic link for default environment
default_env=$HOME"/local_path"

#------------------------------------------------------------
#   All settings in above
#------------------------------------------------------------
#------ PATH for GCC (user installed)
target=${usr_env}'/'${use_gcc}
if [ ${use_gcc} = "default" ]; then
  :
else
  #--- add symbolic link
  export PATH=${default_env}:$PATH

  #--- add user install GCC
  if [ -f ${target}/bin/g++ ]; then
    export PATH=${target}/bin:$PATH
    export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${target}/lib64:$LD_LIBRARY_PATH
    export LD_RUN_PATH=${target}/lib:$LD_RUN_PATH
    export LD_RUN_PATH=${target}/lib64:$LD_RUN_PATH
  fi
fi

#--- PATH for OpenMPI
target=${usr_env}'/'${use_ompi}
if [ -f ${target}/bin/mpicxx ]; then
  export PATH=${target}/bin:$PATH
  export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
fi

#--- PATH for FFTw
target=${usr_env}'/'${use_fftw}
if [ -f ${target}/include/fftw3.h ]; then
  export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
  export FFTW_ROOT=${target}
fi
```

上記のスクリプトファイルを `user_opt.sh` という名前で自分のホームディレクトリに保存しておけば，下記のように `.bashrc` に追記することで設定スクリプトを読み込める
```bash
source $HOME/user_opt.sh
```

一番上の部分の `use_gcc` , `use_ompi` , `use_fftw` に自分がインストールした( `--prefix` に指定した)フォルダ名を，
`usr_env` にそれらのインストールフォルダが存在するディレクトリ (上記の例では `$HOME/local` )を入力するだけでパスの設定が完了する．  
また，違うバージョンを後から追加した場合でも，同じインストールディレクトリにインストールしていれば上記変数のフォルダ名のバージョンを新しくインストールしたものに変えれば適用できる．  
gcc を環境デフォルトのものに戻したくなった場合には，スクリプトの `use_gcc="default"` のコメントアウトを戻し，バージョン指定の方をコメントアウトすればよい．



<a id="gtestのインストール"></a>
<a href="#gtestのインストール"></a>
## おまけ3: Google C++ test framework (google test)のインストール
md_fdpsではテストコードの一部に [googel test](https://github.com/google/googletest) を使用している．  
テストコードのコンパイルおよび実行には google test のイントールが必要である．

まず，上記のGitHubリポジトリからソースコードを入手，展開する

```bash
$ pwd
/home/hogehoge/install

$ ls
googletest-master.zip

$ unzip ./googletest-master.zip
```

展開した google test に移動し，コンパイル用のディレクトリを作成する

```bash
$ cd ./gogletest-master

$ mkdir build
```

cmake を使用しコンパイルする．  
ここで，ユーザー権限でインストールしたGCCを使いたい場合，自動的には探してくれないので環境変数で指定する．  
`cmake` コマンド実行時に `make` で使用するコンパイラのバージョンおよび絶対パスが表示されるので，自分が使用する予定のものになっているか確認すること．

```bash
#--- 環境変数でコンパイラを指定
$ export CC=gcc
$ export CXX=g++

$ cd build
$ cmake .
$ make
```

環境変数 `GTEST_ROOT` を登録する．

```bash
$ cd ..
$ pwd
/home/hogehoge/googletest-master

$ export GTEST_ROOT=$HOME/googletest-master
```

あるいは，ユニットテストのコンパイルに必要な以下のファイルを見つけられるように md_fdps の `Makefile` を修正する．

```
googletest-master/googletest/include
googletest-master/build/googlemock/gtest/libgtest.a
```
