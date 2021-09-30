### note for ./script

md_fdps の支援スクリプトについて

#### VMDmovie_convert.py
md_fdps が出力した　`./pdb/vmd***.pdb` ,`./pos/pos***.tsv` ファイルからVMDで可視化するためのファイル `./vmd_movie.pdb` , `./vmd_movie.crd` を生成する．  
time step の範囲，ステップ間隔は自動で取得するが，実行時オプションでも指定可能．詳細は `-h` オプションのヘルプを参照．  


#### param_file_comment.py
`./model/` フォルダにある `***.param` ファイル内のコメントを一括で更新する．  
スクリプト冒頭で定義されている `comment[]` ディクショナリにキーをキーワード，値をコメントとして定義すると，キーワードの直前のコメントを　`comment[]` に定義されたものに置換する．  
コメントとみなす行は `comment_mark` で指定する．  
行頭に `comment_mark` を含む行が連続する場合，その部分は一塊のコメントとみなされる．

#### convert_model_indicator.py
`./model/` フォルダにある `***.mol2` ファイルを参照してモデル，原子の種類を管理する enum class の定義ファイル `./src/enum_model.hpp` を生成する．  
付属のmMakefileを用いて `make` する時に自動で実行されるが，新しくモデルを追加した際にはもう一度手動で実行するか， `make clean_all` を呼んでから `make` する必要がある．　　
