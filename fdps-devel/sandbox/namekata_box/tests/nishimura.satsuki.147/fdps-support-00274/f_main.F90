!以下コメントにおける"仕様書"とはdoc_specs_ftn_ja.pdfを指す
!マルチプロセス並列で実行した場合(MPIを使用mpirunで実行)、指定しただけのプロセスが立ち上がり、各プロセスが以下のすべての命令を独立に実行することに注意　その際各プロセスに自動的に割り振られる番号をランク数と呼ぶ
!粒子オブジェクトの定義、相互作用計算の関数はuser_defined.F90で定義されている

!-----------------------------------------------------------------------
!/////////////////////// < M A I N  R O U T I N E > ////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
  use fdps_module
  use user_defined_types
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
  use OMP_lib !Including library of OpenMP to use omp_get_wtime()
#endif
  implicit none
  !* Local parameters
  integer, parameter :: ntot=10000 !The number of particles
  double precision, parameter :: rmax=0.5d0,rmin=0.4d0 !Radius of particles
  !-(force parameters)
  real, parameter :: theta = 0.5 ! 見込み角　ツリー法を用いる際必要となる　詳しくは　http://jun.artcompsci.org/kougi/keisan_tenmongakuII/note10/node3.html　等を参照　しかしDEMでは重要でないと考える
  integer, parameter :: n_leaf_limit = 8 !ツリーを切るのをやめる粒子数の上限
  integer, parameter :: n_group_limit = 64 !相互作用リストを共有する粒子数の上限
  !-(domain decomposition)
  real, parameter :: coef_ema= 0.5 !平滑化係数　マルチプロセスで実行した場合の領域分割のパラメータ　０から１の間で設定　値が大きいほど最新の粒子分布の情報が分割に反映される　詳しくは仕様書を参照
  !-(timing parameters)
 !double precision, parameter :: time_end = 10.0d0 !全時間 (original)
  double precision, parameter :: time_end = 1.0d0 !全時間
  double precision, parameter :: dt = 1.0d0/1280.0d0 !時間刻み
  double precision, parameter :: dt_diag = 1.0d0/8.0d0 !標準出力への出力間隔
  double precision, parameter :: dt_snap = 1.0d-1 !ファイルへの出力間隔

  !* Local variables
  integer :: i,j,k,num_loop
  integer :: psys_num,dinfo_num,tree_num
  integer :: nloc,myrank
  logical :: clear
  double precision :: ekin,erot,etot,e_tot(int(time_end/dt_diag))
  double precision :: time_diag,time_snap,time_sys
  type(fdps_f64vec) :: pos_ll,pos_ul
  type(fdps_controller) :: fdps_ctrl
  type(full_particle), dimension(:), pointer :: ptcl
  type(c_funptr) :: pfunc_ep_ep
  integer, dimension(1) :: tmp
  !-(IO)
  character(len=64) :: fname
  !-(time count)
  real :: time1,time2
  real(kind=c_double) :: start_time_whole, end_time_whole
  real(kind=c_double) :: start_time,end_time,etime


  !* Initialize FDPS
  call fdps_ctrl%PS_Initialize()

  if(fdps_ctrl%get_rank() == 0) then !fdps_ctrl%get_rank()：各プロセスのランク数を取得　ランク数が0のプロセスのみ以下を実行
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
     time1 = omp_get_wtime() !計算時間を求めるためにマシンタイムを取得　マルチスレッド並列(OpenMP)を使用する場合omp_get_wtime()で取得する必要がある
#endif
     print * ,'ntot=',ntot
     fname =trim("result") // "/"// "eng_dat.dat"      
     open(unit=10,file=trim(fname),action='write',status='replace') !エネルギー出力用ファイルを開く
  endif
  call fdps_ctrl%barrier()
  start_time_whole = fdps_ctrl%get_wtime()

  !* Set the box size 領域のサイズを指定
  pos_ll%x = -10.0d0
  pos_ll%y = -10.0d0
  pos_ll%z = 0.0d0
  pos_ul%x = 10.0d0
  pos_ul%y = 10.0d0
  pos_ul%z = 10.0d0

  !* Create domain info object 
  call fdps_ctrl%create_dinfo(dinfo_num) !領域情報オブジェクトを作り、識別番号をdinfo_numに格納
  call fdps_ctrl%init_dinfo(dinfo_num,coef_ema) !領域情報オブジェクトを初期化
  call fdps_ctrl%set_boundary_condition(dinfo_num,fdps_bc_periodic_xy) !境界条件を設定する　ここではx,y方向を周期境界、z方向を開放としている　詳しくは仕様書
  call fdps_ctrl%set_pos_root_domain(dinfo_num,pos_ll,pos_ul) !領域の下限と上限を設定する 途中で再びこの命令を呼びだせば、領域のサイズを変更することもできる

  !* Create particle system object
  call fdps_ctrl%create_psys(psys_num,'full_particle') !full_particle型の粒子群オブジェクト(粒子データの集合)をつくり、その識別番号をpys_numに格納する
  call fdps_ctrl%init_psys(psys_num) !粒子群オブジェクトの初期化

  !* Create tree object
  call fdps_ctrl%create_tree(tree_num, &
       "Short,full_particle,full_particle,full_particle,Symmetry") !ツリーオブジェクトを作成し、識別番号をtree_numに格納　
  call fdps_ctrl%init_tree(tree_num,ntot,theta, &
       n_leaf_limit,n_group_limit) !ツリーオブジェクトの初期化

  !* Make an initial condition
  call setup_IC(fdps_ctrl,psys_num,pos_ll,pos_ul,ntot,rmax,rmin) !粒子の初期配置を決めるサブルーチン

  tmp(1) = np !npは壁粒子の数　user_dfined.F90内で定義されている　ランク0のプロセスのみサブルーチンsetup_ICの中でnpを数える　FDPSのバグのため一度配列tmpに代入する必要がある
  call fdps_ctrl%broadcast(tmp,1,0) !ランク0のプロセスがtmp(=np)の値を他のプロセスへ伝達している　詳しくは仕様書　2017年の時点ではバグのため配列データしか伝達できない
  np = tmp(1) !これにより全プロセスがnpの値を共有する 


  !* Adjust the positions of the particles that run over
  !  the computational boundaries.
  call fdps_ctrl%adjust_pos_into_root_domain(psys_num,dinfo_num) !周期境界の場合、領域からはみ出した粒子を適切にもどす   

  !* Domain decomposition and exchange particle
  call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num) !計算領域の分割
  call fdps_ctrl%exchange_particle(psys_num,dinfo_num) !粒子が適切なドメインに配置されるように、粒子の交換を行う

  !* Compute force at the initial time
#if 0
  call fdps_ctrl%barrier()
  start_time = fdps_ctrl%get_wtime()
#endif
  pfunc_ep_ep = c_funloc(calc_force_pp) !user_defined.F90で定義される関数calc_force_ppの関数ポインタ（のアドレス）を代入
  call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
       pfunc_ep_ep, &
       psys_num,    &
       dinfo_num) !粒子群オブジェクトの粒子すべての相互作用を計算し、結果を粒子群オブジェクトに書き込む
#if 0
  call fdps_ctrl%barrier()
  end_time = fdps_ctrl%get_wtime()
  etime = end_time - start_time
  if (fdps_ctrl%get_rank() == 0) then
    write(*,*)'etime      = ',etime
    write(*,*)'start_time = ',start_time
    write(*,*)'end_time   = ',end_time
  end if
  call fdps_ctrl%ps_finalize()
  stop 0
#endif

  !* Time integration
  time_diag = 0.0d0
  time_snap = 0.0d0
  time_sys  = 0.0d0
  num_loop = 0
  do 
     if ( (time_sys >= time_snap) .or. &
          (((time_sys + dt) - time_snap) > (time_snap - time_sys)) ) then
        call output(fdps_ctrl,psys_num) !出力ファイルへ書き出すサブルーチン
        time_snap = time_snap + dt_snap
     end if

     !* Compute energies and output the results
     clear = .true.
     call calc_energy(fdps_ctrl,psys_num,etot,ekin,erot,clear) !全運動エネルギーを計算し、ファイルへ出力するサブルーチン
     if (fdps_ctrl%get_rank() == 0) then !ランク0のプロセスのみ以下を実行
        if ( (time_sys >= time_diag) .or. &
             (((time_sys + dt) - time_diag) > (time_diag - time_sys)) ) then
           write(*,100)time_sys,etot
100        format("time: ",1es20.10e3,", energy: ",1es20.10e3)
           time_diag = time_diag + dt_diag           
           write(10,'(f10.2,f18.12)')time_sys,etot
        end if
     end if

     !* Time evolution
     call nposit(fdps_ctrl,psys_num,dt) !粒子の時間発展を計算するサブルーチン
     time_sys = time_sys + dt


     !* Adjust the positions of the particles that run over
     !  the computational boundaries.
     call fdps_ctrl%adjust_pos_into_root_domain(psys_num,dinfo_num) !周期境界の場合、領域からはみ出した粒子を適切にもどす  

     !* Domain decomposition & exchange particle
     call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num) !計算領域の分割
     call fdps_ctrl%exchange_particle(psys_num,dinfo_num) !粒子が適切なドメインに配置されるように、粒子の交換を行う

     !* Force calculation
     pfunc_ep_ep = c_funloc(calc_force_pp) !user_defined.F90で定義される関数calc_force_ppの関数ポインタ（のアドレス）を代入
     call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
          pfunc_ep_ep, &
          psys_num,    &
          dinfo_num) !粒子群オブジェクトの粒子すべての相互作用を計算し、結果を粒子群オブジェクトに書き込む
   

     !* Update num_loop
     num_loop = num_loop + 1

     !* Termination
     if (time_sys >= time_end) then
        exit
     end if
  end do

  call fdps_ctrl%barrier()
  end_time_whole = fdps_ctrl%get_wtime()
  if(fdps_ctrl%get_rank() == 0) then
     close(unit=10) !エネルギー出力用ファイルを閉じる
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
     time2 = omp_get_wtime() !計算時間を求めるためにマシンタイムを取得
     print *,'time = ', time2-time1, '[s]' !計算時間を出力
#endif
     write(*,*)"time = ",end_time_whole - start_time_whole," [s] (get_wtime)"
  endif
  print*, np

  !* Finalize FDPS
  call fdps_ctrl%PS_Finalize()
end subroutine f_main

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!///////////////////////// < S E T U P _ I C > /////////////////////////
!-----------------------------------------------------------------------
subroutine setup_IC(fdps_ctrl,psys_num,pos_ll,pos_ul,nptcl_glb,rmax,rmin)
  use fdps_vector
  use fdps_module
  use user_defined_types
  implicit none
  type(fdps_controller), intent(IN) :: fdps_ctrl
  integer, intent(IN) :: psys_num,nptcl_glb
  double precision, intent(IN) :: rmax,rmin
  !   integer, intent(out):: np
  !* Local parameters
  double precision, parameter :: m_tot=1.0d0,mass_wall = 4.0d-1,height=5.0d0
  double precision :: pi,search_radius
  !* Local variables
  integer :: i,j,k,n=0
  integer :: nprocs,myrank,np_f
  integer :: ipx,ipy,count=0
  double precision :: r2,cm_mass,rn,delr,r,mass,px,py,pz
  type(fdps_f64vec) :: pos_ll,pos_ul
  type(fdps_f64vec) :: cm_pos,cm_vel,pos
  type(full_particle), dimension(:), pointer :: ptcl
  character(len=64) :: fname

  pi = datan(1.0d0)*4.0d0
  !* Get # of MPI processes and rank number
  nprocs = fdps_ctrl%get_num_procs() !全プロセス数を取得
  myrank = fdps_ctrl%get_rank() !自分のプロセス番号を取得

  !* Make an initial condition at RANK 0
  if (myrank == 0) then !ランク0のプロセスのみが、全粒子の初期状態を設定する

     !fname =trim("result") // "/"// "r.dat" !粒子半径出力用ファイル
     !open(unit=11,file=trim(fname),action='write',status='replace')      

     !* Set # of local particles
     call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_glb) !各プロセスがもつ粒子数を設定　ここではnptcl_glb(=ntot)の全粒子数 

     !* Create an uniform sphere of particles
     !** get the pointer to full particle data
     call fdps_ctrl%get_psys_fptr(psys_num,ptcl) !識別番号pys_numの粒子群オブジェクトで管理される粒子配列(Full_particle型)へのポインタをptclへ格納

     !** initialize Mersenne twister
     call fdps_ctrl%MT_init_genrand(0) !乱数生成器の初期化

     !* Set of particles at bottom surface 底面粒子の設定
     delr = 1.0d-4
     rn = rmax+delr
     ipx = dint((pos_ul%x-pos_ll%x)/(2.0d0*rn))+1
     ipy = dint((pos_ul%y-pos_ll%y)/(2.0d0*rn))+1
     np = ipx*ipy
     do i=1,ipy
        do j=1,ipx
           n=n+1 
           ptcl(n)%id   = n 
           ptcl(n)%mass = mass_wall
           ptcl(n)%pos%x = pos_ll%x+2.0d0*rn*dble(j-1)+1.0d-4
           ptcl(n)%pos%y = pos_ll%y+2.0d0*rn*dble(i-1)+1.0d-4
           ptcl(n)%pos%z = pos_ll%z
           r = rmin+(rmax-rmin)*fdps_ctrl%MT_genrand_res53() ![0.0:1.0)間の浮動小数点を生成するインターフェイス　詳しくは仕様書
           ptcl(n)%r = r
           !write(11,fmt='(f10.5,a)') r,","         
           ptcl(n)%pmi = mass_wall*2.0d0*r*r/5.0d0
           ptcl(n)%smth = 0.0d0
           ptcl(n)%qq = 0.0d0
           ptcl(n)%rw = 0.0d0
           ptcl(n)%u = 0.0d0
           ptcl(n)%f = 0.0d0
           ptcl(n)%force = 0.0d0
           ptcl(n)%moment = 0.0d0
        end do
     end do
     print * ,'np=',np !底面粒子数

     !Set of free particles 自由粒子の設定
     px = pos_ll%x + rn*fdps_ctrl%MT_genrand_res53() 
     py = pos_ll%y + rn*fdps_ctrl%MT_genrand_res53() 
     pz = pos_ll%z + height
     do i=np+1,nptcl_glb
        ptcl(i)%id   = i
        ptcl(i)%pos%x = px
        px = px + 3.0d0*rn + 0.5d0*rn*fdps_ctrl%MT_genrand_res53() 
        ptcl(i)%pos%y = py + 0.5d0*rn*fdps_ctrl%MT_genrand_res53() 
        if (px >= pos_ul%x-2.0d0*rn) then !xが領域を超える時、yをずらしxをもどす
           py = py + 3.0d0*rn + rn*fdps_ctrl%MT_genrand_res53()
           px = pos_ll%x + rn*fdps_ctrl%MT_genrand_res53() 
        endif
        ptcl(i)%pos%z = pz + 0.5d0*rn*fdps_ctrl%MT_genrand_res53()
        if (py >= pos_ul%y-2.0d0*rn) then
           pz = pz + 3.0d0*rn + rn*fdps_ctrl%MT_genrand_res53()
           py = pos_ll%y + rn*fdps_ctrl%MT_genrand_res53() 
        endif
        r = rmin+(rmax-rmin)*fdps_ctrl%MT_genrand_res53()
        ptcl(i)%r = r
        !write(11,fmt='(f10.5,a)') r,","         
        mass = r*r*r*pi*4.0d0/3.0d0
        ptcl(i)%mass = mass
        ptcl(i)%pmi = mass*2.0d0*r*r/5.0d0
        ptcl(i)%smth = r + rmax + 1.0d-4
        ptcl(i)%vel = 0.0d0
        ptcl(i)%qq = 0.0d0
        ptcl(i)%rw = 0.0d0
        ptcl(i)%u = 0.0d0
        ptcl(i)%f = 0.0d0
        ptcl(i)%force = 0.0d0
        ptcl(i)%moment = 0.0d0
        ptcl(i)%veldum = 0.0d0
        ptcl(i)%rwdum = 0.0d0
        ptcl(i)%fdum = 0.0d0
        ptcl(i)%mdum = 0.0d0
        ptcl(i)%je = 0
        ptcl(i)%en = 0.0d0
        ptcl(i)%esy = 0.0d0
        ptcl(i)%esz = 0.0d0
     end do

     !* Output
     !fname = 'initial.dat'
     !open(unit=9,file=trim(fname),action='write',status='replace', &
     !     form='unformatted',access='stream')
     !open(unit=9,file=trim(fname),action='write',status='replace')
     !   do i=1,nptcl_glb
     !     !write(9)ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
     !      write(9,'(3es25.16e3)')ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
     !   end do
     !close(unit=9)

     !* Release the pointer
     nullify( ptcl ) !ポインタの開放
     close(unit=11)

  else
     call fdps_ctrl%set_nptcl_loc(psys_num,0)
  end if

end subroutine setup_IC

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////   < N P O S I T >   /////////////////////////
!-----------------------------------------------------------------------
subroutine nposit(fdps_ctrl,psys_num,dt)
  use fdps_vector
  use fdps_module
  use user_defined_types
  implicit none

  type(fdps_controller), intent(IN) :: fdps_ctrl
  integer, intent(IN) :: psys_num
  double precision, intent(IN) :: dt
  !* Local variables
  integer :: i,nptcl_loc
  double precision :: mass,pmi
  type(full_particle), dimension(:), pointer :: ptcl
  type(fdps_f64vec) :: vel,rw,force,moment,g,u,f

  !重力加速度の設定
  g%x=0.0d0
  g%y=0.0d0
  g%z=-9.8d0

  !* Get # of local particles
  nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num) !各プロセスが持つ粒子数を格納

  !* Get the pointer to full particle data
  call fdps_ctrl%get_psys_fptr(psys_num,ptcl) !識別番号pys_numの粒子群オブジェクトで管理される粒子配列(Full_particle型)へのポインタをptclへ格納
  do i=1,nptcl_loc
     if (ptcl(i)%id <= np) cycle !底面粒子の場合はループから抜ける
     mass = ptcl(i)%mass
     pmi = ptcl(i)%pmi
     force = ptcl(i)%force
     moment = ptcl(i)%moment

     vel = ptcl(i)%vel + dt*((force/mass+g)*1.5d0-ptcl(i)%fdum*0.5d0)
     rw = ptcl(i)%rw + dt*(moment/pmi*1.5d0-ptcl(i)%mdum*0.5d0)
     ptcl(i)%vel = vel 
     ptcl(i)%rw = rw

     u = dt * (vel + ptcl(i)%veldum) * 0.5d0
     f = dt * (rw + ptcl(i)%rwdum) * 0.5d0
     ptcl(i)%u = u
     ptcl(i)%f = f

     ptcl(i)%pos = ptcl(i)%pos + u
     ptcl(i)%qq = ptcl(i)%qq + f

     ptcl(i)%veldum = vel
     ptcl(i)%fdum = force / mass + g
     ptcl(i)%mdum = moment / pmi
     ptcl(i)%rwdum = rw
  enddo

  nullify(ptcl) !ポインタの開放

end subroutine nposit

!-----------------------------------------------------------------------
!//////////////////////    S U B R O U T I N E    //////////////////////
!////////////////////// < C A L C _ E N E R G Y > //////////////////////
!-----------------------------------------------------------------------
subroutine calc_energy(fdps_ctrl,psys_num,etot,ekin,erot,clear)
  use fdps_vector
  use fdps_module
  use user_defined_types
  implicit none
  type(fdps_controller), intent(IN) :: fdps_ctrl
  integer, intent(IN) :: psys_num
  double precision, intent(INOUT) :: etot,ekin,erot
  logical, intent(IN) :: clear
  !* Local variables
  integer :: i,nptcl_loc
  double precision :: etot_loc,ekin_loc,erot_loc
  type(full_particle), dimension(:), pointer :: ptcl

  !* Clear energies
  if (clear .eqv. .true.) then
     etot = 0.0d0
     ekin = 0.0d0
     erot = 0.0d0
  end if

  !* Get # of local particles
  nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num) !各プロセスがもつ粒子数を格納
  call fdps_ctrl%get_psys_fptr(psys_num,ptcl) !識別番号pys_numの粒子群オブジェクトで管理される粒子配列(Full_particle型)へのポインタをptclへ格納

  !* Compute energies
  ekin_loc = 0.0d0
  erot_loc = 0.0d0
  do i=1,nptcl_loc
     if (ptcl(i)%id <= np) cycle
     ekin_loc = ekin_loc + ptcl(i)%mass * (ptcl(i)%vel * ptcl(i)%vel) !各プロセスがもつ粒子の重心の運動エネルギー*2の総和
     erot_loc = erot_loc + ptcl(i)%pmi * (ptcl(i)%rw * ptcl(i)%rw) !各プロセスがもつ粒子の重心周りの回転の運動エネルギー*2の総和
  end do
  ekin_loc = ekin_loc * 0.5d0
  erot_loc = erot_loc * 0.5d0
  etot_loc = ekin_loc + erot_loc
  !call fdps_ctrl%get_sum(ekin_loc,ekin)
  !call fdps_ctrl%get_sum(erot_loc,erot)
  call fdps_ctrl%get_sum(etot_loc,etot) !各プロセスのetot_locを全プロセスで足しあわせ結果をetotに格納

  !* Release the pointer
  nullify(ptcl)

end subroutine calc_energy

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////  < O U T P U T >  /////////////////////////
!-----------------------------------------------------------------------
subroutine output(fdps_ctrl,psys_num)
  use fdps_vector
  use fdps_module
  use user_defined_types
  implicit none
  type(fdps_controller), intent(IN) :: fdps_ctrl
  integer, intent(IN) :: psys_num
  !* Local parameters
  character(len=16), parameter :: root_dir="result"
  character(len=16), parameter :: file_prefix_1st="snap"
  character(len=16), parameter :: file_prefix_2nd="proc"
  !* Local variables
  integer :: i,nptcl_loc
  integer :: myrank
  character(len=5) :: file_num,proc_num
  character(len=64) :: cmd,sub_dir,fname
  type(full_particle), dimension(:), pointer :: ptcl
  !* Static variables
  integer, save :: snap_num=0

  !* Get the rank number
  myrank = fdps_ctrl%get_rank() !自分のランク数を格納

  !* Get # of local particles
  nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num) !各プロセスがもつ粒子数を格納

  !* Get the pointer to full particle data
  call fdps_ctrl%get_psys_fptr(psys_num,ptcl) !識別番号pys_numの粒子群オブジェクトで管理される粒子配列(Full_particle型)へのポインタをptclへ格納

  !* Output
  write(file_num,"(i5.5)")snap_num
  write(proc_num,"(i5.5)")myrank
  fname =  trim(root_dir) // "/" &
       // trim(file_prefix_1st) // file_num // "-" &
       // trim(file_prefix_2nd) // proc_num // ".dat"
  open(unit=9,file=trim(fname),action='write',status='replace')
  do i=1,nptcl_loc
     write(9,'(4(f10.5,a))') ptcl(i)%pos%x,",",ptcl(i)%pos%y,",",ptcl(i)%pos%z,",",ptcl(i)%r,"," !粒子配置と半径を出力　povrayでの動画作成のためカンマ区切りにしている
  end do
  close(unit=9)
  nullify(ptcl) !ポインタの開放

  !* Update snap_num
  snap_num = snap_num + 1

end subroutine output
