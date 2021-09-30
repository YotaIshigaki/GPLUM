!本プログラム内で、粒子および相互作用の定義をしている
!FDPSではメインプログラムが一つのサブルーチンとして定義されているため、グローバル変数を定義したい場合は、下記module内で定義する
!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   implicit none
   integer(c_int) :: np !壁粒子の数　グローバル変数として定義
   integer(c_int), parameter :: nj=24
   real(c_double), parameter :: so=4.0d-1,fri=1.0d0,etan=2.0d2
   real(c_double) :: knn=1.0d4,ktt,etat

   !**** Full particle type 粒子の定義　以下で !$fdps の記号はコメントアウトの記号ではなく,FDPS独自の指示文であることに注意　以降この記号を"指示文"と呼ぶ　粒子の定義についてより理解するためにはチュートリアルp17および仕様書p37以降が参考になる

   type, public, bind(c) :: full_particle !$fdps FP,EPI,EPJ,Force !full_particleという名前のデータ型として,粒子オブジェクトを定義している　指示文以降はこのデータ型が,ユーザーが定義すべき４つのユーザー定義型FP,EPI,EPJ,Forceを担っていることを指示している　４つを別々に定義しても良い　４つの定義型の役割は簡単に述べると, FullParticle型(FP)：粒子データの本体,メインプログラム中からはこのデータ型しか参照できないことに注意　EssentialParticleI型(EPI)：相互作用計算時の力を計算する粒子（力を受ける側の粒子）　EssentialParticleJ型(EPJ)：力を及ぼしている側の粒子　Force型(FP)：EPIに対応して,相互作用の計算結果を格納する型　これらは相互作用計算のコードからより理解できる
      !$fdps copyFromForce full_particle (en(),en()) (esy(),esy()) (esz(),esz()) (je(),je()) (force,force) (moment,moment) !Force型のどの変数をFullParticle型のどの変数にコピーするのかを指示　すなわちどの計算結果をコピーするのか
      !$fdps copyFromFP full_particle (id,id) (mass,mass) (pos,pos) (vel,vel) (smth,smth) (rw,rw) (u,u) (f,f) (r,r) (en(),en()) (esy(),esz()) (esy(),esz()) (je(),je()) (veldum,veldum) (rwdum,rwdum) (fdum,fdum) !FullParticle型からEssentialParticleIおよびJ型にどの変数をコピーするのかを指示　すなわちどの変数が相互作用計算に必要なのか 新しい変数を定義する際はここに加える必要がある
      !$fdps clear  en()=keep, esy()=keep, esz()=keep, je()=keep !Force型の変数について初期化の方法を指示　keepの場合保持,書かなければ0クリア　詳しくは仕様書p43
      !以下粒子が持つ変数　指示文は必須物理量を指定している
      integer(kind=c_long_long) :: id
      real(kind=c_double) :: mass !$fdps charge
      real(kind=c_double) :: pmi
      real(kind=c_double) :: smth !$fdps rsearch
      !探索半径
      type(fdps_f64vec) :: pos !$fdps position
      !type(fdps_f64vec)はFDPS独自の３成分ベクトル型
      type(fdps_f64vec) :: qq
      type(fdps_f64vec) :: vel !$fdps velocity
      type(fdps_f64vec) :: rw
      type(fdps_f64vec) :: u,f
      type(fdps_f64vec) :: force,moment
      type(fdps_f64vec) :: veldum,rwdum,fdum,mdum
      integer(kind=c_long_long) :: je(24)
      real(kind=c_double) :: en(24)
      real(kind=c_double) :: esy(24),esz(24)
      real(kind=c_double) :: r
   end type full_particle

   contains

   !**** Interaction function (particle-particle) 概略は,各プロセスが担当している粒子群がep_i,それらと相互作用しうる粒子群がep_jとして引き渡される　それらが持つ変数を用いてep_iの各粒子が受ける相互作用力を計算し,結果をf(Force型)に格納する　Force型でkeepを指示した変数は前回までの結果が保持されている　各粒子は変数idでしか識別できないことに留意
   subroutine calc_force_pp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp !各粒子数
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j,jj,j0,jk
      integer(c_long_long) :: id_i,id_j
      real(c_double) :: zll,ri,rj,gap,gap_xy
      real(c_double) :: asa,aca,asb,acb
      type(fdps_f64vec) :: xi,xj,ai,gap_vec
      !* Variables of actf
      real(c_double) :: dis,un,usy,usz,veln,velsy,velsz,visn,vissy,vissz,hn,hsy,hsz,hte,en
      real(c_double) :: a11,a12,a13,a21,a22,a23,a24,a25,a26,a31,a32,a33,a34,a35,a36
      type(fdps_f64vec) :: du,dru,dv,drv

      do i=1,n_ip
         id_i = ep_i(i)%id
         if (id_i <= np) cycle !壁粒子が受ける力は計算しない
            xi%x = ep_i(i)%pos%x
            xi%y = ep_i(i)%pos%y
            xi%z = ep_i(i)%pos%z
            ri = ep_i(i)%r
            do j=1,n_jp
               id_j = ep_j(j)%id

               if (id_i == id_j) cycle !同粒子間は計算しない
                  jk = 0
                  do jj=1,nj-6
                     if (ep_i(i)%je(jj) == id_j) then
                        jk = jj
                        exit
                     endif
                  enddo

                  xj%x = ep_j(j)%pos%x
                  xj%y = ep_j(j)%pos%y
                  xj%z = ep_j(j)%pos%z
                  rj = ep_j(j)%r
                  gap_vec = xj - xi
                  gap = dsqrt(gap_vec*gap_vec)


                  if (gap <= (ri+rj)) then
                     gap_xy = dsqrt(gap_vec%x*gap_vec%x+gap_vec%y*gap_vec%y)
                     if (gap_xy == 0.0d0) then
                        if (gap_vec%z <= 0.0d0) then
                           asa = 0.0d0
                           aca = 1.0d0
                        else
                           asa = 0.0d0
                           aca = -1.0d0
                        endif
                     else
                        asa = gap_vec%y/gap_xy
                        aca = gap_vec%x/gap_xy
                     endif
                     asb = -gap_vec%z/gap
                     acb = gap_xy/gap
                     if (jk == 0) then
                        do jj=1,nj-6
                           if (ep_i(i)%je(jj) == 0) then
                              jk = jj
                              f(i)%je(jj) = id_j
                              exit
                           endif
                        enddo
                     endif

                     !************************ Part of actf ****************************
                     !粒子群をサブルーチンに引き渡す方法が不明なためactfを作っていない

                     ktt = knn*so
                     etat = etan*so

                     dis = ri + rj - gap

                     a11=acb*aca
                     a12=asa*acb
                     a13=-asb
                     a21=-asa
                     a22=aca
                     a23=0.0d0
                     a24=asb*aca
                     a25=asa*asb
                     a26=acb
                     a31=asb*aca
                     a32=asa*asb
                     a33=acb
                     a34=asa
                     a35=-aca
                     a36=0.0d0

                     du = ep_i(i)%u - ep_j(j)%u
                     dru = ri*ep_i(i)%f + rj*ep_j(j)%f

                     dv = ep_i(i)%vel - ep_j(j)%vel
                     drv = ri*ep_i(i)%rw + rj*ep_j(j)%rw
                     un = du%x*a11+du%y*a12+du%z*a13
                     usy = du%x*a21+du%y*a22+du%z*a23+dru%x*a24+dru%y*a25+dru%z*a26
                     usz = du%x*a31+du%y*a32+du%z*a33+dru%x*a34+dru%y*a35+dru%z*a36

                     veln  = dv%x*a11+dv%y*a12+dv%z*a13
                     velsy = dv%x*a21+dv%y*a22+dv%z*a23+drv%x*a24+drv%y*a25+drv%z*a26
                     velsz = dv%x*a31+dv%y*a32+dv%z*a33+drv%x*a34+drv%y*a35+drv%z*a36

                     !---------------------------------------
                     if (ep_i(i)%en(jk) == 0.0d0) then
                        if (un /= 0.0d0) then
                           usy = usy*dis/dabs(un)
                           usz = usz*dis/dabs(un)
                        endif
                        un = dis
                     endif
                     !---------------------------------------
                     ! == VISCOUS FORCE ==
                     !  == linear ==
                     !  ** NORMAL **
                     visn = etan*veln
                     !  ** TANGENTIAL **
                     vissy = etat*velsy
                     vissz = etat*velsz

                     ! == ELASTIC FORCE ==
                     !  ** NORMAL **
                     !  == linear ==
                     !f(i)%en(jk) = f(i)%en(jk) + knn*un
                     f(i)%en(jk) = knn*dis
                     !  == non-linear ==
                     ! ***** this is the most important change *****
                     !      f(i)%en(jk) = knn*(dis**1.5d0)
                     !      visn = etan*veln*dsqrt(dis)
                     ! *********************************************
                     !  ** TANGENTIAL **
                     !  == linear ==
                     f(i)%esy(jk) = f(i)%esy(jk) + ktt*usy !Force型に積算
                     f(i)%esz(jk) = f(i)%esz(jk) + ktt*usz
                     !  == non-linear =
                     !     ** overlap **
                     !      nonusy = dabs(f(i)%esy(jk)/ktt)**(2.0d0/3.0d0)
                     !      nonusz = dabs(f(i)%esz(jk)/ktt)**(2.0d0/3.0d0)
                     ! == RESULTANT FORCE ==
                     hn = f(i)%en(jk) +visn

                     if (hn >= 0.0d0) then
                        hsy = f(i)%esy(jk) + vissy
                        hsz = f(i)%esz(jk) + vissz
                     else
                        f(i)%esy(jk) = 0.0d0
                        f(i)%esz(jk) = 0.0d0
                        hn = 0.0d0
                        hsy = 0.0d0
                        hsz = 0.0d0
                     endif

                     ! == COULOMB FRICTION on ELASTIC PART ==
                     hte = hsy*hsy + hsz*hsz
                     hte = dsqrt(hte)
                     if (hte > fri*hn) then
                        hsy = fri*hn*hsy/hte
                        hsz = fri*hn*hsz/hte

                        f(i)%esy(jk) = (hsy - vissy) !/ ktt !kttの必要性は検討の余地あり
                        f(i)%esz(jk) = (hsz - vissz) !/ ktt
                     endif

                     ! == SUM UP Contact forces in Cartecian coordinates ==
                     !計算結果をForce型に格納
                     f(i)%force%x = f(i)%force%x + (-hn*a11-hsy*a21-hsz*a31)
                     f(i)%force%y = f(i)%force%y + (-hn*a12-hsy*a22-hsz*a32)
                     f(i)%force%z = f(i)%force%z + (-hn*a13-hsy*a23-hsz*a33)
                     f(i)%moment%x = f(i)%moment%x + ri*(-hsy*a24-hsz*a34)
                     f(i)%moment%y = f(i)%moment%y + ri*(-hsy*a25-hsz*a35)
                     f(i)%moment%z = f(i)%moment%z + ri*(-hsy*a26-hsz*a36)
!************************ End of actf ****************************
                  else
                     if (jk /= 0) then !接触が切れた場合保持していた接触力をクリア
                        f(i)%en(jk) = 0.0d0
                        f(i)%esy(jk) = 0.0d0
                        f(i)%esz(jk) = 0.0d0
                        f(i)%je(jk) = 0
                     endif
                  endif
             end do
      end do

    end subroutine calc_force_pp

end module user_defined_types
