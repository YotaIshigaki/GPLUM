#ifndef MAX_PAIR
#define MAX_PAIR 5000
#endif

module shift_function
  real(8), save:: r1co0, r1co1, r1co2
  real(8), save:: r6co0, r6co1, r6co2
  real(8), save:: r12co0, r12co1, r12co2
  real(8), save:: r3co1, r3co2
  real(8), save:: r8co1, r8co2
  real(8), save:: r14co1, r14co2
end module shift_function

module ewald_table
  real(8) ::eft(0:4,0:63) = reshape((/ &
  9.985445d-01,-4.343182d-03,-4.271433d-03,-1.319000d-03,6.558229d-05, &
  9.635154d-01,-3.446619d-02,-9.872624d-03,-4.268090d-04,1.415692d-04, &
  8.539463d-01,-7.455800d-02,-9.086888d-03,6.341242d-04,1.076347d-04, &
  6.750769d-01,-1.004364d-01,-3.363002d-03,1.119045d-03,1.082046d-05, &

  4.696700d-01,-1.007025d-01,2.971981d-03,8.750839d-04,-6.299924d-05, &
  2.860746d-01,-8.052218d-02,6.519390d-03,2.915723d-04,-7.324796d-05, &
  1.523237d-01,-5.312552d-02,6.704382d-03,-1.775913d-04,-4.064064d-05, &
  7.090837d-02,-2.948451d-02,4.944704d-03,-3.520935d-04,-4.925659d-06, &

  2.887827d-02,-1.393189d-02,2.881125d-03,-3.075505d-04,1.293026d-05, &
  1.029832d-02,-5.649658d-03,1.379933d-03,-1.904582d-04,1.454575d-05, &
  3.218548d-03,-1.977219d-03,5.546694d-04,-9.263555d-05,9.555838d-06, &
  8.822612d-04,-5.996074d-04,1.894468d-04,-3.688556d-05,4.657932d-06, &

  2.122676d-04,-1.580432d-04,5.543035d-05,-1.228625d-05,1.809561d-06, &
  4.485241d-05,-3.629039d-05,1.397245d-05,-3.468160d-06,5.788463d-07, &
  8.327851d-06,-7.272882d-06,3.046845d-06,-8.367508d-07,1.552541d-07, &
  1.359326d-06,-1.273949d-06,5.765402d-07,-1.735799d-07,3.531452d-08, &

  1.951320d-07,-1.952719d-07,9.489760d-08,-3.109607d-08,6.865060d-09, &
  2.464312d-08,-2.621731d-08,1.361303d-08,-4.826736d-09,1.146881d-09, &
  2.738759d-09,-3.085637d-09,1.704499d-09,-6.508275d-10,1.653371d-10, &
  2.679258d-10,-3.185664d-10,1.865213d-10,-7.639022d-11,2.063405d-11, &

  2.307683d-11,-2.886677d-11,1.785685d-11,-7.818017d-12,2.234888d-12, &
  1.750349d-12,-2.296927d-12,1.496952d-12,-6.986226d-13,2.105081d-13, &
  1.169327d-13,-1.605550d-13,1.099669d-13,-5.457297d-14,1.727229d-14, &
  6.881409d-15,-9.862330d-15,7.083435d-15,-3.730147d-15,1.236245d-15, &

  3.567869d-16,-5.325301d-16,4.003069d-16,-2.232801d-16,7.727578d-17, &
  1.629986d-17,-2.528305d-17,1.985714d-17,-1.171283d-17,4.222811d-18, &
  6.562212d-19,-1.055678d-18,8.649598d-19,-5.388061d-19,2.019083d-19, &
  2.328365d-20,-3.877326d-20,3.309706d-20,-2.174698d-20,8.453291d-21, &

  7.281571d-22,-1.252859d-21,1.112852d-21,-7.704944d-22,3.100994d-22, &
  2.007278d-23,-3.562055d-23,3.288993d-23,-2.397333d-23,9.973076d-24, &
  4.877869d-25,-8.912044d-25,8.546264d-25,-6.552985d-25,2.813401d-25, &
  1.045014d-26,-1.962335d-26,1.952878d-26,-1.574158d-26,6.964744d-27, &

  1.973829d-28,-3.802951d-28,3.925061d-28,-3.324199d-28,1.513650d-28, &
  3.287130d-30,-6.487015d-30,6.940134d-30,-6.172672d-30,2.889018d-30, &
  4.826878d-32,-9.740090d-32,1.079715d-31,-1.008123d-31,4.844197d-32, &
  6.249974d-34,-1.287315d-33,1.478196d-33,-1.448455d-33,7.137873d-34, &

  7.136270d-36,-1.497661d-35,1.781117d-35,-1.831197d-35,9.245026d-36, &
  7.185596d-38,-1.533712d-37,1.889033d-37,-2.037436d-37,1.052801d-37, &
  6.380708d-40,-1.382494d-39,1.763671d-39,-1.995379d-39,1.054338d-39, &
  4.996941d-42,-1.096860d-41,1.449657d-41,-1.720387d-41,9.287501d-42, &

  3.451300d-44,-7.659106d-44,1.049102d-43,-1.306014d-43,7.197545d-44, &
  2.102409d-46,-4.706565d-46,6.685047d-46,-8.730674d-46,4.908083d-46, &
  1.129586d-48,-2.544947d-48,3.751041d-48,-5.140175d-48,2.945444d-48, &
  5.353041d-51,-1.210708d-50,1.853454d-50,-2.665550d-50,1.555841d-50, &

  2.237541d-53,-5.066496d-53,8.065160d-53,-1.217642d-52,7.234606d-53, &
  8.249747d-56,-1.864615d-55,3.090724d-55,-4.900240d-55,2.961794d-55, &
  2.682985d-58,-6.033484d-58,1.043122d-57,-1.737477d-57,1.067669d-57, &
  7.696843d-61,-1.715947d-60,3.100594d-60,-5.428264d-60,3.389284d-60, &

  1.947743d-63,-4.287698d-63,8.116978d-63,-1.494432d-62,9.475746d-63, &
  4.347934d-66,-9.408446d-66,1.871475d-65,-3.625736d-65,2.333431d-65, &
  8.561967d-69,-1.811854d-68,3.800239d-68,-7.752652d-68,5.061641d-68, &
  1.487340d-71,-3.059935d-71,6.796229d-71,-1.461052d-70,9.672512d-71, &

  2.279287d-74,-4.527611d-74,1.070395d-73,-2.426991d-73,1.628445d-73, &
  3.081388d-77,-5.862150d-77,1.484649d-76,-3.553709d-76,2.415601d-76, &
  3.675002d-80,-6.630835d-80,1.813372d-79,-4.587004d-79,3.157371d-79, &
  3.866678d-83,-6.538022d-83,1.950326d-82,-5.219528d-82,3.636664d-82, &

  3.589160d-86,-5.602171d-86,1.846942d-85,-5.236105d-85,3.691346d-85, &
  2.939184d-89,-4.152952d-89,1.539880d-88,-4.631046d-88,3.302141d-88, &
  2.123462d-92,-2.645285d-92,1.130217d-91,-3.611281d-91,2.603515d-91, &
  1.353475d-95,-1.431530d-95,7.301693d-95,-2.482972d-94,1.809258d-94, &

  7.611110d-99,-6.447446d-99,4.151517d-98,-1.505311d-97,1.108252d-97, &
  3.776086d-102,-2.312124d-102,2.076991d-101,-8.047103d-101,5.984041d-101, &
  1.652858d-105,-5.807915d-106,9.141453d-105,-3.793377d-104,2.848314d-104, &
  6.383114d-109,-3.979328d-110,3.538669d-108,-1.576881d-107,1.195189d-107  &
  /),(/5,64/)) 
  real(8) ::ept(0:4,0:63) = reshape((/ &
  8.596838d-01,-1.388595d-01,2.169498d-03,6.965218d-04,-1.649504d-05, &
  5.958830d-01,-1.225437d-01,5.744235d-03,4.570072d-04,-4.043767d-05, &
  3.767591d-01,-9.543768d-02,7.455990d-03,1.097695d-04,-4.288274d-05, &
  2.159250d-01,-6.559378d-02,7.174234d-03,-1.789366d-04,-2.739123d-05, &

  1.116118d-01,-3.978483d-02,5.594680d-03,-3.149562d-04,-6.903178d-06, &
  5.182992d-02,-2.129526d-02,3.660100d-03,-3.073040d-04,7.341110d-06, &
  2.155626d-02,-1.005905d-02,2.043258d-03,-2.242130d-04,1.206972d-05, &
  8.009942d-03,-4.193120d-03,9.827927d-04,-1.321552d-04,1.029767d-05, &

  2.654029d-03,-1.542484d-03,4.097516d-04,-6.500172d-05,6.444329d-06, &
  7.829383d-04,-5.007300d-04,1.486742d-04,-2.713633d-05,3.221155d-06, &
  2.053758d-04,-1.434441d-04,4.707814d-05,-9.712918d-06,1.332560d-06, &
  4.785484d-05,-3.626210d-05,1.303624d-05,-3.000312d-06,4.648605d-07, &

  9.896735d-06,-8.089284d-06,3.161504d-06,-8.034571d-07,1.382916d-07, &
  1.815281d-06,-1.592387d-06,6.722793d-07,-1.871347d-07,3.534286d-08, &
  2.951402d-07,-2.766066d-07,1.254635d-07,-3.800136d-08,7.799510d-09, &
  4.251394d-08,-4.239800d-08,2.056429d-08,-6.740792d-09,1.491821d-09, &

  5.423400d-09,-5.734431d-09,2.962071d-09,-1.045995d-09,2.480137d-10, &
  6.124832d-10,-6.843678d-10,3.751228d-10,-1.421573d-10,3.591675d-11, &
  6.121610d-11,-7.206672d-11,4.178523d-11,-1.693742d-11,4.538757d-12, &
  5.413406d-12,-6.696023d-12,4.095342d-12,-1.770558d-12,5.011958d-13, &

  4.234563d-13,-5.489417d-13,3.532645d-13,-1.624977d-13,4.841853d-14, &
  2.929488d-14,-3.970562d-14,2.682611d-14,-1.310094d-14,4.096081d-15, &
  1.792020d-15,-2.533862d-15,1.793718d-15,-9.282919d-16,3.036909d-16, &
  9.691555d-17,-1.426614d-16,1.056251d-16,-5.783264d-17,1.974705d-17, &

  4.633222d-18,-7.086129d-18,5.478535d-18,-3.169007d-18,1.126778d-18, &
  1.957752d-19,-3.105106d-19,2.503254d-19,-1.527814d-19,5.644996d-20, &
  7.310870d-21,-1.200308d-20,1.007719d-20,-6.482365d-21,2.484127d-21, &
  2.412535d-22,-4.092988d-22,3.574481d-22,-2.421126d-22,9.605947d-23, &

  7.034488d-24,-1.231123d-23,1.117286d-23,-7.961867d-24,3.265231d-24, &
  1.812217d-25,-3.266295d-25,3.077711d-25,-2.305732d-25,9.759539d-26, &
  4.124534d-27,-7.643289d-27,7.471937d-27,-5.881293d-27,2.565701d-27, &
  8.292724d-29,-1.577436d-28,1.598845d-28,-1.321519d-28,5.934055d-29, &

  1.472824d-30,-2.871071d-30,3.015572d-30,-2.616195d-30,1.207709d-30, &
  2.310529d-32,-4.608159d-32,5.013520d-32,-4.563700d-32,2.163341d-32, &
  3.201522d-34,-6.521849d-34,7.347542d-34,-7.015570d-34,3.411290d-34, &
  3.918014d-36,-8.138393d-36,9.492521d-36,-9.505002d-36,4.736036d-36, &

  4.234683d-38,-8.953489d-38,1.081117d-37,-1.135077d-37,5.790011d-38, &
  4.042071d-40,-8.683373d-40,1.085488d-39,-1.194865d-39,6.234065d-40, &
  3.407211d-42,-7.422992d-42,9.608267d-42,-1.108835d-41,5.912152d-42, &
  2.536245d-44,-5.592550d-44,7.497878d-44,-9.071945d-44,4.939153d-44, &

  1.667120d-46,-3.712960d-46,5.158315d-46,-6.544065d-46,3.635282d-46, &
  9.676384d-49,-2.171908d-48,3.128633d-48,-4.162313d-48,2.357464d-48, &
  4.959273d-51,-1.119168d-50,1.672930d-50,-2.334457d-50,1.347141d-50, &
  2.244246d-53,-5.079159d-53,7.886322d-53,-1.154578d-52,6.783885d-53, &

  8.967264d-56,-2.029683d-55,3.277481d-55,-5.035793d-55,3.010753d-55, &
  3.163566d-58,-7.139775d-58,1.200792d-57,-1.937039d-57,1.177700d-57, &
  9.853991d-61,-2.210144d-60,3.878363d-60,-6.571317d-60,4.060564d-60, &
  2.709925d-63,-6.018254d-63,1.104259d-62,-1.966195d-62,1.234121d-62, &

  6.579668d-66,-1.440911d-65,2.771545d-65,-5.188912d-65,3.306541d-65, &
  1.410408d-68,-3.031671d-68,6.131757d-68,-1.207857d-67,7.810147d-68, &
  2.669150d-71,-5.601622d-71,1.195752d-70,-2.480034d-70,1.626432d-70, &
  4.459451d-74,-9.081788d-74,2.055270d-73,-4.491742d-73,2.986240d-73, &

  6.577546d-77,-1.290621d-76,3.113455d-76,-7.176271d-76,4.834432d-76, &
  8.564746d-80,-1.605498d-79,4.156553d-79,-1.011394d-78,6.901076d-79, &
  9.845268d-83,-1.745143d-82,4.889988d-82,-1.257447d-81,8.686722d-82, &
  9.990742d-86,-1.653541d-85,5.069065d-85,-1.379163d-84,9.642269d-85, &

  8.949937d-89,-1.361100d-88,4.629670d-88,-1.334467d-87,9.438484d-88, &
  7.077656d-92,-9.684938d-92,3.724985d-91,-1.139132d-90,8.147784d-91, &
  4.940850d-95,-5.911214d-95,2.639941d-94,-8.578719d-94,6.203041d-94, &
  3.044754d-98,-3.054773d-98,1.647758d-97,-5.699805d-97,4.164960d-97, &

  1.656292d-101,-1.304183d-101,9.056203d-101,-3.341133d-100,2.466438d-100, &
  7.953387d-105,-4.350405d-105,4.381922d-104,-1.727946d-103,1.288234d-103, &
  3.371272d-108,-9.443834d-109,1.866149d-107,-7.884522d-107,5.934647d-107, &
  1.261417d-111,1.916964d-113,6.993126d-111,-3.174205d-110,2.411468d-110 &
  /),(/5,64/)) 

end module ewald_table


subroutine set_shift_coefficient(cutoff)
  use shift_function
  real(8) cutoff
  real(8) cutoff2

  cutoff2 = cutoff*cutoff
  r3co1  =  5.0d0/(cutoff2**2)
  r3co2  =  4.0d0/(cutoff2**2*cutoff)
  r8co1  = 10.0d0/(cutoff2**4*cutoff)
  r8co2  =  9.0d0/(cutoff2**5)
  r14co1 = 16.0d0/(cutoff2**7*cutoff)
  r14co2 = 15.0d0/(cutoff2**8)
  
  r1co1  =  r3co1/3.0d0*1.0d0
  r1co2  =  r3co2/4.0d0*1.0d0
  r6co1  =  r8co1/3.0d0*6.0d0
  r6co2  =  r8co2/4.0d0*6.0d0
  r12co1 = r14co1/3.0d0*12.0d0
  r12co2 = r14co2/4.0d0*12.0d0

  r1co0  = 1.0d0/cutoff      + r1co1*(cutoff2*cutoff)- r1co2*(cutoff2**2)
  r6co0  = 1.0d0/(cutoff2**3)+ r6co1*(cutoff2*cutoff)- r6co2*(cutoff2**2)
  r12co0 = 1.0d0/(cutoff2**6)+r12co1*(cutoff2*cutoff)-r12co2*(cutoff2**2)
end subroutine set_shift_coefficient

subroutine pairlistloopljfe1_i_posc_f(jindex,	&
                        lj,		&
                        npair,		&
                        posi,		&
                        chargei,	&
                        ljmp,           &
                        particlej,	&
                        cutoff2,	&
                        force,		&
                        energy)
#ifdef PRE_CALC_SF
  use shift_function
#endif
  integer jindex(0:MAX_PAIR-1)
  integer lj(0:MAX_PAIR-1)
  integer npair
  real(8) posi(0:2)
  real(8) chargei
  real(8) ljmp(0:3, 0:*)
  real(8) particlej(0:3, 0:*)
  real(8) cutoff2
  real(8) force(0:2)
  real(8) energy

  integer atj

  integer j
  real(8) fx, fy, fz
  real(8) e
  real(8) dx, dy, dz
  real(8) r2, inside
  real(8) pr, pr2, pr3, pr6, pr8, pr12, pr14
  real(8) dp
  real(8) r,r3,r4
#ifndef PRE_CALC_SF
  real(8) cutoff
  real(8) r6co0, r6co1, r6co2
  real(8) r12co0, r12co1, r12co2
  real(8) r8co1, r8co2
  real(8) r14co1, r14co2
  integer jid
  cutoff = sqrt(cutoff2)

  r8co1  = 10.0d0/(cutoff2**4*cutoff)
  r8co2  =  9.0d0/(cutoff2**5)
  r14co1 = 16.0d0/(cutoff2**7*cutoff)
  r14co2 = 15.0d0/(cutoff2**8)

  r6co1  =  r8co1/3.0d0*6.0d0
  r6co2  =  r8co2/4.0d0*6.0d0
  r12co1 = r14co1/3.0d0*12.0d0
  r12co2 = r14co2/4.0d0*12.0d0

  r6co0  = 1.0d0/(cutoff2**3)+ r6co1*(cutoff2*cutoff)- r6co2*(cutoff2**2)
  r12co0 = 1.0d0/(cutoff2**6)+r12co1*(cutoff2*cutoff)-r12co2*(cutoff2**2)
#endif

  fx=force(0)
  fy=force(1)
  fz=force(2)
  e=energy

  do j=0,npair-1
    jid = jindex(j)
    atj = lj(j)
    dx = posi(0)-particlej(0, jid)
    dy = posi(1)-particlej(1, jid)
    dz = posi(2)-particlej(2, jid)
    r2 = dx*dx+dy*dy+dz*dz
    inside = 1.0d0
#ifdef CUTOFF_MASK
    if (r2>cutoff2) then
      inside = 0.0d0
    else
    endif
#else
    r2 = min(r2,cutoff2)
#endif
#ifdef SEPARATE_DIVSQRT
    r  = sqrt(r2)
#else
    pr = 1.0d0/sqrt(r2)
    r = r2*pr
#endif
    r3 = r2*r
#ifdef SEPARATE_DIVSQRT
    pr = 1.0d0/r
#endif
    pr2 = pr*pr
    pr3 = pr2*pr
    pr6 = pr3*pr3
    pr8 = pr6*pr2
    pr12 = pr6*pr6
    pr14 = pr8*pr6
    e = e + (ljmp(0,atj)*(pr12+(r3*(r12co1-r12co2*r)-r12co0)) &
           - ljmp(1,atj)*(pr6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))  &
           )*inside
    dp = (ljmp(2,atj)*(pr14-r*(r14co1-r14co2*r)) &
        - ljmp(3,atj)*(pr8 -r*(r8co1 -r8co2 *r)) &
        )*inside
    fx = fx + dx*dp
    fy = fy + dy*dp
    fz = fz + dz*dp
  end do
  force(0) = fx
  force(1) = fy
  force(2) = fz
  energy = e
end subroutine pairlistloopljfe1_i_posc_f

subroutine pairlistloopljcfe1_i_posc_f(jindex,	&
                        lj,		&
                        npair,		&
                        posi,		&
                        chargei,	&
                        ljmp,           &
                        particlej,	&
                        cutoff2,	&
                        force,		&
                        energy)
#ifdef PRE_CALC_SF
  use shift_function
#endif
  integer jindex(0:MAX_PAIR-1)
  integer lj(0:MAX_PAIR-1)
  integer npair
  real(8) posi(0:2)
  real(8) chargei
  real(8) ljmp(0:3, 0:*)
  real(8) particlej(0:3, 0:*)
  real(8) cutoff2
  real(8) force(0:2)
  real(8) energy

  integer atj

  integer j
  real(8) fx, fy, fz
  real(8) e
  real(8) dx, dy, dz
  real(8) r2, cij, inside
  real(8) pr, pr2, pr3, pr6, pr8, pr12, pr14
  real(8) dp
  real(8) r,r3,r4
#ifndef PRE_CALC_SF
  real(8) cutoff
  real(8) r1co0, r1co1, r1co2
  real(8) r6co0, r6co1, r6co2
  real(8) r12co0, r12co1, r12co2
  real(8) r3co1, r3co2
  real(8) r8co1, r8co2
  real(8) r14co1, r14co2
  integer jid
  cutoff = sqrt(cutoff2)

  r3co1  =  5.0d0/(cutoff2**2)
  r3co2  =  4.0d0/(cutoff2**2*cutoff)
  r8co1  = 10.0d0/(cutoff2**4*cutoff)
  r8co2  =  9.0d0/(cutoff2**5)
  r14co1 = 16.0d0/(cutoff2**7*cutoff)
  r14co2 = 15.0d0/(cutoff2**8)

  r1co1  =  r3co1/3.0d0*1.0d0
  r1co2  =  r3co2/4.0d0*1.0d0
  r6co1  =  r8co1/3.0d0*6.0d0
  r6co2  =  r8co2/4.0d0*6.0d0
  r12co1 = r14co1/3.0d0*12.0d0
  r12co2 = r14co2/4.0d0*12.0d0

  r1co0  = 1.0d0/cutoff      + r1co1*(cutoff2*cutoff)- r1co2*(cutoff2**2)
  r6co0  = 1.0d0/(cutoff2**3)+ r6co1*(cutoff2*cutoff)- r6co2*(cutoff2**2)
  r12co0 = 1.0d0/(cutoff2**6)+r12co1*(cutoff2*cutoff)-r12co2*(cutoff2**2)
#endif

  fx=force(0)
  fy=force(1)
  fz=force(2)
  e=energy

  do j=0,npair-1
    jid = jindex(j)
    atj = lj(j)
    dx = posi(0)-particlej(0, jid)
    dy = posi(1)-particlej(1, jid)
    dz = posi(2)-particlej(2, jid)
    r2 = dx*dx+dy*dy+dz*dz
    inside = 1.0d0
#ifdef CUTOFF_MASK
    if (r2>cutoff2) then
      inside = 0.0d0
    else
    endif
#else
    r2 = min(r2,cutoff2)
#endif
#ifdef SEPARATE_DIVSQRT
    r  = sqrt(r2)
#else
    pr = 1.0d0/sqrt(r2)
    r = r2*pr
#endif
    cij = chargei*particlej(3, jid)
    r3 = r2*r
#ifdef SEPARATE_DIVSQRT
    pr = 1.0d0/r
#endif
    pr2 = pr*pr
    pr3 = pr2*pr
    pr6 = pr3*pr3
    pr8 = pr6*pr2
    pr12 = pr6*pr6
    pr14 = pr8*pr6
    e = e + (ljmp(0,atj)*(pr12+(r3*(r12co1-r12co2*r)-r12co0)) &
           - ljmp(1,atj)*(pr6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))  &
           + cij    *(pr  +(r3*(r1co1 -r1co2 *r)-r1co0 )))*inside
    dp = (ljmp(2,atj)*(pr14-r*(r14co1-r14co2*r)) &
        - ljmp(3,atj)*(pr8 -r*(r8co1 -r8co2 *r)) &
        + cij    *(pr3 -r*(r3co1 -r3co2 *r)))*inside
    fx = fx + dx*dp
    fy = fy + dy*dp
    fz = fz + dz*dp
  end do
  force(0) = fx
  force(1) = fy
  force(2) = fz
  energy = e
end subroutine pairlistloopljcfe1_i_posc_f

subroutine pairlistloopljewaldfe1_i_posc_f(alpha, &
                        jindex,	&
                        lj,		&
                        npair,		&
                        posi,		&
                        chargei,	&
                        ljmp,           &
                        particlej,	&
                        cutoff2,	&
                        force,		&
                        energy)
use ewald_table

  real(8) alpha
  integer jindex(0:MAX_PAIR-1)
  integer lj(0:MAX_PAIR-1)
  integer npair
  real(8) posi(0:2)
  real(8) chargei
  real(8) ljmp(0:3, 0:*)
  real(8) particlej(0:3, 0:*)
  real(8) cutoff2
  real(8) force(0:2)
  real(8) energy

  integer atj

  integer j
  real(8) fx, fy, fz
  real(8) e
  real(8) dx, dy, dz
  real(8) r2, cij, inside
  real(8) pr, pr2, pr3, pr6, pr8, pr12, pr14
  real(8) dp
  real(8) r,r3,r4
  real(8) ar, ddx
  real(8) ee, edp
  integer ix

  fx=force(0)
  fy=force(1)
  fz=force(2)
  e=energy

  do j=0,npair-1
    jid = jindex(j)
    atj = lj(j)
    dx = posi(0)-particlej(0, jid)
    dy = posi(1)-particlej(1, jid)
    dz = posi(2)-particlej(2, jid)
    r2 = dx*dx+dy*dy+dz*dz
    pr = 1.0d0/sqrt(r2)
    r = r2*pr
    ar = r*alpha
    cij = chargei*particlej(3, jid)
    inside = 0.0d0
    if(r2<cutoff2) then
       inside = 1.0d0
    else
    endif
    ix = ar*4.0d0
    ddx = (ar*4.0d0-ix)*2.0d0-1.0d0
    edp = eft(0,ix)+ddx*(eft(1,ix)+ddx*(eft(2,ix)+ddx*(eft(3,ix)+ddx*eft(4,ix))))
    ee = ept(0,ix)+ddx*(ept(1,ix)+ddx*(ept(2,ix)+ddx*(ept(3,ix)+ddx*ept(4,ix))))

    r3 = r2*r
    pr2 = pr*pr
    pr3 = pr2*pr
    pr6 = pr3*pr3
    pr8 = pr6*pr2
    pr12 = pr6*pr6
    pr14 = pr8*pr6
    e = e + (ljmp(0,atj)*(pr12) &
           - ljmp(1,atj)*(pr6)  &
           + cij    *(pr*ee))*inside
    dp = (ljmp(2,atj)*(pr14) &
        - ljmp(3,atj)*(pr8) &
        + cij    *(pr3*edp))*inside
    fx = fx + dx*dp
    fy = fy + dy*dp
    fz = fz + dz*dp
  end do
  force(0) = fx
  force(1) = fy
  force(2) = fz
  energy = e
end subroutine pairlistloopljewaldfe1_i_posc_f

subroutine pairlistloopljcfe1_i_f(jindex,	&
                        lj,		&
                        npair,		&
                        posi,		&
                        chargei,	&
                        particlej,	&
                        cutoff2,	&
                        force,		&
                        energy)
  integer jindex(0:MAX_PAIR-1)
  real(8) lj(0:3,0:MAX_PAIR-1)
  integer npair
  real(8) posi(0:2)
  real(8) chargei
  real(8) particlej(0:12, 0:*)
  real(8) cutoff2
  real(8) force(0:2)
  real(8) energy

  integer j
  real(8) fx, fy, fz
  real(8) e
  real(8) dx, dy, dz
  real(8) r2, cij, inside
  real(8) pr, pr2, pr3, pr6, pr8, pr12, pr14
  real(8) dp
  real(8) r,r3,r4
  real(8) cutoff
  real(8) r1co0, r1co1, r1co2
  real(8) r6co0, r6co1, r6co2
  real(8) r12co0, r12co1, r12co2
  real(8) r3co1, r3co2
  real(8) r8co1, r8co2
  real(8) r14co1, r14co2
  integer jid
  cutoff = sqrt(cutoff2)

  r3co1  =  5.0d0/(cutoff2**2)
  r3co2  =  4.0d0/(cutoff2**2*cutoff)
  r8co1  = 10.0d0/(cutoff2**4*cutoff)
  r8co2  =  9.0d0/(cutoff2**5)
  r14co1 = 16.0d0/(cutoff2**7*cutoff)
  r14co2 = 15.0d0/(cutoff2**8)

  r1co1  =  r3co1/3.0d0*1.0d0
  r1co2  =  r3co2/4.0d0*1.0d0
  r6co1  =  r8co1/3.0d0*6.0d0
  r6co2  =  r8co2/4.0d0*6.0d0
  r12co1 = r14co1/3.0d0*12.0d0
  r12co2 = r14co2/4.0d0*12.0d0

  r1co0  = 1.0d0/cutoff      + r1co1*(cutoff2*cutoff)- r1co2*(cutoff2**2)
  r6co0  = 1.0d0/(cutoff2**3)+ r6co1*(cutoff2*cutoff)- r6co2*(cutoff2**2)
  r12co0 = 1.0d0/(cutoff2**6)+r12co1*(cutoff2*cutoff)-r12co2*(cutoff2**2)

  fx=force(0)
  fy=force(1)
  fz=force(2)
  e=energy
!OCL SERIAL
!OCL SIMD
  do j=0,npair-1
    jid = jindex(j)
    dx = posi(0)-particlej(0, jid)
    dy = posi(1)-particlej(1, jid)
    dz = posi(2)-particlej(2, jid)
    r2 = dx*dx+dy*dy+dz*dz
    inside = 1.0d0
#ifdef CUTOFF_MASK
    if (r2>cutoff2) then
      inside = 0.0d0
    else
    endif
#else
    r2 = min(r2,cutoff2)
#endif
    r  = sqrt(r2)
    r3 = r2*r
    cij = chargei*particlej(9, jid)
    pr = 1.0d0/r
    pr2 = pr*pr
    pr3 = pr2*pr
    pr6 = pr3*pr3
    pr8 = pr6*pr2
    pr12 = pr6*pr6
    pr14 = pr8*pr6
    e = e + (lj(0,j)*(pr12+(r3*(r12co1-r12co2*r)-r12co0)) &
           - lj(1,j)*(pr6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))  &
           + cij    *(pr  +(r3*(r1co1 -r1co2 *r)-r1co0 )))*inside
    dp = (lj(2,j)*(pr14-r*(r14co1-r14co2*r)) &
        - lj(3,j)*(pr8 -r*(r8co1 -r8co2 *r)) &
        + cij    *(pr3 -r*(r3co1 -r3co2 *r)))*inside
    fx = fx + dx*dp
    fy = fy + dy*dp
    fz = fz + dz*dp
  end do
  force(0) = fx
  force(1) = fy
  force(2) = fz
  energy = e
end subroutine pairlistloopljcfe1_i_f

subroutine pairlistloopljcfe1_f(pos,	&
                        charge,		&
                        lj,		&
                        npair,		&
                        posi,		&
                        chargei,	&
                        cutoff2,	&
                        force,		&
                        energy)
  real(8) pos(0:2,0:MAX_PAIR-1)
  real(8) charge(0:MAX_PAIR-1)
  real(8) lj(0:3,0:MAX_PAIR-1)
  integer npair
  real(8) posi(0:2)
  real(8) chargei
  real(8) cutoff2
  real(8) force(0:2)
  real(8) energy

  integer j
  real(8) fx, fy, fz
  real(8) e
  real(8) dx, dy, dz
  real(8) r2, cij, inside
  real(8) pr, pr2, pr3, pr6, pr8, pr12, pr14
  real(8) dp
  real(8) r,r3,r4
  real(8) cutoff
  real(8) r1co0, r1co1, r1co2
  real(8) r6co0, r6co1, r6co2
  real(8) r12co0, r12co1, r12co2
  real(8) r3co1, r3co2
  real(8) r8co1, r8co2
  real(8) r14co1, r14co2
  cutoff = sqrt(cutoff2)

  r3co1  =  5.0d0/(cutoff2**2)
  r3co2  =  4.0d0/(cutoff2**2*cutoff)
  r8co1  = 10.0d0/(cutoff2**4*cutoff)
  r8co2  =  9.0d0/(cutoff2**5)
  r14co1 = 16.0d0/(cutoff2**7*cutoff)
  r14co2 = 15.0d0/(cutoff2**8)

  r1co1  =  r3co1/3.0d0*1.0d0
  r1co2  =  r3co2/4.0d0*1.0d0
  r6co1  =  r8co1/3.0d0*6.0d0
  r6co2  =  r8co2/4.0d0*6.0d0
  r12co1 = r14co1/3.0d0*12.0d0
  r12co2 = r14co2/4.0d0*12.0d0

  r1co0  = 1.0d0/cutoff      + r1co1*(cutoff2*cutoff)- r1co2*(cutoff2**2)
  r6co0  = 1.0d0/(cutoff2**3)+ r6co1*(cutoff2*cutoff)- r6co2*(cutoff2**2)
  r12co0 = 1.0d0/(cutoff2**6)+r12co1*(cutoff2*cutoff)-r12co2*(cutoff2**2)

  fx=force(0)
  fy=force(1)
  fz=force(2)
  e=energy
!OCL SERIAL
  do j=0,npair-1
    dx = posi(0)-pos(0, j)
    dy = posi(1)-pos(1, j)
    dz = posi(2)-pos(2, j)
    r2 = dx*dx+dy*dy+dz*dz
    inside = 1.0d0
#ifdef CUTOFF_MASK
    if (r2>cutoff2) then
      inside = 0.0d0
    else
    endif
#else
    r2 = min(r2,cutoff2)
#endif
    r  = sqrt(r2)
    r3 = r2*r
    cij = chargei*charge(j)
    pr = 1.0d0/r
    pr2 = pr*pr
    pr3 = pr2*pr
    pr6 = pr3*pr3
    pr8 = pr6*pr2
    pr12 = pr6*pr6
    pr14 = pr8*pr6
    e = e + (lj(0,j)*(pr12+(r3*(r12co1-r12co2*r)-r12co0)) &
           - lj(1,j)*(pr6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))  &
           + cij    *(pr  +(r3*(r1co1 -r1co2 *r)-r1co0 )))*inside
    dp = (lj(2,j)*(pr14-r*(r14co1-r14co2*r)) &
        - lj(3,j)*(pr8 -r*(r8co1 -r8co2 *r)) &
        + cij    *(pr3 -r*(r3co1 -r3co2 *r)))*inside
    fx = fx + dx*dp
    fy = fy + dy*dp
    fz = fz + dz*dp
  end do
  force(0) = fx
  force(1) = fy
  force(2) = fz
  energy = e
end subroutine pairlistloopljcfe1_f

