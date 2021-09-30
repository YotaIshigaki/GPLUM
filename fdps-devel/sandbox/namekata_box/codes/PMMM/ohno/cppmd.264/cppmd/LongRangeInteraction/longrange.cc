/* g++ -O3 c.c -lfftw -lm
 * icpc c.c -lfftw -lm
 *  
 * doamin size [0...Lx][0...Ly][0...Lz]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw.h>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/* Ewald parameter */
//#define beta 0.1

#define Lx 62.0
#define Ly Lx
#define Lz Lx

/* water model */
#define bond 0.957
#define angle 104.5
#define theta M_PI/180.0*angle/2.0

/* particle array */
#define Ncp 6
double x[Ncp],y[Ncp],z[Ncp],charge[Ncp];
double fx[Ncp],fy[Ncp],fz[Ncp];

double directsum(int Reciprocal);
double Ewald(double alpha,int Imax,int Jmax,int Kmax);
double PME(double beta,int K1, int K2,int K3);
double realspace(double beta,int K1, int K2,int K3);

double Lagrange(int order, double u, int k);

/**********************************************************************/
double directsum(int Reciprocal){
  for(int i=0; i<Ncp; i++) fx[i] = fy[i] = fz[i] = 0.0;
  double potpl = 0.0, potmi = 0.0;
  for(int i=0; i<Ncp; i++){
    for(int nx = -Reciprocal; nx <= Reciprocal; nx++){
      for(int ny = -Reciprocal; ny <= Reciprocal; ny++){
        for(int nz = -Reciprocal; nz <= Reciprocal; nz++){
          for(int j=0; j<Ncp; j++){
            double xj = x[i] - x[j] + Lx*nx;
            double yj = y[i] - y[j] + Ly*ny;
            double zj = z[i] - z[j] + Lz*nz;
            double r2 = SQR(xj) + SQR(yj) + SQR(zj);
            double qq = charge[i] * charge[j];
            if(r2>0.0){
              double r = sqrt(r2);
              if(qq > 0.0) potpl += 0.5 * qq/r;
              if(qq < 0.0) potmi += 0.5 * qq/r;
              fx[i] += qq * xj/(r*r*r);
              fy[i] += qq * yj/(r*r*r);
              fz[i] += qq * zj/(r*r*r);
            }
          }
        }
      }
    }
  }
  printf("# f = %e %e %e\n",fx[1],fy[1],fz[1]);
  //printf("abs(f) = %e\n",sqrt(SQR(fx[1])+SQR(fy[1])+SQR(fz[1])));
  return potpl + potmi;
}

/**********************************************************************/
double Ewald(double alpha,int Imax,int Jmax,int Kmax){
  for(int i=0; i<Ncp; i++) fx[i] = fy[i] = fz[i] = 0.0;

  double pot = 0;
  /* direct part */
  for(int i=0; i<Ncp; i++){
    for(int j=0; j<Ncp; j++){
      double xj = x[i] - x[j];
      double yj = y[i] - y[j];
      double zj = z[i] - z[j];
      double r2 = SQR(xj) + SQR(yj) + SQR(zj);
      if(r2 > 0.0){
        double r = sqrt(r2);
        double ar = alpha * r;
        double qq = charge[i] * charge[j];
        pot += 0.5 * qq * erfc(ar)/r;
        double coeff = qq 
          * (erfc(ar) + 2.0*ar*exp(-ar*ar)/sqrt(M_PI))/(r*r*r);
        fx[i] += xj * coeff;
        fy[i] += yj * coeff;
        fz[i] += zj * coeff;
      }
    }
  }
  /* Fourier part */
  for(int k=-Kmax; k<=Kmax; k++){
    for(int j=-Jmax; j<=Jmax; j++){
      for(int i=-Imax; i<=Imax; i++){
        double mx = (double)i/Lx;
        double my = (double)j/Ly;
        double mz = (double)k/Lz;
        double m2 = SQR(mx) + SQR(my) + SQR(mz);
        if(m2 > 0.0){
          double s = 0.0, c = 0.0;
          for(int p=0; p<Ncp; p++){
            double phase = 2.0 * M_PI * (mx * x[p] + my * y[p] + mz * z[p]);
            s += charge[p] * sin(phase);
            c += charge[p] * cos(phase); 
          }
          pot += 1.0/(2.0*M_PI*Lx*Ly*Lz) * exp(-SQR(M_PI/alpha) * m2)/m2 
            *(SQR(s) + SQR(c));
          for(int p=0; p<Ncp; p++){
            double phase = 2.0 * M_PI * (mx * x[p] + my * y[p] + mz * z[p]);
            double coeff = 2.0 * charge[p]/(Lx * Ly * Lz)
                * exp(-SQR(M_PI/alpha)*m2)/m2 
                * (sin(phase) * c - cos(phase) * s);
            fx[p] += coeff * mx;
            fy[p] += coeff * my;
            fz[p] += coeff * mz;
          }
        }
      }
    }
  }
  /* self energy part */
  for(int i=0; i<Ncp; i++) pot -=alpha/sqrt(M_PI)*SQR(charge[i]);
  printf("# f = %e %e %e\n",fx[1],fy[1],fz[1]);
  return pot;
}

/**********************************************************************/
double Lagrange(int order1, double u, int k){
  int p = order1/2;
  if(u <= -p || u >= p){
    //printf("bingo! u=%e k=%d\n",u,k);
    return 0.0;
  }else{
    double W2p = 1.0;
    for(int j=-p; j<=p-1; j++){
      if(j != k) W2p *=(u + (double)(j-k))/(double)(j-k);
    }
    //if(W2p>1.0)printf("bingo!");
    return (W2p);
  }
}

/**********************************************************************/
/* Particle mesh Ewald : Essmann et al (1995) 
*/
double PME(double beta, int K1, int K2, int K3){
#define PMEorder 2 // this must be even number 
#define PMEhalf (PMEorder/2)
  int i,j,k,p,k1,k2,k3,n1,n2,n3,m1,m2,m3,k1n1K1,k2n2K2,k3n3K3;
  double r2,r,pot,u1,u2,u3,qq,mm,Wx,Wy,Wz,u1k1,u2k2,u3k3;
  double C[K1][K2][K3];
  fftw_complex Q[K1][K2][K3];
  fftwnd_plan plan;

  plan = fftw3d_create_plan(K1,K2,K3,
      FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);

  pot = 0.0;
  /* direct part */
  for(i=0; i<Ncp; i++){
    for(j=0; j<Ncp; j++){
      r2 = SQR(x[i]-x[j]) + SQR(y[i]-y[j]) + SQR(z[i]-z[j]);
      if(r2 > 0.0){
        r = sqrt(r2);
        qq = charge[i]*charge[j];
        pot += 0.5 * qq * erfc(beta*r)/r;
      }
    }
  }

  /* charge assign */

  for(i=0; i<K1; i++){
    for(j=0; j<K2; j++){
      for(k=0; k<K3; k++){
        Q[i][j][k].re = 0.0; 
        Q[i][j][k].im = 0.0; 
      }
    }
  }


  for(p=0; p<Ncp; p++){
    u1 = x[p]/Lx * K1; 
    u2 = y[p]/Ly * K2; 
    u3 = z[p]/Lz * K3; 

    for(k1n1K1 = -PMEhalf + u1; k1n1K1 < PMEhalf + u1; k1n1K1++){
      for(k2n2K2 = -PMEhalf + u2; k2n2K2 < PMEhalf + u2; k2n2K2++){
        for(k3n3K3 = -PMEhalf + u3; k3n3K3 < PMEhalf + u3; k3n3K3++){

          for(n1 = -1; n1 <= 1; n1++){
            for(n2 = -1; n2 <= 1; n2++){
              for(n3 = -1; n3 <= 1; n3++){
                k1 = k1n1K1 - n1*K1;
                k2 = k2n2K2 - n2*K2;
                k3 = k3n3K3 - n3*K3;
                if(k1>=0 && k1<K1 && k2>=0 && k2<K2 && k3>=0 && k3<K3){
                  Wx = Lagrange(PMEorder, u1 - k1n1K1, k1); 
                  Wy = Lagrange(PMEorder, u2 - k2n2K2, k2);
                  Wz = Lagrange(PMEorder, u3 - k3n3K3, k3);
                  Q[k1][k2][k3].re += charge[p] * Wx * Wy * Wz; 
                  printf("p=%d u=%e %e %e k=%d %d %d W=%e\n",p,u1,u2,u3,k1,k2,k3,Wx*Wy*Wz);
                }

              }
            }
          }
        }
      }
    }

  }

  /*  
      for(k1=0; k1<K1; k1++){
      k2 = K2/2;
      for(k3=0; k3<K3; k3++){
      printf("%d %d %e\n",k1,k3,Q[k1][k2][k3].re);
      }
      printf("\n");
      }
      */

  fftwnd_one(plan, &Q[0][0][0], NULL); /*  compute F(Q)  */

  /* kernel */

  for(m1=0; m1<K1; m1++){
    for(m2=0; m2<K2; m2++){
      for(m3=0; m3<K3; m3++){

        i= (m1 > K1/2) ? (m1 - K1) : m1;
        j= (m2 > K2/2) ? (m2 - K2) : m2;
        k= (m3 > K3/2) ? (m3 - K3) : m3;

        i = m1; j = m2; k = m3;
        mm = SQR((double)i/Lx) + SQR((double)j/Ly) + SQR((double)k/Lz);
        if(mm > 0.0){
          C[m1][m2][m3] =1.0/(M_PI*Lx*Ly*Lz)
            * exp(-SQR(M_PI/beta) * mm)/mm; 
        }else{
          C[m1][m2][m3] =0.0;
        }
      }
    }
  }

  /* convolution */

  for(m1=0; m1<K1; m1++){
    for(m2=0; m2<K2; m2++){
      for(m3=0; m3<K3; m3++){
        pot += 0.5 * C[m1][m2][m3] 
          * (SQR(Q[m1][m2][m3].re) + SQR(Q[m1][m2][m3].im));
        //Q[m1][m2][m3].re *=C[m1][m2][m3];
        //Q[m1][m2][m3].im *=C[m1][m2][m3];
      }
    }
  }

  /* self energy part */
  for(i=0;i<Ncp;i++){
    pot -=beta/sqrt(M_PI)*SQR(charge[i]);
  }

  fftwnd_destroy_plan(plan);
  return (pot);
}

/**********************************************************************/
double realspace(double beta, int K1, int K2, int K3){
  int i,j,k,p,x1,x2,x3,nx,ny,nz;
  double r2,r,pot,kk,mm,dx1,dx2,dx3,coeff,xj,yj,zj,qq,ar;
  double dkx,dky,dkz,dx1sq,dx2sq,dx3sq,pcoeff,chargeNeutrality;
  double C[K1][K2][K3];
  fftw_complex Q[K1][K2][K3];
  fftwnd_plan plan;

  plan = fftw3d_create_plan(K1,K2,K3,
      FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);

  pot = 0.0;
  for(i=0; i<Ncp; i++){
    fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
  }
#if 0
  /* direct part */
  for(i=0; i<Ncp; i++){
    for(j=0; j<Ncp; j++){
      xj = x[i] - x[j];
      yj = y[i] - y[j];
      zj = z[i] - z[j];
      r2 = SQR(xj) + SQR(yj) + SQR(zj);
      if(r2 > 0.0){
        r = sqrt(r2);
        //ar = beta/sqrt(2.0)*r;
        ar = beta * r;
        qq = charge[i] * charge[j];
        pot += 0.5 * qq * erfc(ar)/r;
        coeff = qq * 
          (erfc(ar) + 2.0*ar*exp(-ar*ar)/sqrt(M_PI))/(r*r*r);
        fx[i] += xj * coeff;
        fy[i] += yj * coeff;
        fz[i] += zj * coeff;
      }
    }
  }
#endif
  //  printf("%e\n",pot);
  //coeff = CUBE(beta/sqrt(2.0*M_PI));
  coeff = CUBE(beta/sqrt(M_PI));

#if 0
int cnt = 0;
#endif
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        Q[x1][x2][x3].re = 0.0;
        Q[x1][x2][x3].im = 0.0;
        xj = (double)x1/K1 * Lx;
        yj = (double)x2/K2 * Ly;
        zj = (double)x3/K3 * Lz;
        for(p=0; p<Ncp; p++){
#if 0
          const double mm = SQR(xj-x[p]) + SQR(yj-y[p]) + SQR(zj-z[p]);
          const double chg = coeff*charge[p]*exp(-mm*SQR(beta));
          Q[x1][x2][x3].re += chg;
          printf("--- %d, %d, %d ---\n",x1,x2,x3);
          printf("  [%f,%f,%f] - (%f,%f,%f)\n", xj,yj,zj,x[p],y[p],z[p]);
          printf("P: (%f,%f,%f),%f\tr2:%f,chg:%g,Q:%g,%d\n", x[p],y[p],z[p],charge[p],mm,chg,Q[x1][x2][x3].re,++cnt);
#else
          mm = SQR(xj-x[p]) + SQR(yj-y[p]) + SQR(zj-z[p]);
          //Q[x1][x2][x3].re += coeff*charge[p]*exp(-mm*SQR(beta)/2.0);
          Q[x1][x2][x3].re += coeff*charge[p]*exp(-mm*SQR(beta));
#endif
        }
      }
    }
  }

  chargeNeutrality = 0.0;
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        chargeNeutrality +=Q[x1][x2][x3].re;
      }
    }
  }
  chargeNeutrality /= (double)(K1*K2*K3);
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        Q[x1][x2][x3].re -=chargeNeutrality;
        C[x1][x2][x3] = Q[x1][x2][x3].re;
      }
    }
  }

#if 0
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        printf("  % f", C[x1][x2][x3]);
      }
      putchar('\n');
    }
    puts("\n");
  }
#endif

  fftwnd_one(plan, &Q[0][0][0], NULL); /*  compute  */
#if 0
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
#if 0
        printf("  % f", Q[x1][x2][x3].re);
#else
        printf("  % g", Q[x1][x2][x3].im);
#endif
      }
      putchar('\n');
    }
    puts("\n");
  }
#endif


  dkx = 2.0 * M_PI/(double)K1;
  dky = 2.0 * M_PI/(double)K2;
  dkz = 2.0 * M_PI/(double)K3;
  dx1sq = SQR(Lx/(double)K1);
  dx2sq = SQR(Ly/(double)K2);
  dx3sq = SQR(Lz/(double)K3);

  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        i = (x1 > K1/2) ? (x1 - K1) : x1;
        j = (x2 > K2/2) ? (x2 - K2) : x2;
        k = (x3 > K3/2) ? (x3 - K3) : x3;
        pcoeff =-(2.0 * cos((double)i * dkx) - 2.0)/dx1sq 
          -(2.0 * cos((double)j * dky) - 2.0)/dx2sq
          -(2.0 * cos((double)k * dkz) - 2.0)/dx3sq;

#if 0
        printf("k: %d, pcoeff.z: %g\n", k, -(2.0 * cos((double)k * dkz) - 2.0)/dx3sq);
#endif
        //pcoeff = SQR(2.0*M_PI*i/Lx)+SQR(2.0*M_PI*j/Ly)+SQR(2.0*M_PI*k/Lz);
        if(pcoeff != 0.0){
          Q[x1][x2][x3].re *= 4.0*M_PI/pcoeff; 
          Q[x1][x2][x3].im *= 4.0*M_PI/pcoeff; 
        }else{
          Q[x1][x2][x3].re = 0.0; 
          Q[x1][x2][x3].im = 0.0; 
        }
      }
    }
  }

  fftwnd_one(plan, &Q[0][0][0], NULL); /*  compute   */
#if 0
  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
#if 1
        printf("  % f", Q[x1][x2][x3].re);
#else
        printf("  % g", Q[x1][x2][x3].im);
#endif
      }
      putchar('\n');
    }
    puts("\n");
  }
#endif

  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
        Q[x1][x2][x3].re /=(double)(K1 * K2 * K3);
      }
    }
  }
  dx1 = Lx/(double)K1;
  dx2 = Ly/(double)K2;
  dx3 = Lz/(double)K3;

  for(x1=0; x1<K1; x1++){
    for(x2=0; x2<K2; x2++){
      for(x3=0; x3<K3; x3++){
#if 0
        printf("(%d,%d,%d)\t%f\t%f\n",x1,x2,x3,C[x1][x2][x3],Q[x1][x2][x3].re);
#endif
        pot += 0.5 * C[x1][x2][x3] * Q[x1][x2][x3].re * dx1 * dx2 * dx3;
        coeff = 2.0 * CUBE(beta/sqrt(M_PI)) * SQR(beta) 
          * Q[x1][x2][x3].re * dx1 * dx2 * dx3;
        xj = (double)x1 * dx1;
        yj = (double)x2 * dx2;
        zj = (double)x3 * dx3;
        for(p=0; p<Ncp; p++){
          r2 = SQR(xj - x[p]) + SQR(yj - y[p]) + SQR(zj - z[p]);
          fx[p] -= coeff * charge[p] * exp(-SQR(beta)*r2) * (xj - x[p]);
          fy[p] -= coeff * charge[p] * exp(-SQR(beta)*r2) * (yj - y[p]);
          fz[p] -= coeff * charge[p] * exp(-SQR(beta)*r2) * (zj - z[p]);
        }
      }
    }
  }
#if 1
  for(int p=0; p<Ncp; p++){
    printf("%d: (% g,% g,% g)\n",p,fx[p],fy[p],fz[p]);
  }
  printf("potential: %g\n", pot);
#endif

#if 0
  /* self energy part */
  for(i=0; i<Ncp; i++){
    //    pot -=beta/sqrt(2.0*M_PI)*SQR(charge[i]);
    pot -=beta/sqrt(M_PI)*SQR(charge[i]);
  }
#endif

  printf("# f = %e %e %e\n",fx[1],fy[1],fz[1]);



  /*  
      for(x1=0; x1<K1; x1++){
      x2 = K2/2;
      for(x3=0; x3<K3; x3++){
      printf("%d %d %e %e %e\n",x1,x3,
      C[x1][x2][x3],Q[x1][x2][x3].re,Q[x1][x2][x3].im);
      }
      printf("\n");
      }
      */
  fftwnd_destroy_plan(plan);

  return (pot);
}

/**********************************************************************/

int main(){
  /* water 1 */
  /* Oxygen         Hydrogen       Hydrogen */
  charge[0] = -0.82; charge[1] = 0.41; charge[2] = 0.41; 
  x[0]=Lx/2;      x[1]=x[0]+bond*sin(theta);  x[2]=x[0]-bond*sin(theta);
  y[0]=Ly/2;      y[1]=y[0];                  y[2]=y[0];
  z[0]=Lz/2+1.0;  z[1]=z[0]+bond*cos(theta);  z[2]=z[0]+bond*cos(theta);

  /* water 2 */
  /* Oxygen         Hydrogen       Hydrogen */
  charge[3] = -0.82; charge[4] = 0.41; charge[5] = 0.41; 
  x[3]= Lx/2;     x[4]=x[3]+bond*sin(theta);  x[5]=x[3]-bond*sin(theta);
  y[3]= Ly/2;     y[4]=y[3];                  y[5]=y[3];
  z[3]= Lz/2-1.0; z[4]=z[3]-bond*cos(theta);  z[5]=z[3]-bond*cos(theta);

  //printf("#directsum pot=%.15e\n", directsum(64));
  //printf("#Ewaldsum  pot=%.15e\n",  Ewald(0.6, 64, 64, 64));
  //printf("#PME       pot=%.15e\n",  PME(0.25,32,32,32));
  printf("#realspace pot=%.15e\n",  realspace(0.25,32,32,32));
  return 0;
}
