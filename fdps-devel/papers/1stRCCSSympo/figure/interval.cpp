#include<cstdio>
#include<cmath>

int main(){
  for(int i=0;i<16;i++){
    float x = powf(2,i-7);
    printf("%lf %lf\n",x,sqrt(erfc(sqrt(x))/x + 2*exp(-x)/sqrt(x)/M_PI));
  }
}
