#include <iostream>
#include <iomanip>

typedef union {
    float f;
    int i;
    unsigned int ui;
} UN32;

typedef union {
    double d;
    long long int l;
    unsigned long long int ul;
} UN64;

template <class T>
void calc_bitseq(T data, char *bitseq) {
    const int size = sizeof(T)*8;
    for (int i=0; i<size; i++) bitseq[i]='0';
    bitseq[size] = '\0';
    for (int i=0; i<size; i++) {
        if ((i != 0) && (data == 0)) break;
        if (data % 2 == 0) {
            bitseq[(size-1)-i] = '0';
        } else {
            bitseq[(size-1)-i] = '1';
        }
        data /= 2;
    }
}

int main (int argc, char *argv[]) {
    //* Local parameters
    const double eps = 1.23456789012345e-10;
    //* Local variables
    UN64 x,y,xrel,yrel;
    UN32 xf,yf,xrelf,yrelf;
    UN64 cen;
    unsigned long long int mask;
    char bufc[256],bufx[256],bufy[256];

    //* Some check
    printf("sizeof(int)           = %d\n",sizeof(int));
    printf("sizeof(long int)      = %d\n",sizeof(long int));
    printf("sizeof(long long int) = %d\n",sizeof(long long int));

    x.d = 1.0e0;
    y.d = x.d + eps;
    xf.f = x.d;
    yf.f = y.d;
    cen.d = 0.5*(x.d+y.d);
    calc_bitseq(cen.ul,bufc);
    printf("cen (org)     = %lf\n",cen.d);
    printf("cen (org,bit) = %s\n",bufc);
    mask = 63;
    cen.ul &= ~mask; // mask the lowest 6 bits
    calc_bitseq(cen.ul,bufc);
    printf("cen (new)     = %lf\n",cen.d);
    printf("cen (new,bit) = %s\n",bufc);
    xrel.d = x.d - cen.d;
    yrel.d = y.d - cen.d;
    xrelf.f = xrel.d;
    yrelf.f = yrel.d;
 
    //* Output
    printf("-------\n");
    printf("x         = %-25.16lf\n",x.d);
    printf("y         = %-25.16lf\n",y.d);
    printf("x  (hex)  = %lx\n",x.ul); 
    printf("y  (hex)  = %lx\n",y.ul); 
    calc_bitseq(x.ul,bufx);
    calc_bitseq(y.ul,bufy);
    printf("x  (bin)  = %s\n",bufx);
    printf("y  (bin)  = %s\n",bufy);
    printf("diff.     = %-25.16lf\n",y.d-x.d);
    printf("-------\n");
    printf("xf        = %-25.16f\n",xf.f);
    printf("yf        = %-25.16f\n",yf.f);
    printf("xf (hex)  = %x\n",xf.ui); 
    printf("yf (hex)  = %x\n",yf.ui); 
    calc_bitseq(xf.ui,bufx);
    calc_bitseq(yf.ui,bufy);
    printf("xf (bin)  = %s\n",bufx); 
    printf("yf (bin)  = %s\n",bufy); 
    printf("diff.     = %-25.16lf\n",yf.f-xf.f);
    printf("-------\n");
    printf("xr        = %-25.16lf\n",xrel.d);
    printf("yr        = %-25.16lf\n",yrel.d);
    printf("xr (hex)  = %lx\n",xrel.ul); 
    printf("yr (hex)  = %lx\n",yrel.ul); 
    calc_bitseq(xrel.ul,bufx);
    calc_bitseq(yrel.ul,bufy);
    printf("xr (bin)  = %s\n",bufx);
    printf("yr (bin)  = %s\n",bufy);
    printf("diff.     = %-25.16lf\n",yrel.d-xrel.d);
    printf("-------\n");
    printf("xrf       = %-25.16f\n",xrelf.f);
    printf("yrf       = %-25.16f\n",yrelf.f);
    printf("xrf (hex) = %x\n",xrelf.ui); 
    printf("yrf (hex) = %x\n",yrelf.ui); 
    calc_bitseq(xrelf.ui,bufx);
    calc_bitseq(yrelf.ui,bufy);
    printf("xrf (bin) = %s\n",bufx); 
    printf("yrf (bin) = %s\n",bufy); 
    printf("diff.     = %-25.16lf\n",yrelf.f-xrelf.f);
    printf("-------\n");

    return 0;
}


