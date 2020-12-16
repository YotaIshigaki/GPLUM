#include <iostream>
#include <cmath>
#include "discpp.h"

#define SQ(x) ((x)*(x))

constexpr int N_lo = 5;
constexpr int N_hi = 256;
constexpr int NX = 256;
constexpr int NY = 256;

template <class T>
void locate(const T xx[], const int n, const T x, int & j) {
   int jl = -1;
   int ju = n;
   while( (ju -jl) > 1) {
      const int jm = (ju+jl)/2;
      if (x > xx[jm]) jl = jm;
      else ju = jm;
   }
   j=jl;
}

template <int n>
class ColorTable {
public:
    double s[n]; // position between [0,1]
    double r[n]; // R of RGB
    double g[n]; // G of RGB
    double b[n]; // B of RGB

    int size() const { return n; }
    void set(double s[], double r[], double g[], double b[]) {
        for (int i = 0; i < n; i++) {
            this->s[i] = s[i];
            this->r[i] = r[i];
            this->g[i] = g[i];
            this->b[i] = b[i];
        }
    }
    template <int m>
    void interpTo(ColorTable<m> & dest) const {
        const double ds = 1.0/static_cast<double>(dest.size()-1);
        for (int i = 0; i < dest.size(); i++) {
            const double s = ds * i;
            dest.s[i] = s;
            int j;
            locate(this->s, this->size(), s, j);
            if ((0 <= j) && (j <= this->size() - 2)) {
                const double offset = s - this->s[j];
                const double intvl = this->s[j+1] - this->s[j];
                dest.r[i] = this->r[j] + (this->r[j+1] - this->r[j])/intvl * offset;
                dest.g[i] = this->g[j] + (this->g[j+1] - this->g[j])/intvl * offset;
                dest.b[i] = this->b[j] + (this->b[j+1] - this->b[j])/intvl * offset;
            } else if (j < 0) {
                dest.r[i] = this->r[0];
                dest.g[i] = this->g[0];
                dest.b[i] = this->b[0];
            } else if (this->size() - 1 <= j) {
                dest.r[i] = this->r[this->size()-1];
                dest.g[i] = this->g[this->size()-1];
                dest.b[i] = this->b[this->size()-1];
            }
        }
    }
};

int main(int argc, char *argv[]) {

    // Make a color table
    ColorTable<N_lo> ctbl_lo;
    {
        double s[] = {0.0, 0.3, 0.55, 0.8, 1.0};
        double r[] = {0.0, 0.5,  1.0, 1.0, 1.0};
        double g[] = {0.0, 0.0,  0.5, 1.0, 1.0};
        double b[] = {0.0, 0.0,  0.0, 0.3, 1.0};
        ctbl_lo.set(s, r, g, b);
    }
    ColorTable<N_hi> ctbl_hi;
    ctbl_lo.interpTo(ctbl_hi);

    // Make a figure using the color table above
    Dislin g;
    g.metafl("PS");
    int page_size[] = {2500, 2500};
    g.page(page_size[0], page_size[1]);
    g.pagmod("LAND");
    char filename[] = "test.ps";
    g.setfil(filename);
    g.disini();
    g.texmod("ON");

    g.labdig(2,"X");
    g.labdig(2,"Y");
    g.labdig(1,"Z");
    g.titlin("Test title", 3);
    g.name("Test label (x)","X");
    g.name("Test label (y)","Y");
    g.name("Test label (cb)","Z");
    g.autres(NX, NY);
    int pos[] = {750, 1750};
    g.axspos(pos[0], pos[1]);
    int len[] = {1000, 1000, 1000};
    g.ax3len(len[0], len[1], len[2]);
    g.setrgb(0.0, 0.0, 0.0);
    g.complx();
    //g.myvlt(ctbl_lo.r, ctbl_lo.g, ctbl_lo.b, N_lo); // use coarse color table
    g.myvlt(ctbl_hi.r, ctbl_hi.g, ctbl_hi.b, N_hi); // use fine color table
    double xmin[] = {0.0, 0.0, 0.0};
    double xmax[] = {1.0, 1.0, 1.0}; 
    double xtics[] = {0.25, 0.25, 0.1};
    g.graf3(xmin[0], xmax[0], xmin[0], xtics[0],
            xmin[1], xmax[1], xmin[1], xtics[1],
            xmin[2], xmax[2], xmin[2], xtics[2]);
    double data[NX*NY]; 
    int idx=0;
    for (int ix = 0; ix < NX; ix++) {
        for (int iy = 0; iy < NY; iy++) {
            const double pi = 4.0*std::atan(1.0);
            const double x = (double)ix / (double) NX;
            const double y = (double)iy / (double) NY;
            const double phase_x = 2.0 * pi * x;
            const double phase_y = 4.0 * pi * y;
            data[idx] = SQ(std::sin(phase_x))
                      * SQ(std::sin(phase_y));
            idx++;
        }
    }
    g.crvmat(data, NX, NY, 1, 1);
    g.vkytit(-150);
    g.title();
    g.endgrf();

    g.texmod("OFF");
    g.disfin();

    return 0;
}
