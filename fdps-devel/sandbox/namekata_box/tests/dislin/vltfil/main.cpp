/* C++ headers */
#include <iostream>
/* DISLIN header */
#include "discpp.h"

constexpr int N_lo = 5;
constexpr int N_hi = 256;

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
    // Make color tables
    ColorTable<N_lo> ctbl_lo;
    {
        double s[] = {0.0, 0.3, 0.55, 0.8, 1.0};
        double r[] = {0.0, 0.5,  1.0, 1.0, 1.0};
        double g[] = {0.0, 0.0,  0.5, 1.0, 1.0};
        double b[] = {0.0, 0.0,  0.0, 0.3, 1.0};
        ctbl_lo.set(s,r,g,b);
    }
    ColorTable<N_hi> ctbl_hi;
    ctbl_lo.interpTo(ctbl_hi);
 
    Dislin g;
    g.disini();
    // Save
    {
        g.myvlt(ctbl_hi.r, ctbl_hi.g, ctbl_hi.b, N_hi);
        const char copt[] = "SAVE";
        char filename[] = "HotMetal.tbl";
        g.vltfil(filename, copt);
    }
    // Load
    {
        const char copt[] = "LOAD";
        char filename[] = "HotMetal.tbl";
        g.vltfil(filename, copt);
    }
    // Save again to check the RGB values of the loaded color table
    {
        const char copt[] = "SAVE";
        char filename[] = "HotMetal_chk.tbl";
        g.vltfil(filename, copt);
    }

    // Finalize DISLIN
    g.disfin();

    return 0;
}
