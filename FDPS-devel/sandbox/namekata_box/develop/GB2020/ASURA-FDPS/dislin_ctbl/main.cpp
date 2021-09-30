/* C++ headers */
#include <iostream>
/* DISLIN header */
#include "discpp.h"

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
    void checkRGBValue() const {
        for (int i = 0; i < n; i++) {
            if ((this->r[i] < 0.0) || (this->r[i] > 1.0)) {
                std::cout << "Incorrect value in r[]." << std::endl;
                std::cout << "i = " << i << " r[i] = " << this->r[i] << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if ((this->g[i] < 0.0) || (this->g[i] > 1.0)) {
                std::cout << "Incorrect value in g[]." << std::endl;
                std::cout << "i = " << i << " g[i] = " << this->g[i] << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if ((this->b[i] < 0.0) || (this->b[i] > 1.0)) {
                std::cout << "Incorrect value in b[]." << std::endl;
                std::cout << "i = " << i << " b[i] = " << this->b[i] << std::endl;
                std::exit(EXIT_FAILURE);
            }
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
    // [1] Prepare low-resolution color tables
    // Hot Metal
    ColorTable<5> ctbl_HotMetal;
    {
        double s[] = {0.0, 0.3, 0.55, 0.8, 1.0};
        double r[] = {0.0, 0.5,  1.0, 1.0, 1.0};
        double g[] = {0.0, 0.0,  0.5, 1.0, 1.0};
        double b[] = {0.0, 0.0,  0.0, 0.3, 1.0};
        ctbl_HotMetal.set(s,r,g,b);
    }

    // weird IRAF
    ColorTable<10> ctbl_weirdIRAF;
    {
        double s[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
        double r[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
        double g[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
        double b[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
        ctbl_weirdIRAF.set(s,r,g,b);
    }

    // AIPS
    ColorTable<20> ctbl_AIPS;
    {
        double s[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                      0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
        double r[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        double g[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
                      0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
        double b[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        ctbl_AIPS.set(s,r,g,b);
    }

    // Blue Gradation
    ColorTable<2> ctbl_BlueGradation;
    {
        double s[] = {0.0, 1.0};
        double r[] = {1.0, 0.0};
        double g[] = {1.0, 0.0};
        double b[] = {1.0, 1.0};
        ctbl_BlueGradation.set(s,r,g,b);
    }

    // Seismic
    ColorTable<5> ctbl_Seismic;
    {
        double s[] = {0.0, 0.15, 0.5, 0.85, 1.0};
        double r[] = {0.0,  0.0, 1.0,  1.0, 0.5};
        double g[] = {0.0,  0.0, 1.0,  0.0, 0.0};
        double b[] = {0.5,  1.0, 1.0,  0.0, 0.0};
        ctbl_Seismic.set(s,r,g,b);
    }

    // True Rainbow
    ColorTable<11> ctbl_TrueRainbow;
    {
        double s[] = {0.0, 0.1, 0.15, 0.25, 0.425, 0.5, 0.65, 0.675, 0.825, 0.925, 1.0};
        double r[] = {0.0, 0.4,  0.4,  0.0,   0.0, 0.0,  0.0,   0.5,   1.0,   1.0, 1.0};
        double g[] = {0.0, 0.0,  0.0,  0.0,  0.85, 0.9,  1.0,   1.0,   1.0,   0.5, 0.0};
        double b[] = {0.0, 0.4,  0.7,  0.9,   1.0, 0.6,  0.0,   0.0,   0.0,   0.0, 0.0};
        ctbl_TrueRainbow.set(s,r,g,b);
    }

    // Brown Gradation
    ColorTable<6> ctbl_BrownGradation;
    {
        double s[] = {0.0, 0.5,     0.7,  0.825,    0.95, 1.0};
        double r[] = {1.0, 1.0,     1.0, 0.5451,  0.5529, 0.0};
        double g[] = {1.0, 0.9, 0.62745, 0.2706, 0.22352, 0.0};
        double b[] = {1.0, 0.8,     0.0, 0.0745,  0.1529, 0.0};
        ctbl_BrownGradation.set(s,r,g,b);
    }
   
    // Green Gradation Spring-Green
    // (http://noz.day-break.net/webcolor/green.html)
    ColorTable<9> ctbl_GreenGradation_SpringGreen;
    {
        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
        double r[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
        double g[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
        double b[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
        ctbl_GreenGradation_SpringGreen.set(s,r,g,b);
    }

    // Green Gradation Teal-Green
    // (http://noz.day-break.net/webcolor/green.html)
    ColorTable<9> ctbl_GreenGradation_TealGreen;
    {
        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
        double r[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
        double g[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
        double b[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
        ctbl_GreenGradation_TealGreen.set(s,r,g,b);
    }

    // Magenta
    // (http://noz.day-break.net/webcolor/magenta.html)
    ColorTable<10> ctbl_Magenta;
    {
        double s[] = {0.0, 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 1.0};
        double r[] = {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.8,   0.6,   0.4, 0.2};
        double g[] = {1.0,   0.8,   0.6,   0.4,   0.2,   0.0,   0.0,   0.0,   0.0, 0.2};
        double b[] = {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.8,   0.6,   0.4, 0.2};
        ctbl_Magenta.set(s,r,g,b);
    }

    // Violet Magenta
    // (http://noz.day-break.net/webcolor/magenta.html)
    ColorTable<9> ctbl_VioletMagenta;
    {
        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
        double r[] = {1.0,   0.8,  0.8,   0.8, 0.8,   0.8,  0.6,   0.4, 0.2};
        double g[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.2};
        double b[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
        ctbl_VioletMagenta.set(s,r,g,b);
    }

    // Azure Blue
    // (http://noz.day-break.net/webcolor/azure.html)
    ColorTable<9> ctbl_AzureBlue;
    {
        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
        double r[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
        double g[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
        double b[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
        ctbl_AzureBlue.set(s,r,g,b);
    }

    // Gnuplot Default
    // (http://gnuplot.sourceforge.net/demo/world2.html)
    ColorTable<9> ctbl_GnuplotDefault;
    {
        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,  1.0};
        double r[] = {0.05,  0.4,  0.6,   0.8, 0.6,   0.8,  1.0,   1.0,  1.0};
        double g[] = {0.0,   0.0,  0.0,   0.0, 0.2,   0.4,  0.8,   1.0,  1.0};
        double b[] = {0.1,   0.6,  1.0,   0.6, 0.0,   0.2,  0.0,   0.0, 0.75};
        ctbl_GnuplotDefault.set(s,r,g,b);
    }

    // [2] Make high-resolution color table for DISLIN
    constexpr int N = 256;
    const char copt[] = "SAVE";
 
    Dislin g;
    // [2-1] Initialize DISLIN
    g.disini();
    // [2-2] Output color tables
    // Hot Metal
    {
        ColorTable<N> ctbl;
        ctbl_HotMetal.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "HotMetal.tbl";
        g.vltfil(filename, copt);
    }
    // weird IRAF
    {
        ColorTable<N> ctbl;
        ctbl_weirdIRAF.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "WeirdIRAF.tbl";
        g.vltfil(filename, copt);
    }
    // AIPS
    {
        ColorTable<N> ctbl;
        ctbl_AIPS.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "AIPS.tbl";
        g.vltfil(filename, copt);
    }
    // Blue Gradation
    {
        ColorTable<N> ctbl;
        ctbl_BlueGradation.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "BlueGradation.tbl";
        g.vltfil(filename, copt);
    }
    // Seismic
    {
        ColorTable<N> ctbl;
        ctbl_Seismic.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "Seismic.tbl";
        g.vltfil(filename, copt);
    }
    // True Rainbow
    {
        ColorTable<N> ctbl;
        ctbl_TrueRainbow.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "TrueRainbow.tbl";
        g.vltfil(filename, copt);
    }
    // Brown Gradation
    {
        ColorTable<N> ctbl;
        ctbl_BrownGradation.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "BrownGradation.tbl";
        g.vltfil(filename, copt);
    }
    // Green Gradation (Spring-Green)
    {
        ColorTable<N> ctbl;
        ctbl_GreenGradation_SpringGreen.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "SpringGreen.tbl";
        g.vltfil(filename, copt);
    }
    // Green Gradation (Teal-Green)
    {
        ColorTable<N> ctbl;
        ctbl_GreenGradation_TealGreen.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "TealGreen.tbl";
        g.vltfil(filename, copt);
    }
    // Magenta
    {
        ColorTable<N> ctbl;
        ctbl_Magenta.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "Magenta.tbl";
        g.vltfil(filename, copt);
    }
    // VioletMagenta
    {
        ColorTable<N> ctbl;
        ctbl_VioletMagenta.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "VioletMagenta.tbl";
        g.vltfil(filename, copt);
    }
    // Azure Blue
    {
        ColorTable<N> ctbl;
        ctbl_AzureBlue.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "AzureBlue.tbl";
        g.vltfil(filename, copt);
    }
    // Gnuplot Default
    {
        ColorTable<N> ctbl;
        ctbl_GnuplotDefault.interpTo(ctbl);
        g.myvlt(ctbl.r, ctbl.g, ctbl.b, N);
        char filename[] = "GnuplotDefault.tbl";
        g.vltfil(filename, copt);
    }

    // Finalize DISLIN
    g.disfin();

    return 0;
}
