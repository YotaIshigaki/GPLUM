#include <iostream>
#include <cmath>
#include "discpp.h"

#define SQ(x) ((x)*(x))

#define NX (256)
#define NY (256)

template <class T>
class Plotter {
private:
    Dislin g;
    T dummy;

private:
    void setFont(const int font_height,
                 const char * font_name) {
        this->g.height(font_height);
        this->g.psfont(font_name);
    }
public:
    void drawFig() {
        this->g.metafl("PS");
        int size_page[] = {2500, 2500};
        this->g.page(size_page[0], size_page[1]);
        this->g.pagmod("LAND");
        char filename[] = "test_cpp_ver2.ps";
        this->g.setfil(filename);
        this->g.disini();
        this->g.texmod("ON");

        this->g.labdig(2,"X");
        this->g.labdig(2,"Y");
        this->g.labdig(1,"Z");
        this->g.titlin("Test title",3);
        this->g.name("Test label (x)","X");
        this->g.name("Test label (y)","Y");
#if 0
        this->g.name("Test label (cb)","Z");
#else
        char label_cb[] = {"$\\log_{10} n({\\rm H}) \\;[{\\rm cm}^{-3}]$"};
        this->g.name(label_cb, "Z");
#endif
        int N_grid[] = {NX, NY};
        this->g.autres(N_grid[0], N_grid[1]);
        int xpos[] = {750, 1750};
        this->g.axspos(xpos[0],xpos[1]);
        int xlen[] = {1000, 1000, 1000};
        this->g.ax3len(xlen[0],xlen[1], xlen[2]);
        this->g.setrgb(0.0, 0.0, 0.0);
        this->setFont(50, "Helvetica-Bold");
        this->g.complx();
        double xmin[] = {0.0, 0.0, 0.0};
        double xmax[] = {1.0, 1.0, 1.0}; 
        double xtics[] = {0.25, 0.25, 0.1};
        this->g.graf3(xmin[0], xmax[0], xmin[0], xtics[0],
                      xmin[1], xmax[1], xmin[1], xtics[1],
                      xmin[2], xmax[2], xmin[2], xtics[2]);
        double data[NX*NY]; 
        int idx=0;
        for (int ix = 0; ix < N_grid[0]; ix++) {
            for (int iy = 0; iy < N_grid[1]; iy++) {
                const double pi = 4.0*std::atan(1.0);
                const double x = (double)ix / (double) N_grid[0];
                const double y = (double)iy / (double) N_grid[1];
                const double phase_x = 2.0 * pi * x;
                const double phase_y = 4.0 * pi * y;
                data[idx] = SQ(std::sin(phase_x))
                          * SQ(std::sin(phase_y));
                idx++;
            }
        }
        this->g.crvmat(data, N_grid[0], N_grid[1], 1, 1);
        int htext_title = 65;
        this->setFont(htext_title, "Palatino-Bold");
        this->g.vkytit(-150);
        this->g.title();
        this->g.endgrf();

        this->g.texmod("OFF");
        this->g.disfin();
    }
};

int main(int argc, char *argv[]) {

    Plotter<double> plotter;
    plotter.drawFig();

    return 0;
}
