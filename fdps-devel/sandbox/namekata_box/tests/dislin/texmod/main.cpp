#include <iostream>
#include <cmath>
#include "discpp.h"

#define SQ(x) ((x)*(x))

#define NX (256)
#define NY (256)

void setFont(Dislin & g,
             const int font_height,
             const char * font_name) {
    g.height(font_height);
    g.psfont(font_name);
}

int main(int argc, char *argv[]) {

    Dislin g;
    g.metafl("PS");
    int size_page[] = {2500, 2500};
    g.page(size_page[0],size_page[1]);
    g.pagmod("LAND");
    char filename[] = "test_cpp.ps";
    g.setfil(filename);
    g.disini();
    g.texmod("ON");

    g.labdig(2,"X");
    g.labdig(2,"Y");
    g.labdig(1,"Z");
    g.titlin("Test title",3);
    g.name("Test label (x)","X");
    g.name("Test label (y)","Y");
#if 0
    g.name("Test label (cb)","Z");
#else
    char label_cb[] = {"$\\log_{10} n({\\rm H}) \\;[{\\rm cm}^{-3}]$"};
    g.name(label_cb, "Z");
#endif
    int N_grid[] = {NX, NY};
    g.autres(N_grid[0], N_grid[1]);
    int xpos[] = {750, 1750};
    g.axspos(xpos[0],xpos[1]);
    int xlen[] = {1000, 1000, 1000};
    g.ax3len(xlen[0],xlen[1], xlen[2]);
    g.setrgb(0.0, 0.0, 0.0);
    setFont(g, 50, "Helvetica-Bold");
    g.complx();
    double xmin[] = {0.0, 0.0, 0.0};
    double xmax[] = {1.0, 1.0, 1.0}; 
    double xtics[] = {0.25, 0.25, 0.1};
    g.graf3(xmin[0], xmax[0], xmin[0], xtics[0],
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
    g.crvmat(data, N_grid[0], N_grid[1], 1, 1);
    int htext_title = 65;
    setFont(g, htext_title, "Palatino-Bold");
    g.vkytit(-150);
    g.title();
    g.endgrf();

    g.texmod("OFF");
    g.disfin();

    return 0;
}
