#include <stdio.h>
#include <math.h>
#include "dislin.h"

#define SQ(x) ((x)*(x))

void setFont(const int font_height,
             const char * font_name) {
    height(font_height);
    psfont(font_name);
}

int main(int argc, char *argv[]) {

    metafl("PS");
    int size_page[] = {2500, 2500};
    page(size_page[0],size_page[1]);
    pagmod("LAND");
    char filename[] = "test_c.ps";
    setfil(filename);
    disini();
    texmod("ON");

    labdig(2,"X");
    labdig(2,"Y");
    labdig(1,"Z");
    titlin("Test title",3);
    name("Test label (x)","X");
    name("Test label (y)","Y");
#if 0
    name("Test label (cb)","Z");
#else
    char label_cb[] = {"$\\log_{10} n({\\rm H}) \\;[{\\rm cm}^{-3}]$"};
    name(label_cb, "Z");
#endif
    int N_grid[] = {256, 256};
    autres(N_grid[0], N_grid[1]);
    int xpos[] = {750, 1750};
    axspos(xpos[0],xpos[1]);
    int xlen[] = {1000, 1000, 1000};
    ax3len(xlen[0],xlen[1], xlen[2]);
    setrgb(0.0, 0.0, 0.0);
    setFont(50, "Helvetica-Bold");
    complx();
    double xmin[] = {0.0, 0.0, 0.0};
    double xmax[] = {1.0, 1.0, 1.0}; 
    double xtics[] = {0.25, 0.25, 0.1};
    graf3(xmin[0], xmax[0], xmin[0], xtics[0],
          xmin[1], xmax[1], xmin[1], xtics[1],
          xmin[2], xmax[2], xmin[2], xtics[2]);
    double plotted_data[256*256]; 
    int idx=0;
    for (int ix = 0; ix < N_grid[0]; ix++) {
        for (int iy = 0; iy < N_grid[1]; iy++) {
            const double pi = 4.0*atan(1.0);
            const double x = (double)ix / (double) N_grid[0];
            const double y = (double)iy / (double) N_grid[1];
            const double phase_x = 2.0 * pi * x;
            const double phase_y = 4.0 * pi * y;
            plotted_data[idx] = SQ(sin(phase_x))
                              * SQ(sin(phase_y));
            idx++;
        }
    }
    crvmat(plotted_data, N_grid[0], N_grid[1], 1, 1);
    int htext_title = 65;
    setFont(htext_title, "Palatino-Bold");
    vkytit(-150);
    title();
    endgrf();

    texmod("OFF");
    disfin();

    return 0;
}
