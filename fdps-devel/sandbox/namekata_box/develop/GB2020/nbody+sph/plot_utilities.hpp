#pragma once
/* C++ headers */
#include "common.h"
/* MPI header */
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
#include "mpi.h"
#endif
/* FDPS header */
#include <particle_simulator.hpp>
/* DISLIN header */
#include "discpp.h"
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "SPH_kernel.h"

namespace dislin_utilities {
    // [Notes]
    //     (1) we use the term `plot unit` to represent the unit system (mass, length, etc.)
    //         used for the plot. The plot unit is generally different from the simulation
    //         unit.

    enum class ScalarQuantity {NoPlot, NumberDensity, Temperature};
    enum class VectorQuantity {NoPlot, Velocity};
    enum class ColorTable {Dislin_RAIN, Dislin_SPEC, Dislin_GREY, Dislin_RRAIN,
                           Dislin_RSPEC, Dislin_RGREY, Dislin_TEMP,
                           HotMetal, weirdIRAF, AIPS, BlueGradation, Seismic,
                           TrueRainbow, BrownGradation, SpringGreen, TealGreen,
                           Magenta, VioletMagenta, AzureBlue, GnuplotDefault};

#if 0
    template <class T>
    class VectorField {
    public:
        int type; // four digits `wxyt` representing the type of vector.
        T norm_min_cgs; // the minimum of the norms of vectors in the cgs unit.
        T norm_max_cgs; // the maximum of the norms of vectors in the cgs unit.
        T maxlen; // the maximum length of vector in the plot unit.
        std::vector<T> xsta, ysta, xend, yend;
        // the start and end positions of vectors in the plot unit.
    };
#endif

    template <class T>
    class Panel {
    public: 
        /* Member variables related to the geometric position in the page  */
        int pos[2]; // the position of the panel (axis system) in the PAGE unit.
        int len[3]; // len[0:1] stores the size of the panel in the PAGE unit,
                    // while len[2] stores the length of the longer side of 
                    // the color bar in the page unit.
        T xmin[3]; // the lower limits of the corresponding axis in the plot unit.
        T xmax[3]; // the upper limits of the corresponding axis in the plot unit.
        T xtics[3]; // the intervals of the ticks.
        std::string label[4]; // 0 := the name of the X-axis,
                              // 1 := the name of the Y-axis,
                              // 2 := the name of the color bar,
                              // 3 := the tile of the plot.
        ColorTable ctbl; // Full PATH of a color table.

        /* Plot unit */
        T unit_leng;
        T unit_time;
        T unit_mass;
        T unit_velc;
        T unit_data; // the unit of plotted data 
        std::string label_leng;
        std::string label_time;
        std::string label_mass;
        std::string label_velc;

        /* Member variables related to data processing and plot */
        ScalarQuantity cplot_qty; // which quantity is plotted for the color plot.
        bool take_log10; // if true, take log10.
        bool replace_by_zero; // if true, the data value 0.0 is replaced by
                              // the pre-defined lower limit value.
        T data_minval; // the minimum value of the data in the CGS unit.
        T data_maxval; // the maximum value of the data in the CGS unit.
        T data_lowerlimit; // the lower limit of the data in the CGS unit.
        T data_upperlimit; // the upper limit of the data in the CGS unit.
        std::vector<T> data; // the flatted, two-dimensional data used for color plot.
                             // It must be given in the CGS unit.
    };

    template <class T>
    class Page {
    public:
        int n_panel; // the total number of panels in the page.
        int n_panel_x; // the number of panels along the horizontal direction.
        int n_panel_y; // the number of panels along the vertical direction.
        int size_page[2]; // the size of the page in the X,Y directions.
                          // This is determine by size_panel, n_panel_x, n_panel_y.
        int size_panel[2]; // the size of panels (`panel` is defined as color plot retion)
                           // in the page unit.
        std::vector<int> size_mgn_top,   // the size of margins around each panel in the
                         size_mgn_btm,   // page unit. In the margins, axis names, color bars,
                         size_mgn_left,  // and so on will be plotted.
                         size_mgn_right;
        std::vector< std::array<int,2> > size_lot;
            // the size of each lot, which we define as a region consists of a panel and
            // its margins. Normally, these are calculated as
            // size_lot[i][0] = size_panel[0]
            //                + size_mgn_left[i]
            //                + size_mgn_right[i];
            // size_lot[i][1] = size_panel[1]
            //                + size_mgn_btm[i]
            //                + size_mgn_top[i];
        std::vector< std::array<int,2> > pos_lot;
            // the position of the left bottom corner of each lot in the page unit.
        int npx_page[2];
            // the resolution of figure; only used for the cases that the save format
            // is an image file format.
        Panel<T> panel; // this data structure is shared by all the panels

        Page() {
            npx_page[0] = 0;
            npx_page[1] = 0;
        }
    };

    template <class T>
    class Mesh {
    public:
        // The variables below defines the configuration of plotted region
        // in the simulation box. The plotted region is defined as
        // a rectangle region in the x'''-y'''-plane of the x'''y'''z'''
        // coordinate, which is obtained by the coordinate rotation using
        // the Euler angles (\alpha, \beta, \gamma) with `origin` being 
        // the origin of the rotation.
        // All the values must be given in the simulation unit.
        PS::Vector3<T> origin; // the origin of plot region
        T alpha, beta, gamma; // the Euler angles (z-x'-z'' rotation)
        PS::Orthotope2<T> bbox; // the boundary box size on the x'''-y''' plane.

        // The variables below defines two-dimensional grid points, 
        // which are used for plot.
        int nx,ny;
        T dx, dy;
        std::vector<T> x,y; 

        template <class U>
        void define(const PS::Vector3<U> origin,
                    const U alpha, const U beta, const U gamma, // Euler angles
                    const PS::Orthotope2<U> bbox,
                    const int nx, const int ny) {
            // Copy the arguments to the member variables
            this->origin = origin;
            this->alpha = alpha;
            this->beta = beta;
            this->gamma = gamma;
            this->bbox = bbox;
            this->nx = nx;
            this->ny = ny;

            // Define x''' coordinate
            {
                this->x.resize(nx);
                const T length = (bbox.high_.x - bbox.low_.x);
                const T center = 0.5 * (bbox.low_.x + bbox.high_.x);
                const T dx = length / static_cast<T>(nx);
                const int cx = nx/2;
                this->dx = dx;
                if (nx % 2 == 1) { // # of the grid points is odd
                    for (int i = 0; i < nx; i++) {
                        if (i < cx) this->x[i] = center - dx * (cx - i);
                        else this->x[i] = center + dx * (i - cx);
                    }
                } else { // # of the grid points is even
                    for (int i = 0; i < nx; i++) {
                        if (i < cx) this->x[i] = center + dx * i - 0.5 * (length - dx);
                        else this->x[i] = center - dx * (nx - 1 - i) + 0.5 * (length - dx);
                    }
                }
            }

            // Define y''' coordinate
            {
                this->y.resize(ny);
                const T length = (bbox.high_.y - bbox.low_.y);
                const T center = 0.5 * (bbox.low_.y + bbox.high_.y);
                const T dy = length / static_cast<T>(ny);
                const int cy = ny/2;
                this->dy = dy;
                if (ny % 2 == 1) { // # of the grid points is odd
                    for (int i = 0; i < ny; i++) {
                        if (i < cy) this->y[i] = center - dy * (cy - i);
                        else this->y[i] = center + dy * (i - cy);
                    }
                } else { // # of the grid points is even
                    for (int i = 0; i < ny; i++) {
                        if (i < cy) this->y[i] = center + dy * i - 0.5 * (length - dy);
                        else this->y[i] = center - dy * (ny - 1 - i) + 0.5 * (length - dy);
                    }
                }
            }

#if 0
            // Check
            if (PS::Comm::getRank() == 0) {
                std::ofstream ofs;
                ofs.open("mesh.txt", std::ios::trunc);
                ofs << "origin = " << this->origin << std::endl;
                ofs << "alpha = " << this->alpha << std::endl;
                ofs << "beta  = " << this->beta << std::endl;
                ofs << "gamma = " << this->gamma << std::endl;
                ofs << "bbox = " << this->bbox << std::endl;
                ofs << "nx = " << this->nx << std::endl;
                ofs << "ny = " << this->ny << std::endl;
                ofs << "dx = " << this->dx << std::endl;
                ofs << "dy = " << this->dy << std::endl;
                ofs << "Content of (i, x[i]):" << std::endl;
                for (int i = 0; i < this->nx; i++)
                    ofs << i << "    " << this->x.at(i) << std::endl;
                ofs << "Content of (j, y[j]):" << std::endl;
                for (int j = 0; j < this->ny; j++)
                    ofs << j << "    " << this->y.at(j) << std::endl;
                ofs.close();
            }
#endif

        }

        template <class U>
        PS::Vector3<U> getPosTransformed(const PS::Vector3<U> & pos_ptcl) const {
            // This returns the position of a particle in the dislin coordinate.
            const T x = pos_ptcl.x - this->origin.x;
            const T y = pos_ptcl.y - this->origin.y;
            const T z = pos_ptcl.z - this->origin.z;
            const T a = this->alpha;
            const T b = this->beta;
            const T g = this->gamma;
            PS::Vector3<U> pos;
            pos.x = x * (  std::cos(g) * std::cos(a) - std::sin(g) * std::cos(b) * std::sin(a))
                  + y * (  std::cos(g) * std::sin(a) + std::sin(g) * std::cos(b) * std::cos(a))
                  + z * (  std::sin(g) * std::sin(b));
            pos.y = x * (- std::sin(g) * std::cos(a) - std::cos(g) * std::cos(b) * std::sin(a))
                  + y * (- std::sin(g) * std::sin(a) + std::cos(g) * std::cos(b) * std::cos(a))
                  + z * (  std::cos(g) * std::sin(b));
            pos.z = x * (  std::sin(b) * std::sin(a))
                  + y * (- std::sin(b) * std::cos(a))
                  + z * (  std::cos(b));
            return pos;
        }

        PS::Vector3<T> getPosOriginal(const int i, const int j) const {
            // This returns the position of a grid point (i,j) in the simulation coordinate.
            const T x = this->x.at(i); // x'''
            const T y = this->y.at(j); // y'''
            const T a = this->alpha;
            const T b = this->beta;
            const T g = this->gamma;
            PS::Vector3<T> pos;
            pos.x = x * (  std::cos(a) * std::cos(g) - std::sin(a) * std::cos(b) * std::sin(g))
                  + y * (- std::cos(a) * std::sin(g) - std::sin(a) * std::cos(b) * std::cos(g));
            pos.y = x * (  std::sin(a) * std::cos(g) + std::cos(a) * std::cos(b) * std::sin(g))
                  + y * (- std::sin(a) * std::sin(g) + std::cos(a) * std::cos(b) * std::cos(g));
            pos.z = x * (  std::sin(b) * std::sin(g))
                  + y * (  std::sin(b) * std::cos(g));
            return this->origin + pos;
        }

        template <class U>
        PS::S32ort2 getInteractionRegion(const PS::Vector3<U> & pos_ptcl,
                                         const U h) const {
            const PS::Vector3<T> pos = this->getPosTransformed(pos_ptcl);
            PS::S32ort2 box;
            if (std::fabs(pos.z) <= h) {
                const U safety = 1.001;
                box.low_.x  = std::max(0, static_cast<int>((pos.x - safety * h - this->getXMin()) / this->dx));
                box.low_.y  = std::max(0, static_cast<int>((pos.y - safety * h - this->getYMin()) / this->dy));
                box.high_.x = std::min(static_cast<int>(this->x.size() - 1),
                                       static_cast<int>((pos.x + safety * h - this->getXMin()) / this->dx));
                box.high_.y = std::min(static_cast<int>(this->y.size() - 1),
                                       static_cast<int>((pos.y + safety * h - this->getYMin()) / this->dy));
            }
            return box;
        }

        int getNumberOfMeshAlongXAxis() const {
            return this->nx;
        }
        int getNumberOfMeshAlongYAxis() const {
            return this->ny;
        }
        T getXMin() const {
            return this->x.front();
        }
        T getXMax() const {
            return this->x.back();
        }
        T getYMin() const {
            return this->y.front();
        }
        T getYMax() const {
            return this->y.back();
        }
    };

    template <class T>
    class RGBTable {
    public:
        int n;
        std::vector<T> s; // position between [0,1]
        std::vector<T> r; // R of RGB
        std::vector<T> g; // G of RGB
        std::vector<T> b; // B of RGB
    
        int size() const { return this->n; }

        template <class U> // float or double
        void set(const int n,
                 const U s[],
                 const U r[],
                 const U g[],
                 const U b[]) {
            this->n = n;
            this->s.resize(n);
            this->r.resize(n);
            this->g.resize(n);
            this->b.resize(n);
            for (int i = 0; i < n; i++) {
                this->s[i] = s[i];
                this->r[i] = r[i];
                this->g[i] = g[i];
                this->b[i] = b[i];
            }
        }

        void check() const {
            for (int i = 0; i < this->size(); i++) {
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

        template <class U>
        void interpTo(RGBTable<U> & dest, const int n_dest) const {
            // Reinitialize dest
            dest.n = n_dest;
            dest.s.resize(n_dest);
            dest.r.resize(n_dest);
            dest.g.resize(n_dest);
            dest.b.resize(n_dest);
            // Interpolation
            const T ds = 1.0/static_cast<T>(dest.size()-1);
            for (int i = 0; i < dest.size(); i++) {
                const T s = ds * i;
                dest.s[i] = s;
                // Binary search to find the first array index at which s[] > s
                typename std::vector<T>::const_iterator it
                    = std::lower_bound(this->s.begin(), this->s.end(), s);
                int j;
                if (it != this->s.end()) j = std::distance(this->s.begin(), it) - 1;
                else j = this->size() - 1; 
                // Compute dest.r[i], dest.g[i], dest.b[i]
                if ((0 <= j) && (j <= this->size() - 2)) {
                    const T offset = s - this->s[j];
                    const T intvl = this->s[j+1] - this->s[j];
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

    template <class T>
    class Plotter {
    private:
        Dislin g;

        std::string save_fmt;
        std::string filename;

        Page<T> page;
        Mesh<T> mesh;
    public:

        Plotter() : save_fmt(""), filename("") {}

    private:
        bool isSupportedFormat(const std::string save_fmt) {
            if ((save_fmt == "GKSL") ||
                (save_fmt == "CGM")  ||
                (save_fmt == "PS")   ||
                (save_fmt == "EPS")  ||
                (save_fmt == "PDF")  ||
                (save_fmt == "HPGL") ||
                (save_fmt == "SVG")  ||
                (save_fmt == "IPE")  ||
                (save_fmt == "WMF")  ||
                (save_fmt == "GIF")  ||
                (save_fmt == "TIFF") ||
                (save_fmt == "PNG")  ||
                (save_fmt == "PPM")  ||
                (save_fmt == "IMAG") ||
                (save_fmt == "BMP")  ||
                (save_fmt == "VIRT") ||
                (save_fmt == "CONS") ||
                (save_fmt == "XWIN") ||
                (save_fmt == "GL")  ) {
                return true;
            } else {
                return false;
            }
        }

        bool isOutputDeviceImageFile() {
            assert(this->save_fmt != "");
            if ((this->save_fmt == "GIF")  ||
                (this->save_fmt == "TIFF") ||
                (this->save_fmt == "PNG")  ||
                (this->save_fmt == "PPM")  ||
                (this->save_fmt == "IMAG") ||
                (this->save_fmt == "BMP")) {
                return true;
            } else {
                return false;
            }
        }

        bool isOutputDeviceVectorFile() {
            assert(this->save_fmt != "");
            if ((this->save_fmt == "PS")  ||
                (this->save_fmt == "EPS") ||
                (this->save_fmt == "PDF") ||
                (this->save_fmt == "SVG")) {
                return true;
            } else {
                return false;
            }
        }

        void initDislin(const double npx_in_page_unit[2]) {
            this->g.metafl(this->save_fmt.c_str());
            this->g.page(this->page.size_page[0],
                         this->page.size_page[1]);
            this->g.pagmod("NONE");
            this->g.setfil(this->filename.c_str());
            if (this->isOutputDeviceImageFile()) {
                // If the format is an image file,
                // the background and foreground colors
                // are set to white and black, respectively.
                this->g.scrmod("REVERS");
                // Set the resulution of figure
                const int npx_page[2]
                    = {npx_in_page_unit[0] * this->page.size_page[0],
                       npx_in_page_unit[1] * this->page.size_page[1]};
                this->g.winsiz(npx_page[0], npx_page[1]);
                std::cout << "Resolution: "
                          << npx_page[0] << " (X), "
                          << npx_page[1] << " (Y) (in pixels)."
                          << std::endl;
                // Set the coloring mode for TIFF, PNG, BMP, and IMAGE
                this->g.imgfmt("RGB");
            }
            this->g.disini();
        }

        void setColorTable(ColorTable tbl_name) {
            RGBTable<T> tbl_raw;
            switch (tbl_name) {
                case ColorTable::Dislin_RAIN:
                    this->g.setvlt("RAIN");
                    return;
                case ColorTable::Dislin_SPEC:
                    this->g.setvlt("SPEC");
                    return;
                case ColorTable::Dislin_GREY:
                    this->g.setvlt("GREY");
                    return;
                case ColorTable::Dislin_RRAIN:
                    this->g.setvlt("RRAIN");
                    return;
                case ColorTable::Dislin_RSPEC:
                    this->g.setvlt("RSPEC");
                    return;
                case ColorTable::Dislin_RGREY:
                    this->g.setvlt("RGREY");
                    return;
                case ColorTable::Dislin_TEMP:
                    this->g.setvlt("TEMP");
                    return;
                case ColorTable::HotMetal:
                    {
                        double s[] = {0.0, 0.3, 0.55, 0.8, 1.0};
                        double r[] = {0.0, 0.5,  1.0, 1.0, 1.0};
                        double g[] = {0.0, 0.0,  0.5, 1.0, 1.0};
                        double b[] = {0.0, 0.0,  0.0, 0.3, 1.0};
                        tbl_raw.set(5,s,r,g,b);
                    }
                    break;
                case ColorTable::weirdIRAF:
                    {
                        double s[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
                        double r[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
                        double g[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
                        double b[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
                        tbl_raw.set(10,s,r,g,b);
                    }
                    break;
                case ColorTable::AIPS:
                    {
                        double s[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                                      0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
                        double r[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 
                                      0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
                        double g[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
                                      0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
                        double b[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                        tbl_raw.set(20,s,r,g,b);
                    }
                    break;
                case ColorTable::BlueGradation:
                    {
                        double s[] = {0.0, 1.0};
                        double r[] = {1.0, 0.0};
                        double g[] = {1.0, 0.0};
                        double b[] = {1.0, 1.0};
                        tbl_raw.set(2,s,r,g,b);
                    }
                    break;
                case ColorTable::Seismic:
                    {
                        double s[] = {0.0, 0.15, 0.5, 0.85, 1.0};
                        double r[] = {0.0,  0.0, 1.0,  1.0, 0.5};
                        double g[] = {0.0,  0.0, 1.0,  0.0, 0.0};
                        double b[] = {0.5,  1.0, 1.0,  0.0, 0.0};
                        tbl_raw.set(5,s,r,g,b);
                    }
                    break;
                case ColorTable::TrueRainbow:
                    {
                        double s[] = {0.0, 0.1, 0.15, 0.25, 0.425, 0.5, 0.65, 0.675, 0.825, 0.925, 1.0};
                        double r[] = {0.0, 0.4,  0.4,  0.0,   0.0, 0.0,  0.0,   0.5,   1.0,   1.0, 1.0};
                        double g[] = {0.0, 0.0,  0.0,  0.0,  0.85, 0.9,  1.0,   1.0,   1.0,   0.5, 0.0};
                        double b[] = {0.0, 0.4,  0.7,  0.9,   1.0, 0.6,  0.0,   0.0,   0.0,   0.0, 0.0};
                        tbl_raw.set(11,s,r,g,b);
                    }
                    break;
                case ColorTable::BrownGradation:
                    {
                        double s[] = {0.0, 0.5,     0.7,  0.825,    0.95, 1.0};
                        double r[] = {1.0, 1.0,     1.0, 0.5451,  0.5529, 0.0};
                        double g[] = {1.0, 0.9, 0.62745, 0.2706, 0.22352, 0.0};
                        double b[] = {1.0, 0.8,     0.0, 0.0745,  0.1529, 0.0};
                        tbl_raw.set(6,s,r,g,b);
                    }
                    break;
                case ColorTable::SpringGreen:
                    // (http://noz.day-break.net/webcolor/green.html)
                    {
                        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
                        double r[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
                        double g[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
                        double b[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
                        tbl_raw.set(9,s,r,g,b);
                    }
                    break;
                case ColorTable::TealGreen:
                    // (http://noz.day-break.net/webcolor/green.html)
                    {
                        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
                        double r[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
                        double g[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
                        double b[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
                        tbl_raw.set(9,s,r,g,b);
                    }
                    break;
                case ColorTable::Magenta:
                    // (http://noz.day-break.net/webcolor/magenta.html)
                    {
                        double s[] = {0.0, 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 1.0};
                        double r[] = {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.8,   0.6,   0.4, 0.2};
                        double g[] = {1.0,   0.8,   0.6,   0.4,   0.2,   0.0,   0.0,   0.0,   0.0, 0.2};
                        double b[] = {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.8,   0.6,   0.4, 0.2};
                        tbl_raw.set(10,s,r,g,b);
                    }
                    break;
                case ColorTable::VioletMagenta:
                    // (http://noz.day-break.net/webcolor/magenta.html)
                    {
                        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
                        double r[] = {1.0,   0.8,  0.8,   0.8, 0.8,   0.8,  0.6,   0.4, 0.2};
                        double g[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.2};
                        double b[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
                        tbl_raw.set(9,s,r,g,b);
                    }
                    break;
                case ColorTable::AzureBlue:
                    {
                        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
                        double r[] = {1.0,   0.6,  0.4,   0.2, 0.0,   0.0,  0.0,   0.0, 0.0};
                        double g[] = {1.0,   0.8,  0.6,   0.4, 0.2,   0.0,  0.0,   0.0, 0.0};
                        double b[] = {1.0,   1.0,  1.0,   1.0, 1.0,   0.8,  0.6,   0.4, 0.2};
                        tbl_raw.set(9,s,r,g,b);
                    }
                    break;
                case ColorTable::GnuplotDefault: 
                    // (http://gnuplot.sourceforge.net/demo/world2.html)
                    {
                        double s[] = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,  1.0};
                        double r[] = {0.05,  0.4,  0.6,   0.8, 0.6,   0.8,  1.0,   1.0,  1.0};
                        double g[] = {0.0,   0.0,  0.0,   0.0, 0.2,   0.4,  0.8,   1.0,  1.0};
                        double b[] = {0.1,   0.6,  1.0,   0.6, 0.0,   0.2,  0.0,   0.0, 0.75};
                        tbl_raw.set(9,s,r,g,b);
                    }
                    break;
                default:
                    assert(false);
            }
            // Make high-resolution color table and use it
            constexpr int N = 256;
            RGBTable<T> tbl;
            tbl_raw.interpTo(tbl, N);
            this->g.myvlt(tbl.r.data(), tbl.g.data(), tbl.b.data(), N);
        }
        void setFont(const int font_height,
                     const std::string font_name) {
            this->g.height(font_height);
            if (this->isOutputDeviceVectorFile()) {
                this->g.psfont(font_name.c_str());
            }
        }
        void setColorBarProp(const T xmin,
                             const T xmax,
                             const int ntics,
                             const int len,
                             const ColorTable ctbl,
                             const std::string label) {
            this->page.panel.xmin[2] = xmin;
            this->page.panel.xmax[2] = xmax;
            this->page.panel.xtics[2] = (xmax-xmin)/static_cast<T>(ntics-1);
            this->page.panel.len[2] = len;
            this->page.panel.ctbl = ctbl;
            this->page.panel.label[2] = label;
        }
    public:
        /* Member functions */
        void setPageLayoutAndOpenPage(const int nrows,
                                      const int ncols,
                                      const int size_panel[2],
                                      const double npx_in_page_unit[2],
                                      const char * save_fmt_in_char,
                                      const char * filename_in_char) {
            const std::string save_fmt(save_fmt_in_char);
            const std::string filename(filename_in_char);
            this->setPageLayoutAndOpenPage(nrows, ncols,
                                           size_panel, npx_in_page_unit,
                                           save_fmt, filename);
        }
        void setPageLayoutAndOpenPage(const int nrows,
                                      const int ncols,
                                      const int size_panel[2],
                                      const double npx_in_page_unit[2],
                                      const std::string save_fmt,
                                      const std::string filename) {
            // Error check
            if (!this->isSupportedFormat(save_fmt)) {
                if (PS::Comm::getRank() == 0) {
                    std::cout << "Save format " << save_fmt << " is not supported in DISLIN." << std::endl;
                }
                PS::Abort(-1);
                std::exit(EXIT_FAILURE);
            }

            // Set the configuration of lots
            assert(nrows > 0);
            assert(ncols > 0);
            this->page.n_panel   = nrows * ncols;
            this->page.n_panel_x = ncols;
            this->page.n_panel_y = nrows;
            
            // Set the sizes of each panel and margins
            this->page.size_mgn_top.resize(this->page.n_panel);
            this->page.size_mgn_btm.resize(this->page.n_panel);
            this->page.size_mgn_left.resize(this->page.n_panel);
            this->page.size_mgn_right.resize(this->page.n_panel);
            this->page.size_lot.resize(this->page.n_panel);
            this->page.pos_lot.resize(this->page.n_panel);
            assert(size_panel[0] > 0);
            assert(size_panel[1] > 0);
            this->page.size_panel[0] = size_panel[0];
            this->page.size_panel[1] = size_panel[1];
            for (int j = 0; j < nrows; j++) {
                for (int i = 0; i < ncols; i++) {
                    const int idx = i + ncols * j;
                    // Set top margin
                    if (idx < ncols) this->page.size_mgn_top[idx] = 500; // at the top row
                    else this->page.size_mgn_top[idx] = 75;
                    // Set bottom margin
                    if (idx >= ncols*(nrows-1) ) this->page.size_mgn_btm[idx] = 500; // at the bottom row
                    else this->page.size_mgn_btm[idx] = 75;
                    // Set right margin
                    if ((idx+1) % ncols == 0) this->page.size_mgn_right[idx] = 1000; // in the rightmost column
                    else this->page.size_mgn_right[idx] = 300; // 100 org.
                    // Set left margin
                    if (idx % ncols == 0) this->page.size_mgn_left[idx] = 500;
                    else this->page.size_mgn_left[idx] = 300; // 100 org.
                    // Set the size of a lot
                    this->page.size_lot[idx][0] = this->page.size_panel[0]
                                                + this->page.size_mgn_left[idx]
                                                + this->page.size_mgn_right[idx];
                    this->page.size_lot[idx][1] = this->page.size_panel[1]
                                                + this->page.size_mgn_top[idx]
                                                + this->page.size_mgn_btm[idx];
                }
            }
            // Set the size of page and the positions of lots
            this->page.size_page[0] = 0;
            this->page.size_page[1] = 0;
            for (int j = 0; j < nrows; j++) {
                for (int i = 0; i < ncols; i++) {
                    const int idx = i + ncols*j;
                    // Position & page size (X)
                    if (idx % ncols == 0) { // In the leftmost column
                        this->page.pos_lot[idx][0] = 0;
                        this->page.size_page[1] += this->page.size_lot[idx][1];
                    } else {
                       this->page.pos_lot[idx][0] = this->page.pos_lot[idx-1][0]
                                                  + this->page.size_lot[idx-1][0];
                    }
                    // Position & page size (Y)
                    if (idx < ncols) { // In the top row
                       this->page.pos_lot[idx][1] = this->page.size_lot[idx][1];
                       this->page.size_page[0] += this->page.size_lot[idx][0];
                    } else {
                       this->page.pos_lot[idx][1] = this->page.pos_lot[idx-ncols][1]
                                                  + this->page.size_lot[idx][1];
                    }
                }
            }
            if (PS::Comm::getRank() == 0) {
                std::cout << "Page sizes: "
                          << this->page.size_page[0] << " (X), "
                          << this->page.size_page[1] << " (Y)"
                          << std::endl;
            }

            // Initialize Dislin
            this->save_fmt = save_fmt;
            this->filename = filename;
            if (PS::Comm::getRank() == 0) {
                this->initDislin(npx_in_page_unit);
            }

        }

        template <class U>
        void setMesh(const PS::Vector3<U> origin,
                     const U alpha, const U beta, const U gamma,
                     const PS::Orthotope2<U> bbox,
                     const int nx, const int ny) {
            this->mesh.define(origin, alpha, beta, gamma, bbox, nx, ny);
        }

        void setPanelConfig(const int panel_id,
                            const ScalarQuantity cplot_qty) { 
            // Error check
            assert(0<= panel_id && panel_id < this->page.n_panel);

            // Set the physical quantities to be plotted
            this->page.panel.cplot_qty   = cplot_qty;
     
            // Set the plot unit and corresponding labels
            this->page.panel.unit_leng = phys_const::kpc;
            this->page.panel.unit_time = phys_const::Myr;
            this->page.panel.unit_mass = phys_const::Msolar;
            this->page.panel.unit_velc = phys_const::km;
            this->page.panel.label_leng = " [kpc]";
            this->page.panel.label_time = " [Myr]";
            this->page.panel.label_mass = " [Msun]";
            this->page.panel.label_velc = " [km/s]";

            // Set the unit of data
            switch (cplot_qty) {
                case ScalarQuantity::NumberDensity:
                    this->page.panel.unit_data = 1.0; // i.e. [cm^{-3}]
                    break;
                case ScalarQuantity::Temperature:
                    this->page.panel.unit_data = 1.0; // i.e. [K]
                    break;
                default:
                    this->page.panel.unit_data = 1.0;
            }

            // Calculate the position of a target panel
            this->page.panel.pos[0] = this->page.pos_lot[panel_id][0]
                                    + this->page.size_mgn_left[panel_id];
            this->page.panel.pos[1] = this->page.pos_lot[panel_id][1]
                                    - this->page.size_mgn_btm[panel_id];
            this->page.panel.len[0] = this->page.size_panel[0];
            this->page.panel.len[1] = this->page.size_panel[1];


            // Set the resolution of the plot
            const int nx = mesh.getNumberOfMeshAlongXAxis(); 
            const int ny = mesh.getNumberOfMeshAlongYAxis();
            const int ntot = nx * ny;
            this->page.panel.data.resize(ntot);
            for (int i = 0; i < ntot; i++)
                this->page.panel.data[i] = 0.0;
            
            // Set 
            int ntics[] = {4, 4};
            const T xlen = mesh.getXMax() - mesh.getXMin();
            const T ylen = mesh.getYMax() - mesh.getYMin();
            const T coeff = (run_param::unit::leng/this->page.panel.unit_leng);
            this->page.panel.xmin[0]  = mesh.getXMin() * coeff;
            this->page.panel.xmax[0]  = mesh.getXMax() * coeff;
            this->page.panel.xmin[1]  = mesh.getYMin() * coeff;
            this->page.panel.xmax[1]  = mesh.getYMax() * coeff;
            this->page.panel.xtics[0] = xlen * coeff / ntics[0];
            this->page.panel.xtics[1] = ylen * coeff / ntics[1];
            this->page.panel.label[0] = "x" + this->page.panel.label_leng;
            this->page.panel.label[1] = "y" + this->page.panel.label_leng;

        }

        template <class Tpsys>
        void setupPlotData(const Tpsys & psys) {
            const PS::S32 nx = this->mesh.getNumberOfMeshAlongXAxis();
            const PS::S32 ny = this->mesh.getNumberOfMeshAlongYAxis();
            const PS::S32 ntot = nx * ny;
            // Initialze data and normalization constants
            bool do_normalize {true};
            std::vector<T> data_loc, wsum_loc, wsum_glb;
            data_loc.resize(ntot);
            wsum_loc.resize(ntot);
            wsum_glb.resize(ntot);
            for (auto & val : data_loc) val = 0.0;
            for (auto & val : wsum_loc) val = 0.0;
            // Calculate plotted data 
            switch (this->page.panel.cplot_qty) {
                // Here we assumed that the density is already calculated somewhere else.
                case ScalarQuantity::NumberDensity:
                    do_normalize = false;
                    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
                        const PS::F64 m_i = psys[i].getMass();
                        const PS::F64vec pos_i = psys[i].getPos();
                        const PS::F64 h_i = psys[i].getKernelSupportRadius();
                        PS::S32ort2 rect = this->mesh.getInteractionRegion(pos_i, h_i);
                        if (rect.isValid()) {
                            for (PS::S32 ii = rect.low_.x; ii <= rect.high_.x; ii++) {
                                for (PS::S32 jj = rect.low_.y; jj <= rect.high_.y; jj++) {
                                    const PS::F64vec pos_j = this->mesh.getPosOriginal(ii,jj);
                                    const PS::F64vec dr = pos_i - pos_j;
                                    const PS::F64 r = std::sqrt(dr * dr);
                                    const PS::S32 idx = ii * ny + jj;
                                    const PS::F64 coeff = run_param::unit::dens 
                                                        * run_param::ism::Xhydrogen
                                                        / phys_const::Mhydrogen;
                                    const PS::F64 part = m_i * W(r, h_i) * coeff;
                                    data_loc[idx] += m_i * W(r, h_i) * coeff;
                                }
                            }
                        }
                    }
                    break;
                case ScalarQuantity::Temperature:
                    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
                        const PS::F64 m_i = psys[i].getMass();
                        const PS::F64vec pos_i = psys[i].getPos();
                        const PS::F64 h_i = psys[i].getKernelSupportRadius();
                        const PS::F64 dens_i = psys[i].getMassDensity();
                        const PS::F64 temp_i = psys[i].getTemperature();
                        PS::S32ort2 rect = this->mesh.getInteractionRegion(pos_i, h_i);
                        if (rect.isValid()) {
                            for (PS::S32 ii = rect.low_.x; ii <= rect.high_.x; ii++) {
                                for (PS::S32 jj = rect.low_.y; jj <= rect.high_.y; jj++) {
                                    const PS::F64vec pos_j = this->mesh.getPosOriginal(ii,jj);
                                    const PS::F64vec dr = pos_i - pos_j;
                                    const PS::F64 r = std::sqrt(dr * dr);
                                    const PS::F64 w = m_i * W(r, h_i) / dens_i;
                                    const PS::S32 idx = ii * ny + jj;
                                    const PS::F64 coeff = run_param::unit::temp;
                                    data_loc[idx] += temp_i * w * coeff;
                                    wsum_loc[idx] += w;
                                }
                            }
                        }
                    }
                    break;
            }

            // Reduction
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
            if (do_normalize) {
                MPI_Reduce(wsum_loc.data(), wsum_glb.data(), ntot,
                           PS::GetDataType<T>(), MPI_SUM, 0,
                           MPI_COMM_WORLD);
            }
            MPI_Reduce(data_loc.data(), this->page.panel.data.data(), ntot,
                       PS::GetDataType<T>(), MPI_SUM, 0,
                       MPI_COMM_WORLD);
#else
            if (do_normalize) {
                for (PS::S32 i = 0; i < ntot; i++) {
                    if (wsum_loc[i] > 0.0) this->page.panel.data[i] = data_loc[i] / wsum_loc[i];
                    else this->page.panel.data[i] = data_loc[i];
                }
            } else {
                for (PS::S32 i = 0; i < ntot; i++) this->page.panel.data[i] = data_loc[i];
            }
#endif

#if 0
            // Check
            std::ofstream ofs;
            ofs.open("check.txt", std::ios::trunc);
            for (PS::S32 i = 0; i < nx; i++) {
                for (PS::S32 j = 0; j < ny; j++) {
                    const PS::S32 idx = i * ny + j;
                    ofs << i << "    "
                        << j << "    "
                        << this->page.panel.data[idx] << std::endl;
                }
            }
            ofs.close();
#endif

        }

        void drawPanel(const int panel_id,
                       const std::string panel_text) {
            // Common settings
            const int htext_title = 65;
            this->g.reset("ALL"); 
        
            // Calculate minimum/maximum
            this->page.panel.data_minval =  std::numeric_limits<T>::max();
            this->page.panel.data_maxval = -std::numeric_limits<T>::max();
            for (auto val : this->page.panel.data) {
                if (this->page.panel.data_minval > val) this->page.panel.data_minval = val;
                if (this->page.panel.data_maxval < val) this->page.panel.data_maxval = val;
            }
        
            // Set lower limit/upper limit (CGS unit)
            switch (this->page.panel.cplot_qty) {
                case ScalarQuantity::NumberDensity:
                    this->page.panel.take_log10 = true;
                    this->page.panel.data_lowerlimit = 1.0e-4;
                    this->page.panel.data_upperlimit = 1.0e4;
                    break;
                case ScalarQuantity::Temperature:
                    this->page.panel.take_log10 = true;
                    this->page.panel.data_lowerlimit = 1.0e0;
                    this->page.panel.data_upperlimit = 1.0e6;
                    break;
            }
        
            // Fill data with data_lowerlimit/data_upperlimit (Option)
            for (auto & val : this->page.panel.data) {
                if (val < this->page.panel.data_lowerlimit) val = this->page.panel.data_lowerlimit;
                if (val > this->page.panel.data_upperlimit) val = this->page.panel.data_upperlimit;
                if (this->page.panel.take_log10) val = std::log10(val/this->page.panel.unit_data);
                else val = val/this->page.panel.unit_data;
            }
  
            // Setup labels, etc. 
            this->g.texmod("ON");
            int ntics;
            std::string label_tmp;
            switch (this->page.panel.cplot_qty) {
                case ScalarQuantity::NumberDensity:
                    this->page.panel.label[3] = "Hydrogen number density";
                    label_tmp = "$\\log_{10} n({\\rm H}) \\;[{\\rm cm}^{-3}]$";
                    ntics = 9;
                    break;
                case ScalarQuantity::Temperature:
                    this->page.panel.label[3] = "Gas temperature";
                    label_tmp = "$\\log_{10} T_{\\rm gas} \\;[{\\rm K}]$";
                    ntics = 7;
                    break;
            }
  
            // Additional processes 
            T xmin, xmax;
            if (this->page.panel.take_log10) {
               xmin = std::log10(this->page.panel.data_lowerlimit
                                /this->page.panel.unit_data);
               xmax = std::log10(this->page.panel.data_upperlimit
                                /this->page.panel.unit_data);
            } else {
               xmin = this->page.panel.data_lowerlimit/this->page.panel.unit_data;
               xmax = this->page.panel.data_upperlimit/this->page.panel.unit_data;
            }
            const int len = this->page.panel.len[1];
            ColorTable ctbl_name;
            switch (this->page.panel.cplot_qty) {
                case ScalarQuantity::Temperature:
                    ctbl_name = ColorTable::HotMetal;
                    setColorBarProp(xmin, xmax, ntics, len, ctbl_name, label_tmp);
                    break;
                default: 
                    ctbl_name = ColorTable::TrueRainbow;
                    setColorBarProp(xmin, xmax, ntics, len, ctbl_name, label_tmp);
                    break;
            }

            // Draw panel
            this->g.labdig(2,"X");
            this->g.labdig(2,"Y");
            this->g.labdig(1,"Z");
#if 0
            if ((panel_id + 1) % this->page.n_panel_x != 0) {
                // We do not plot the color bar except for
                // the panels at the rightmost column.
                this->g.nobar();
            }
#endif
            // We plot X-label only for the panels at the bottom raw.
            std::string btm_x_axis;
            if (panel_id >= this->page.n_panel_x * (this->page.n_panel_y - 1)) {
                this->g.name(this->page.panel.label[0].c_str(), "X");
                btm_x_axis = "NAME";
            } else {
                this->g.name("","X");
                btm_x_axis = "TICKS";
            }
            // We plot Y-label only for the panels at the leftmost column
            std::string left_y_axis;
            if (panel_id % this->page.n_panel_x == 0) {
               this->g.name(this->page.panel.label[1].c_str(), "Y");
               left_y_axis = "NAME";
            } else {
               this->g.name("", "Y");
               left_y_axis = "TICKS";
            }
            this->g.name(this->page.panel.label[2].c_str(), "Z");
            this->g.namdis(60,"Z"); 
            this->g.hname(65);
            this->g.setgrf(btm_x_axis.c_str(), left_y_axis.c_str(), "TICKS", "TICKS");
            // [Suppl.] ------------------------------------------------------------
            //   (1) NAMEDIS() sets the distance between axis names and labels.
            //   (2) HNAME() defines the character height for axis names.
            //   (3) SETGRF() removes a part of an axis or a complete axis
            //       from an axis system.
            //
            //       The call is:   CALL SETGRF(C1, C2, C3, C4)  # level 1, 2, 3
            //
            //       Ci are character strings corresponding to the four axes of
            //       an axis system. C1 corresponds to the lower X-axis,
            //       C2 to the left Y-axis, C3 to the upper X-axis and C4 to the
            //       right Y-axis. The parameters can have the values 'NONE',
            //       'LINE', 'TICKS', 'LABELS' and 'NAME'. With 'NONE',
            //       complete axes will be suppressed, with 'LINE', only axis lines
            //       will be plotted, with 'TICKS', axis lines and ticks will be
            //       plotted, with 'LABELS' axis lines, ticks and labels will be
            //       plotted and with 'NAME', all axis elements will be displayed.
            //       Default: ('NAME', 'NAME', 'TICKS', 'TICKS').
            //----------------------------------------------------------------------
            const int nx = this->mesh.getNumberOfMeshAlongXAxis();
            const int ny = this->mesh.getNumberOfMeshAlongYAxis();
            this->g.autres(nx, ny);
            this->g.axspos(this->page.panel.pos[0],
                           this->page.panel.pos[1]);
            this->g.ax3len(this->page.panel.len[0],
                           this->page.panel.len[1],
                           this->page.panel.len[2]);
            this->setColorTable(this->page.panel.ctbl);
            this->g.setrgb(0.0, 0.0, 0.0);
            this->setFont(35, "Palatino-Bold");
            this->g.complx();
            this->g.graf3(this->page.panel.xmin[0], this->page.panel.xmax[0],
                          this->page.panel.xmin[0], this->page.panel.xtics[0],
                          this->page.panel.xmin[1], this->page.panel.xmax[1],
                          this->page.panel.xmin[1], this->page.panel.xtics[1],
                          this->page.panel.xmin[2], this->page.panel.xmax[2],
                          this->page.panel.xmin[2], this->page.panel.xtics[2]);
            this->g.crvmat(this->page.panel.data.data(), nx, ny, 1, 1);
            this->setFont(htext_title, "Palatino-Bold");
            this->g.vkytit(-150);
            // [Suppl.] ------------------------------------------------------------
            //   The space between titles and axis systems can be enlarged or
            //   reduced with VKYTIT. By default, the space is 2 * character height.
            //----------------------------------------------------------------------
            this->g.title();
            this->g.endgrf();
            // Plot calculation time
            {
                char text_time[64];
                sprintf(text_time, "T = %6.3f %s",
                        run_param::basic::time * run_param::unit::time / this->page.panel.unit_time,
                        this->page.panel.label_time.c_str());
                int htext = static_cast<int>(this->page.size_panel[0]/20.0);
                int xpos  = this->page.pos_lot[panel_id][0]
                          + this->page.size_mgn_left[panel_id]
                          + static_cast<int>(0.75 * htext);
                int ypos  = this->page.pos_lot[panel_id][1]
                          - this->page.size_mgn_btm[panel_id]
                          - this->page.size_panel[1]
                          + static_cast<int>(0.5 * htext);
                switch (this->page.panel.cplot_qty) {
                    //case ScalarQuantity::Temperature:
                    //    this->g.setrgb(0.0, 0.0, 0.0); // black
                    //    break;
                    default:
                        this->g.setrgb(1.0, 1.0, 1.0); // white
                        break;
                }
                this->setFont(htext, "Palatino-Bold");
                this->g.messag(text_time, xpos, ypos);
            }
            // [Suppl.] ------------------------------------------------------
            //   MESSAG(cstr,nx,ny) plots a text, where 
            //   cstr is a character string (<= 256 characters).
            //   nx,ny are the plot coordinates of the upper left corner.
            //----------------------------------------------------------------
            // Plot panel text
            {
                int htext = static_cast<int>(this->page.size_panel[1]/10.0);
                int xpos  = this->page.pos_lot[panel_id][0]
                          + this->page.size_mgn_left[panel_id]
                          + static_cast<int>(0.25 * htext);
                int ypos  = this->page.pos_lot[panel_id][1]
                          - this->page.size_mgn_btm[panel_id]
                          - this->page.size_panel[1]
                          + static_cast<int>(1.25 * htext);
                switch (this->page.panel.cplot_qty) {
                    //case ScalarQuantity::Temperature:
                    //    this->g.setrgb(0.0, 0.0, 0.0);
                    //    break;
                    default:
                        this->g.setrgb(1.0, 1.0, 1.0);
                        break;
                }
                this->setFont(htext, "Palatino-Bold");
                this->g.messag(panel_text.c_str(), xpos, ypos);
            }
            this->g.setrgb(0.0, 0.0, 0.0); // reset to BLACK

#if 0
            // Color bar
            if (panel_id == this->page.n_panel - 1) {
            //if ((panel_id + 1) % this->page.n_panel_x == 0) {
                int htext = 35;
                this->setFont(htext, "Palatino-Bold");
                this->g.complx();
                int len = this->page.size_page[1]
                        - this->page.size_mgn_top[0]
                        - this->page.size_mgn_btm[panel_id];
                int xpos = this->page.size_page[0]
                         - this->page.size_mgn_right[panel_id]
                         + 0.075 * this->page.size_panel[0];
                int ypos = this->page.size_page[1]
                         - this->page.size_mgn_btm[panel_id];
                this->g.labdig(1, "XYZ");
                const int dir_cb = 0;
#if 1
                this->g.zaxis(this->page.panel.xmin[2],
                              this->page.panel.xmax[2],
                              this->page.panel.xmin[2],
                              this->page.panel.xtics[2],
                              len, this->page.panel.label[2].c_str(),
                              0, dir_cb, xpos, ypos);
#else
                this->g.zaxis(this->page.panel.xmin[2],
                              this->page.panel.xmax[2],
                              this->page.panel.xmin[2],
                              this->page.panel.xtics[2],
                              len, "", 0, dir_cb, xpos, ypos);
                xpos = this->page.pos_lot[panel_id][0]
                     + this->page.size_mgn_left[panel_id]
                     + this->page.size_panel[0]
                     + 0.5 * this->page.size_panel[0];
                ypos = this->page.pos_lot[panel_id][1]
                     - this->page.size_mgn_btm[panel_id]
                     - this->page.size_panel[1];
                     + 0.5 * this->g.nlmess(this->page.panel.label[2].c_str());
                this->g.angle(-90); // rorate the text 90 degrees clock-wise.
                this->g.messag(this->page.panel.label[2].c_str(), xpos, ypos);
                this->g.angle(90); // reset
#endif
            }
#endif
            
            this->g.texmod("OFF");
              
        }
     
        void closePage() {
            this->g.disfin();
        }
     
    };

}
namespace dislin_util = dislin_utilities;
