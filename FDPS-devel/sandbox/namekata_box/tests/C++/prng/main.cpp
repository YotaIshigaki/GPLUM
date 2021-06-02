/* C++ headers */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
/* user-defined headers */
#include "prng.h"

int main(int argc, char *argv[]) {
    // Test of RestorablePRNG
    {
        std::size_t seed = 0;
        std::size_t count = 0;
        RestorablePRNG<std::mt19937_64> gen;
        gen.init(seed, count);
        std::uniform_real_distribution<double> dist_r {0, 1.0};
        double r = dist_r(gen);
        gen.printCount();
        std::uniform_int_distribution<long> dist_l(0,10);
        long l = dist_l(gen);
        gen.printCount();
    }

    // Test restoring from binary file
    {
        std::mt19937_64 gen_a;
        std::uniform_int_distribution<int> dist(0,256);
        const int n = 8;
        for (int i=0; i<n; i++) std::cout << dist(gen_a) << std::endl;

        std::ofstream ofs;
        const std::string filename = "gen.dat";
        ofs.open(filename.c_str(),std::ios::trunc | std::ios::binary);
        ofs.write((char *)&gen_a, sizeof(gen_a));
        ofs.close();


        std::mt19937_64 gen_b;
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
        ifs.read((char *)&gen_b, sizeof(gen_b));
        ifs.close();

        for (int i=0; i<n; i++) {
            std::cout << "[a] i = " << i << dist(gen_a) << std::endl;
            std::cout << "[b] i = " << i << dist(gen_b) << std::endl;
        }
    }

    // Test restoring from ascii file
    {
        std::mt19937_64 gen_a;
        std::uniform_int_distribution<int> dist(0,256);
        const int n = 8;
        for (int i=0; i<n; i++) std::cout << dist(gen_a) << std::endl;

        std::ofstream ofs;
        const std::string filename = "gen.dat";
        ofs.open(filename.c_str(),std::ios::trunc);
        ofs << gen_a;
        ofs.close();


        std::mt19937_64 gen_b;
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in);
        ifs >> gen_b;
        ifs.close();

        for (int i=0; i<n; i++) {
            std::cout << "[a] i = " << i << dist(gen_a) << std::endl;
            std::cout << "[b] i = " << i << dist(gen_b) << std::endl;
        }
    }

    return 0;
}
