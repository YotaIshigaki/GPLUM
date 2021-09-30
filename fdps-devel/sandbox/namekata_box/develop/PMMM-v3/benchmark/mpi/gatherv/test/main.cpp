// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>

int main(int argc, char *argv[]) {

    const int n_proc = 1024;
    const long seed=19810614;
    srand48(seed);

#if 1
    const int data_size_max=1024*1024;
    for (int data_size=8; data_size<=data_size_max; data_size *= 4) {
        const int recvcount_tot = data_size * n_proc;
        // Set fractions using a pseudo-random number generator
        std::vector<double> fractions;
        double fsum {0.0};
        for (int k=0; k<n_proc; k++) {
            double f = drand48();
            fractions.push_back(f);
            fsum += f;
            
        }
        for (int k=0; k<n_proc; k++) fractions[k] /= fsum; // scale
        // Set recvcounts
        std::vector<int> recvcounts;
        int csum {0};
        for (int k=0; k<n_proc; k++) {
            int cnt = recvcount_tot * fractions[k];
            recvcounts.push_back(cnt);
            csum += cnt;
        }
        if (csum != recvcount_tot) {
            int diff = recvcount_tot - csum;
            for (int k=0; k<diff; k++) recvcounts[k] += 1;
        }
        // Set recvdispls
        std::vector<int> recvdispls;
        recvdispls.resize(n_proc);
        recvdispls[0]=0;
        for (int k=1; k<n_proc-1; k++) 
            recvdispls[k] = recvdispls[k-1] + recvcounts[k-1];
        // Check the setting of recvcounts
        int sum = 0;
        for (int k=0; k<n_proc; k++) sum += recvcounts[k];
        if (sum != recvcount_tot) {
            std::cout << "recvcounts is wrong:" << std::endl;
            std::cout << "    data_size     = " << data_size << std::endl;
            std::cout << "    sum           = " << sum << std::endl;
            std::cout << "    recvcount_tot = " << recvcount_tot << std::endl;
        }
    }
#else
    std::vector<int> data;
    data.push_back(5);
    data.push_back(1);
    data.push_back(2);
    data.push_back(3);
    data.push_back(4);

    std::sort(data.begin(),data.end());
    {for(int i=0;i<data.size();i++) printf("%d\n",data[i]);}
#endif

    return 0;
}
