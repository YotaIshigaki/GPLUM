#include <iostream>
#include <typeinfo>
#include <type_traits>

int main()
{
    const int n_proc = 16;
    const int my_rank = 2;
    for (int jj=0; jj<n_proc; jj++) {
        const int rank_recv = (my_rank + jj) % n_proc;
        const int rank_send = (n_proc + my_rank - jj) % n_proc;
        std::cout << "rank_recv = " << rank_recv 
                  << " rank_send = " << rank_send
                  << std::endl;
    }

    return 0;
}
