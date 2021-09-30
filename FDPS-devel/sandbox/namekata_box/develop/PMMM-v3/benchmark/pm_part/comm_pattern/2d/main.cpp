#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cassert>

// DO NOT CHANGE THE VALUES OF THE MACROS BELOW:
#define RINGLIKE                       (0)
#define NESTED_SQUARE                  (1)
#define NESTED_BY_MANHATTAN_DISTANCE   (2)

#define MODEL (2)

enum class SequencerType {Ringlike, NestedSquare, NestedByManhattanDistance};

class Sequencer {
private:
    bool is_initialized_;
    int n_proc_x_, n_proc_y_, n_proc_;
    SequencerType type_;
    
    void calcPosition(const int rank, int & x, int & y) {
        assert(is_initialized_);
        y = rank/n_proc_x_;
        x = rank - n_proc_x_ * y;
    }

public:

    Sequencer(): is_initialized_(false) {}

    void initialize(const int n_proc_x,
                    const int n_proc_y,
                    const SequencerType type) {
        assert(n_proc_x > 0 && n_proc_y > 0);
        n_proc_x_ = n_proc_x;
        n_proc_y_ = n_proc_y;
        n_proc_   = n_proc_x_ * n_proc_y_;
        type_ = type;
        is_initialized_ = true;
    }


    void calcRelativePositionOfCommPartner(const int kth, const int my_rank,
                                           int & x, int & y, int & vx, int & vy,
                                           const bool debug = false) {
        // Calculate the position of my_rank
        calcPosition(my_rank, x, y);

        if (type_ == SequencerType::Ringlike) {
            // Calculate the destination rank
            const int dest = (my_rank + kth) % n_proc_;
            int x_dest, y_dest;
            calcPosition(dest, x_dest, y_dest);
            // Calculate the relative position of destination rank (minimum image)
            int vx = x_dest - x;
            int vy = y_dest - y;
            if (vx >  n_proc_x_/2) vx -= n_proc_x_;
            if (vx < -n_proc_x_/2) vx += n_proc_x_;
            if (vy >  n_proc_y_/2) vy -= n_proc_y_;
            if (vy < -n_proc_y_/2) vy += n_proc_y_;

        } else if (type_ == SequencerType::NestedSquare) {
            // Calculate the current level
            int lev {0};
            int k_start {0}, k_end {0};
            for(;;) {
                if (k_start <= kth && kth <= k_end) break;
                lev++;
                k_start = k_end + 1;
                k_end = (2*lev+1) * (2*lev+1) - 1;
            }
            //std::cout << "lev = " << lev << " k_start = " << k_start << " k_end = " << k_end << std::endl;
            // Set the default start point (correspond to k_start)
            vx = lev;
            vy = 0;
            // Set the initial shift of the starting point
            int n_shift_offset {0};
#if 1
#if 0
            // Pattern #1
            n_shift_offset = (my_rank % 4) * lev;
#else
            // Pattern #2
            {
                const int rem_x = x % 2;
                const int rem_y = y % 2;
                if (rem_x == 1 && rem_y == 0) n_shift_offset = 2*lev;
                else if (rem_x == 1 && rem_y == 1) n_shift_offset = 2*(2*lev);
                else if (rem_x == 0 && rem_y == 0) n_shift_offset = 3*(2*lev);
            }
#endif
#endif
            // Calculate the relative position of destination
            int n_shift = kth - k_start + n_shift_offset;
            for(;;) {
                if (vx == lev && vy < lev) {
                    while (n_shift) {
                        if (vy < lev) {
                            vy++; n_shift--;
                        }
                        if (vy == lev) break;
                    }
                } else if (-lev < vx && vy == lev) {
                    while (n_shift) {
                        if (vx > -lev) {
                            vx--; n_shift--;
                        }
                        if (vx == -lev) break;
                    }
                } else if (vx == -lev && -lev < vy) {
                    while (n_shift) {
                        if (vy > -lev) {
                            vy--; n_shift--;
                        }
                        if (vy == -lev) break;
                    }
                } else if (vx < lev && vy == -lev) {
                    while (n_shift) {
                        if (vx < lev) {
                            vx++; n_shift--;
                        }
                        if (vx == lev) break;
                    }
                }
                if (n_shift == 0) break;
            }
        } else if (type_ == SequencerType::NestedByManhattanDistance) {
            // Calculate the current distance
            // Note that the number of cells at Manhattan distance d is
            // given by 4*d for d >= 1 in two dimension. For d = 0,
            // the number of cells is 1. 
            int dis {0};
            int k_start {0}, k_end {0};
            for (;;) {
                if (k_start <= kth && kth <= k_end) break;
                dis++;
                k_start = k_end + 1;
                k_end = k_start + 4*dis - 1;
            }
            if (debug) std::cout << "dis = " << dis << " k_start = " << k_start << " k_end = " << k_end << std::endl;
            // Set the default starting posint (correspond to k_start)
            vx = dis;
            vy = 0;
            if (debug) std::cout << "vx = " << vx << " vy = " << vy << std::endl;
            // Set the initial shift of the starting point
            int n_shift_offset {0};
#if 1
#if 0
            // Pattern #1
            n_shift_offset = (my_rank % 4) * dis;
#else
            // Pattern #2
            {
                const int rem_x = x % 2;
                const int rem_y = y % 2;
                if (rem_x == 1 && rem_y == 0) n_shift_offset = dis;
                else if (rem_x == 1 && rem_y == 1) n_shift_offset = 2*dis;
                else if (rem_x == 0 && rem_y == 1) n_shift_offset = 3*dis;
            }
#endif
#endif
            // Calculate the relative position of destination rank
            int n_shift = kth - k_start + n_shift_offset;
            if (debug)
                std::cout << "n_shift = " << n_shift
                          << " n_shift_offset = " << n_shift_offset
                          << std::endl;
            for (;;) {
                if (debug) std::cout << "vx = " << vx << " vy = " << vy << std::endl;
                // Determine the direction of shift
                if ((1 <= vx && vx <= dis) && vy >= 0) {
                    while (n_shift) {
                        if (vy < dis) {
                            vx--; vy++; n_shift--;
                        }
                        if (vy == dis) break;
                    }
                } else if (vx <= 0 && (1 <= vy && vy <= dis)) {
                    while (n_shift) {
                        if (vx > -dis) {
                            vx--; vy--; n_shift--;
                        }
                        if (vx == -dis) break;
                    }
                } else if ((-dis <= vx && vx <= -1) && vy <= 0) {
                    while (n_shift) {
                        if (vy > -dis) {
                            vx++; vy--; n_shift--;
                        }
                        if (vy == -dis) break;
                    }
                } else if (vx >= 0 && (-dis <= vy && vy <= -1)) {
                    while (n_shift) {
                        if (vx < dis) {
                            vx++; vy++; n_shift--;
                        }
                        if (vx == dis) break;
                    }
                }
                if (n_shift == 0) break;
            }
        }
    }

};


int main(int argc, char *argv[]) {

    // Node geometry
    constexpr int n_proc_x = 8;
    constexpr int n_proc_y = 8;
    constexpr int n_proc = n_proc_x * n_proc_y;

    // List of destination ranks
    bool is_registered[n_proc][n_proc];
    for (int i=0; i<n_proc; i++)
        for (int j=0; j<n_proc; j++)
            is_registered[i][j]=false;
    std::vector<int> x_list[n_proc], y_list[n_proc], vx_list[n_proc], vy_list[n_proc];

    // Compute the communication pattern
    bool debug {false};
    for (int k = 0;; k++) { // corerspond to kth MPI_Send
        for (int my_rank = 0; my_rank < n_proc; my_rank++) { // my_rank is a sender
            // Output information
            debug = false;
            //if (my_rank == 0 && (k == 25 || k == 33)) debug = true;
            if (debug) std::cout << "k = " << k << " my_rank = " << my_rank << std::endl;
            // Compute x,y,vx,vy
            int x, y, vx, vy;
            Sequencer seq;
#if MODEL == RINGLIKE
            seq.initialize(n_proc_x, n_proc_y, SequencerType::Ringlike);
#elif MODEL == NESTED_SQUARE
            seq.initialize(n_proc_x, n_proc_y, SequencerType::NestedSquare);
#elif MODEL == NESTED_BY_MANHATTAN_DISTANCE
            seq.initialize(n_proc_x, n_proc_y, SequencerType::NestedByManhattanDistance);
#else
#error The macro `MODEL` has an invalid value.
#endif
            seq.calcRelativePositionOfCommPartner(k, my_rank, x, y, vx, vy, debug);
            // Register if the rank indicated by (x,y,vx,vy) is first appearance
            int x_dest = x + vx;
            int y_dest = y + vy;
            while (x_dest < 0) { x_dest += n_proc_x; }
            while (x_dest >= n_proc_x) { x_dest -= n_proc_x; }
            while (y_dest < 0) { y_dest += n_proc_y; }
            while (y_dest >= n_proc_y) { y_dest -= n_proc_y; }
            int dest = x_dest + y_dest * n_proc_x;
            if (dest < 0 || n_proc <= dest) {
                std::cout << "k = " << k
                          << " dest = "  << dest
                          << " x_dest = " << x_dest
                          << " y_dest = " << y_dest
                          << " x = " << x
                          << " y = " << y
                          << " vx = " << vx
                          << " vy = " << vy
                          << std::endl;
                assert(0 <= dest && dest < n_proc);
            }
            if (is_registered[my_rank][dest] == false) {
                x_list[my_rank].push_back(x);
                y_list[my_rank].push_back(y);
                vx_list[my_rank].push_back(vx);
                vy_list[my_rank].push_back(vy);
                is_registered[my_rank][dest] = true;
            }
            if (debug)
                std::cout << "k = " << k
                          << " dest = "  << dest
                          << " x_dest = " << x_dest
                          << " y_dest = " << y_dest
                          << " x = " << x
                          << " y = " << y
                          << " vx = " << vx
                          << " vy = " << vy
                          << std::endl;
        }
        // Espace from the infinite loop if possible
        bool all_clear {true};
        for (int i=0; i<n_proc; i++)
            for (int j=0; j<n_proc; j++)
                if (is_registered[i][j] == false) 
                    all_clear = false;
        if (all_clear) break;
    }

    // Output to a file
    const std::string filename = "comm_pattern_2d.txt";
    std::ofstream output_file;
    output_file.open(filename.c_str(),std::ios::trunc);
    for (int k=0; k < n_proc; k++) {
        for (int my_rank=0; my_rank < n_proc; my_rank++) {
            output_file << x_list[my_rank][k] << "    "
                        << y_list[my_rank][k] << "    "
                        << vx_list[my_rank][k] << "    "
                        << vy_list[my_rank][k]
                        << std::endl;
        }
        output_file << std::endl << std::endl; // separator 
    } 
    output_file.close();

    return 0;
}
