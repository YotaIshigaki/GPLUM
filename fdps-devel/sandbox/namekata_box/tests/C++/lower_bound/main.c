#include <iostream>
#include <vector>
#include <algorithm>

int main(int argc, char *argv[]) {

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    // Find 0
    {
        double val = 0;
        decltype(v)::iterator it = std::lower_bound(v.begin(), v.end(), val);
        if (it != v.end()) {
            std::size_t pos = std::distance(v.begin(), it);
            std::cout << "val = " << val
                      << ": " << *it << " pos=" << pos << std::endl;
        } else {
            std::cout << val << " not found." << std::endl;
        }
    }

    // Find 3.5
    {
        double val = 3.5;
        decltype(v)::iterator it = std::lower_bound(v.begin(), v.end(), val);
        if (it != v.end()) {
            std::size_t pos = std::distance(v.begin(), it);
            std::cout << "val = " << val 
                      << ": "  << *it << " pos=" << pos << std::endl;
        } else {
            std::cout << val << " not found." << std::endl;
        }
    }

    // Find 7
    {
        double val = 7.0;
        decltype(v)::iterator it = std::lower_bound(v.begin(), v.end(), val);
        if (it != v.end()) {
            std::size_t pos = std::distance(v.begin(), it);
            std::cout << "val = " << val 
                      << ": "  << *it << " pos=" << pos << std::endl;
        } else {
            std::cout << val << " not found." << std::endl;
        }
    }

    return 0;
}
