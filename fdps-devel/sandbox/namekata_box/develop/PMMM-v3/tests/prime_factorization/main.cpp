#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include <unordered_map>


std::map<int, int> primeFactorization(const int val) {
    assert(val > 1);
    std::map<int, int> ret;
    int tmp = val;
    int prime_number = 2;
    do {
        // Repeat a division until tmp cannot be divided
        // by the current prime number
        do {
            if (tmp % prime_number == 0) {
               tmp /= prime_number;
               ret[prime_number]++;
            } else break;
        } while (1);
        // Update prime number
        if (tmp > 1) {
            prime_number++;
            if (prime_number == tmp) {
                ret[prime_number]++;
                break;
            }
        } else break;
    } while (1);
    return ret;
}

int getMaximumFactorLessThanOrEqualTo(const int val, const int UL) {
    // Perform prime factorization
    std::map<int, int> mp = primeFactorization(val);
    // Copy information of prime factorization
    std::vector<int> numbers, exponents;
    for (auto itr = mp.begin(); itr != mp.end(); ++itr) {
        numbers.push_back(itr->first);
        exponents.push_back(itr->second);
    }
    // Calculate # of combinations
    int n_combo = 1;
    for (int i = 0; i < exponents.size(); i++)
        n_combo *= (exponents[i] + 1); // +1 comes from the case of 0.
#if 0
    // Check
    std::cout << "n_combo = " << n_combo << std::endl;
    for (int i = 0; i < exponents.size(); i++) {
        std::cout << "i = " << i 
                  << " number = " << numbers[i]
                  << " exponent = " << exponents[i]
                  << std::endl;
    }
#endif
    // Calculate 
    std::vector<int> strides;
    strides.resize(exponents.size());
    strides[0] = 1;
    for (int i = 1; i < exponents.size(); i++) {
        strides[i] = strides[i-1] * (exponents[i-1] + 1);
    }
#if 0
    // Check
    for (int i = 0; i < strides.size(); i++) {
        std::cout << "i = " << i
                  << " stride = " << strides[i]
                  << std::endl;
    }
#endif
    // Find the maximum factor
    int fact_max = -1;
    for (int n = 0; n < n_combo; n++) {
        std::vector<int> powers;
        powers.resize(exponents.size());
        //for (int i = 0; i < powers.size(); i++) powers[i] = 0; // clear
        int tmp = n;
        for (int i = powers.size() - 1; i >= 0; i--) {
            powers[i] = tmp / strides[i];
            tmp -= strides[i] * powers[i];
        }
        assert(tmp == 0);
#if 0
        // Check
        for (int i = 0; i < powers.size(); i++) {
            std::cout << "i = " << i
                      << " power = " << powers[i]
                      << std::endl;
        }
#endif
        // Calculate a factor corresponding to powers
        int fact = 1; 
        for (int i = 0; i < numbers.size(); i++) {
#if 1
            int fact_tmp = 1;
            for (int k = 0; k < powers[i]; k++) fact_tmp *= numbers[i];
            fact *= fact_tmp;
#else
            fact *= std::pow(numbers[i], powers[i]);
#endif
        }
        // Compare
        if (fact == UL) {
            return fact;
        } else if (fact < UL) {
            if (fact > fact_max) fact_max = fact;
        }
    }
    return fact_max;
}

int main(int argc, char *argv[]) {

    // Test for primeFactorization
    const int n_group = 2*3*4*5*6*7;
    std::cout << "org. = " << n_group << std::endl;
    std::map<int, int> mp = primeFactorization(n_group);

    int rslt = 1;
    for (auto itr = mp.begin(); itr != mp.end(); ++itr) {
        std::cout << "key = " << itr->first
                  << ", val = " << itr->second
                  << std::endl;
        rslt *= std::pow(itr->first, itr->second);
    }
    std::cout << "rslt = " << rslt << std::endl;


    // Test for getMaximumFactorLessThanOrEqualTo
    const int val = 2*2*2*3*4*5*6*7*8*9;
    const int UL = 2567;
    std::cout << "val = " << val 
              << " UL = " << UL
              << " rslt = " << getMaximumFactorLessThanOrEqualTo(val, UL)
              << std::endl;

    return 0;
}
