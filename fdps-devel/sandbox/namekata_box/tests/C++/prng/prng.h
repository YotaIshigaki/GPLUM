#pragma once
/* C++ headers */
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

template <class Tengine>
class RestorablePRNG {
public:
    using result_type = typename Tengine::result_type;
    /* member variables */
private:
    std::size_t seed_;
    std::size_t count_;
    Tengine engine_;

    /* member functions */
public:
    void init(std::size_t seed, std::size_t cnt) {
        this->seed_ = seed;
        this->count_ = cnt;
        this->engine_.seed(this->seed_);
    }

    void printCount() const {
        std::cout << "count_ = " << this->count_ << std::endl;
    }

    static result_type min() { return Tengine::min(); }
    static result_type max() { return Tengine::max(); }
    result_type operator () () {
        const result_type res = engine_();
        ++count_;
        return res;
    }
};
