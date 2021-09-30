#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
// Header file of FDPS
#include <particle_simulator.hpp>

template <class T>
void MergeOrthotopeList(std::vector<T> & from_list,
                        std::vector<T> & to_list) {
}

template <>
void MergeOrthotopeList<PS::F64ort>(std::vector<PS::F64ort> & from_list,
                                    std::vector<PS::F64ort> & to_list) {
    while (!from_list.empty()) {
        // Extract vertex from `from_list`
        PS::F64ort f = from_list.back();
        from_list.pop_back();
        // Find elements of dest which intersects with vertex above.
        std::vector<PS::S32> delete_list;
        PS::F64ort rslt = f;
        for (PS::S32 i = 0; i < to_list.size(); i++) {
            if (f.overlapped(to_list[i])) {
                delete_list.push_back(i);
                rslt.merge(to_list[i]);
            }
        }
        // Delete
        for (const auto & idx: delete_list) {
            to_list.erase(to_list.begin() + idx);
        }
        // Add
        to_list.push_back(rslt);
    };

}

int main(int argc, char *argv[]) {

    // Make lists
#if 0
    std::vector<PS::F64ort> from_list;
    from_list.push_back(PS::F64ort(PS::F64vec(0.0, 0.0, 0.0),
                                   PS::F64vec(1.0, 1.0, 1.0)));
    from_list.push_back(PS::F64ort(PS::F64vec(1.0, 1.0, 1.0),
                                   PS::F64vec(2.0, 2.0, 2.0)));

    std::vector<PS::F64ort> to_list;
    to_list.push_back(PS::F64ort(PS::F64vec(0.5, 0.5, 0.5),
                                 PS::F64vec(1.5, 1.5, 1.5)));
    to_list.push_back(PS::F64ort(PS::F64vec(2.1, 2.1, 2.1),
                                 PS::F64vec(2.5, 2.5, 2.5)));
#else
    std::vector<PS::F64ort> from_list;
    from_list.push_back(PS::F64ort(PS::F64vec(0.0, 0.0, 0.0),
                                   PS::F64vec(1.0, 1.0, 1.0)));
    from_list.push_back(PS::F64ort(PS::F64vec(1.0, 1.0, 1.0),
                                   PS::F64vec(2.0, 2.0, 2.0)));

    std::vector<PS::F64ort> to_list;
    to_list.push_back(PS::F64ort(PS::F64vec(0.0, 0.0, 0.0),
                                 PS::F64vec(1.0, 1.0, 1.0)));
    to_list.push_back(PS::F64ort(PS::F64vec(1.0, 1.0, 1.0),
                                 PS::F64vec(2.0, 2.0, 2.0)));
#endif

    // Merge
    MergeOrthotopeList(from_list, to_list);


    // Check the result
    for (const auto & e : to_list) {
        std::cout << e << std::endl;
    }



    return 0;
}
