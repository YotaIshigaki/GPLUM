#include <iostream>
#include <typeinfo>
#include <type_traits>

template <class T>
class TreeForForce {
public:
    void initialize() {
        static bool first_call = true;
        if (first_call) {
            std::cout << "route (1)" << std::endl;
            first_call = false;
        }
        else std::cout << "route (2)" << std::endl;
    }
};
 
int main()
{
    TreeForForce<int> tree0, tree;
    tree0.initialize();
    tree.initialize();

    return 0;
}
