#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {

    std::vector< std::vector<int> > lol;
    const int n_list = 2;
    lol.resize(n_list);
    lol[0].push_back(0);
    lol[0].push_back(2);
    lol[0].push_back(4);
    lol[1].push_back(1);
    lol[1].push_back(3);
    lol[1].push_back(5);

    //lol.erase(lol.begin() + 1);

    return 0;
}
