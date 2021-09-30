#include <iostream>

int main(int argc, char *argv[]) {

    const int N = 8;
    int *buf = new int[N*N*N];
    int ***cell = reinterpret_cast<int (*)[N][N]>(buf);

    return 0;
}
