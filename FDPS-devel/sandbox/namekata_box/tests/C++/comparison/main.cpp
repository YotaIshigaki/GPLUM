#include <iostream>


int main(int argc, char *argv[]) {
    int x = 1;
    char ch = '1';
    char ch2 = *reinterpret_cast<char *>(&x);
    std::cout << "test 1: " << true << std::endl;
    std::cout << "test 2: " << false << std::endl;
    std::cout << "test 3: " << (x == ch) << std::endl;
    std::cout << "test 4: " << (x == *reinterpret_cast<char *>(&x)) << std::endl;
    std::cout << "test 5: " << (x == ch2) << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "ch = " << ch << std::endl;
    std::cout << "*reinterpret_cast<char *>(&x) = " << ch2 << std::endl;
    std::cout << "(int) *reinterpret_cast<char *>(&x) = " << (int)ch2 << std::endl;
    if (NULL) std::cout << "branch (1)" << std::endl;
    else std::cout << "branch (2)" << std::endl;
    //for (int i = 0; i < 256; i++) 
    //    std::cout << "i = " << i << " (char)i = " << static_cast<char>(i) << std::endl;
    return 0;
}
