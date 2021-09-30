#include <iostream>

typedef struct {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
} MAGI_Tipsy_Header;

int main(int argc, char *argv[]) {

    MAGI_Tipsy_Header header;
    
    printf("sizeof(MAGI_Tipsy_Header) = %d\n",sizeof(header));

    return 0;
}
