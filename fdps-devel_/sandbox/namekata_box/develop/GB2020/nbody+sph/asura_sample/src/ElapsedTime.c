#include "config.h"

//#include <time.h>

double GetElapsedTime(void){
	return (MPI_Wtime());
}

//getrusage() 関数を使った時間計測
//プロセスが利用したCPU時間を計測する(高精度)。
//clock の代わりに用いられるが、かなり短い時間になると計測できない。
//その場合は、gettimeofday() を使う。 
#include <sys/resource.h>
#include <sys/time.h>

double GetCPUTime(void){
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec*1e-6);
}

//gettimeofday() 関数を使った時間計測
//経過時間を計測する。CPU時間ではないため、他のプロセスが実行している場合には注意が必要。
//機種依存しないため時間計測ではほとんどの場合 gettimeofday() が使われている。
//例えば、MPICH の MPI_Wtime()関数では gettimeofday() を使っている。 
#include <time.h>

double GetElapsedTime_(void){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + (double)tv.tv_usec*1e-6);
}
