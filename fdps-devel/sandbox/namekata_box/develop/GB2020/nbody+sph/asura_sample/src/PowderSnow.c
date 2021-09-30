#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include "config.h"

/*
 * PowderSnow is a program which can compute 64-bit unique IDs based on time,
 * shard ID, and sequential numbers.  PowderSnow is based on SnowFlake (the
 * twitter ID generator) and the instagram ID generator.  To fit with particle
 * simulations, the time baseline has been changed from 41 years (instagram ID
 * generator) to 4 years. This program can generate 256 unique IDs per
 * milliseconds. 
 */

#define BASE_BIT (64)
#define TIME_BIT (37)   // at least 4 years
#define SHARD_BIT (19)  
#define SEQUENCE_BIT (8)
#define SEQUENCE_BASE (1<<8)

//19:01:18 (JST) 9th May 2018
static unsigned long Offset = 1525860078344;

static bool FirstCall = true;
static int ShardID = 0;

/*
 * This function sets ShardID in accordance with the local parameters, i.e., the
 * process ID and the host ID.  The ShardID has 19 bits. The first 12 bits are
 * filled with the bits come from the host ID and the other 7 bits are that
 * comes from the process ID.
 */
static void PowderSnowSetShardIDbyLocalEnv(void){
    pid_t ProcessID = getpid();
    long int HostID = gethostid();

    HostID &= 0xFFFF;
    long int HostIDLower = (HostID&0xF);
    HostID = (HostID>>4)|HostIDLower;
    ProcessID &= 0x7F;

    ShardID = (HostID<<7)|ProcessID;
    FirstCall = false;
    return ;
}

/*
 * With this function, you can specify the ShardID value.  The first 19 bits are
 * adopted.
 */
void PowderSnowSetShardID(const int id){
    const int ShardBase = (1<<SHARD_BIT);
    ShardID = id%ShardBase;
    FirstCall = false;
    return ;
}

/*
 * You can specify the lower 7 bits of the ShardID value. 
 */
void PowderSnowSetShardIDLowerSevenBits(const int id){
    long int HostID = gethostid();

    HostID &= 0xFFFF;
    long int HostIDLower = (HostID&0xF);
    HostID = (HostID>>4)|HostIDLower;
    int LocalID = id & 0x7F;

    ShardID = (HostID<<7)|LocalID;
    FirstCall = false;

    return ;
}

/*
 * This function returns milliseconds.
 */
static unsigned long GetMilliseconds(void){
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME,&ts);
    return (ts.tv_sec*1000 + ts.tv_nsec/1000000);
}

static unsigned long LastTime = 0;
static unsigned long Sequence = 0;

/*
 * This function returns a unique ID based on the current time, shard ID, and
 * the sequential number.  The return value is unsigned long (64 bits).  This
 * function can generate 256 unique IDs per milliseconds.
 */
unsigned long PowderSnowIDGenerator(void){

    if(FirstCall){
        PowderSnowSetShardIDbyLocalEnv();
    }

    unsigned long CurrentTime = GetMilliseconds();
    while(LastTime > CurrentTime){
        usleep(100);
        CurrentTime = GetMilliseconds();
    }

    if(LastTime!=CurrentTime){
        Sequence = 0;
    }else{
        Sequence = (Sequence+1)%SEQUENCE_BASE;
        if(Sequence == 0){
            CurrentTime = GetMilliseconds();
            while(LastTime >= CurrentTime){
                usleep(1);
                CurrentTime = GetMilliseconds();
            }
        }
    }
    LastTime = CurrentTime;
    return (((CurrentTime-Offset) << (BASE_BIT-TIME_BIT)) | (ShardID<<(BASE_BIT-TIME_BIT-SHARD_BIT)) | Sequence);
}
