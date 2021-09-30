/*
 * This file includes allocation functions and memory release functions.
 */
#ifndef	__ALLOCATION_H_INCLUDED__
#define	__ALLOCATION_H_INCLUDED__
#include "config.h"

void *Allocate(size_t nSize){

    void *p;

    p = malloc(nSize);

    if(p == NULL){
        fprintf(stderr,"Memory Allocation Error:%s:line%d:%s()\n",
                __FILE__,__LINE__,__FUNCTION__);
        exit(AllocationError);
    }

    return (p);
}

void *ReAllocate(void *ptr, size_t nSize){

    void *p;

    p = realloc(ptr,nSize);

    if(p == NULL){
        fprintf(stderr,"Memory Allocation Error:%s:line%d:%s()\n",
                __FILE__,__LINE__,__FUNCTION__);
        exit(AllocationError);
    }

    return (p);
}

void Delete(void *p){
    free(p);
    return;
}
#endif
