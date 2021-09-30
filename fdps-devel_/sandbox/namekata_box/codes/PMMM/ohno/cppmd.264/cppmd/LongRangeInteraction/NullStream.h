#ifndef NULLSTREAM_H
#define NULLSTREAM_H

#include <iostream>

typedef struct nullstream: std::ostream {
    nullstream(): std::ios(0), std::ostream(0) {}
} NullStream;

#endif
