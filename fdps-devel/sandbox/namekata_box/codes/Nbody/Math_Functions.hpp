/* Standard headers */
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cfloat>
/* User-defined headers */
#include "preprocess_keywords.h"

namespace math_functions {

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
///////////////////////////   < S I G N >   ///////////////////////////
/*-------------------------------------------------------------------*/
template <typename T>
inline T sign(T val){
    return static_cast<T> (1 - (val <= 0) - (val < 0));
}


} // END of math_functions
