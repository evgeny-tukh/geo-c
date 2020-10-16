#pragma once

#include "geodefs.h"

#ifdef __cplusplus
namespace geo {
#endif

bool calcRhumblinePos (bool useWgs84, Pos *origin, double range, double bearing, Pos *dest);
bool calcRhumblineDistAndBrg (bool useWgs84, Pos *origin, Pos *dest, double *range, double *bearing);

inline double valToRad (double val) { return val * RAD_IN_DEG; }
inline double valToDeg (double val) { return val / RAD_IN_DEG; }
inline void degToRad (double *val) { *val *= RAD_IN_DEG; }
inline void radToDeg (double *val) { *val /= RAD_IN_DEG; }

inline void ptToRad (Pos *pos) {
    pos->lat *= RAD_IN_DEG; 
    pos->lon *= RAD_IN_DEG; 
}
inline void ptToDeg (Pos *pos) {
    pos->lat /= RAD_IN_DEG;
    pos->lon /= RAD_IN_DEG;
}

#ifdef __cplusplus
}
#endif
