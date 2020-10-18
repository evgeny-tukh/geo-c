#define _INTERNAL_

#include <cstdint>
#include <math.h>
#include "geodefs.h"

#ifndef __cplusplus
#include <stdbool.h>
#endif

#ifdef __cplusplus
namespace geo {
#endif

bool calcRhumblinePos (bool useWgs84, Pos *origin, double range, double bearing, Pos *dest);
bool calcRhumblineDistAndBrg (bool useWgs84, Pos *origin, Pos *dest, double *range, double *bearing);

inline bool _calcRhumblinePos (bool useWgs84, double lat, double lon, double range, double bearing, double& destLat, double& destLon) {
    Pos origin { lat, lon }, dest { 0.0, 0.0 };
    
    bool result = calcRhumblinePos (useWgs84, & origin, range, bearing, & dest);

    if (result) {
        destLat = dest.lat;
        destLon = dest.lon;
    }

    return result;
}

inline bool _calcRhumblineDistAndBrg (bool useWgs84, double lat1, double lon1, double lat2, double lon2, double *range, double *bearing) {
    Pos origin { lat1, lon1 }, dest { lat2, lon2 };

    return calcRhumblineDistAndBrg (useWgs84, & origin, & dest, range, bearing);
}

bool calcRhumblineDistAndBrg (bool useWgs84, Pos *origin, Pos *dest, double *range, double *bearing) {
    if (!origin || !dest) return false;

    auto begLat = origin->lat;
    auto begLon = origin->lon;
    auto endLat = dest->lat;
    auto endLon = dest->lon;

    if (range) *range = 0.0;
    if (bearing) *bearing = 0.0;

    if (invalidVal (begLat) || invalidVal (begLon) || invalidVal (endLat) || invalidVal (endLon) ||
        fabs (begLat) > 10000. || fabs (begLon) > 10000. || fabs (endLat) > 10000. || fabs (endLon) > 10000.) return false;

    normalizeLon (& begLon);
    normalizeLon (& endLon);

    // Calculate difference in latitude and longitude. Beware wrap round in longitude
    auto lonDiff = endLon - begLon;

    if (lonDiff <= - PI) lonDiff += TWO_PI;
    else if (lonDiff > PI) lonDiff -= TWO_PI;

    auto latDif = endLat - begLat;
    auto eccentricity = sqrt (WGS84_FLATTENING * 2.0 - WGS84_FLATTENING * WGS84_FLATTENING);

    // Test special case of same meridian
	if (fabs (latDif) < DEFPRECISION) {
        auto partial = eccentricity * sin (endLat);
		
        if (bearing) *bearing = lonDiff > 0.0 ? HALF_PI : HALF_OF_THREE_PI;
		if (range) *range = fabs (WGS84_EQUAT_RAD_M / METERS_IN_NM * lonDiff * cos (endLat) / sqrt (1.0 - partial * partial));

        return true;
	}

    // Special case of parallel sailings
	if (fabs (lonDiff) < DEFPRECISION) {
        auto meridDist = meridionalDist (endLat, eccentricity, WGS84_EQUAT_RAD_M) - meridionalDist(begLat, eccentricity, WGS84_EQUAT_RAD_M);
		
        if (range) *range = fabs (meridDist);
		if (bearing) *bearing = latDif > 0.0 ? 0.0 : 180.0;

        return true;
	}	

    // General case
	auto meridionalDiff = meridionalPart (endLat, eccentricity) - meridionalPart(begLat, eccentricity);
	auto meridDistDiff  = meridionalDist (endLat, eccentricity, WGS84_EQUAT_RAD_M) - meridionalDist (begLat, eccentricity, WGS84_EQUAT_RAD_M);
	
    auto brg = atan2 (DEG_IN_RAD * lonDiff, meridionalDiff);
    auto rng = meridDistDiff / cos (brg);

    normalizeAngle (& brg);

    if (bearing) *bearing = brg;
    if (range) *range = rng;

    return true;
}

bool calcRhumblineDistAndBrg2 (bool useWgs84, Pos *origin, Pos *dest, double *range, double *bearing) {
    if (range) *range = 0.0f;
    if (bearing) *bearing = 0.0f;

    if (!geo::checkPos (origin) || !geo::checkPos (dest)) return false;

    double begLat = origin->lat;
    double begLon = origin->lon;
    double endLat = dest->lat;
    double endLon = dest->lon;

    if (!useWgs84) {
        // Spherical algorythm based on inversed sign for West-East (West positive, East negative)
        begLon = -begLon;
        endLon = -endLon;
    }

    normalizeLon (& begLon);
    normalizeLon (& endLon);

    if (!checkSegmentRag (begLat, begLon, endLat, endLon, true, false)) return false;

    auto latDiff     = endLat - begLat;
    auto westLonDiff = fmod (endLon - begLon, TWO_PI);
    auto eastLonDiff = fmod (begLon - endLon, TWO_PI);
    auto tangentRate = tan (endLat * 0.5 + QUARTER_PI) / tan (begLat * 0.5 + QUARTER_PI);
    
    double dirCosine, deltaPhi;

    if (tangentRate <= 0.0) return false;
    
    if (!useWgs84) {
        if (begLat == endLat) {
            // Rhumb line lays on the equator
            dirCosine = cos (begLat); 
            deltaPhi  = 0.0;
        } else {
            dirCosine = latDiff / (deltaPhi = log (tangentRate));
        }

        // Search the shortest rhumb line
        if (westLonDiff < eastLonDiff) {
            // Westerly rhumb line is the shortest
            *bearing = fmod (atan2 (-westLonDiff, deltaPhi), TWO_PI);
            *range = hypoLen (dirCosine * westLonDiff, latDiff) * MILES_IN_RAD;
        } else {
            // Easterly rhumb line is the shortest
            *bearing  = fmod (atan2 (eastLonDiff, deltaPhi), TWO_PI);
            *range = hypoLen (dirCosine * eastLonDiff, latDiff) * MILES_IN_RAD;
        }
    } else {
        // Use known ellipsoid
        double  deltaEasting,
                deltaNorthing,
                scaleCoef,
                stepDistance = STEP_RANGE,
                middleLat    = (begLat + endLat) * 0.5,
                eccentricity = sqrt (WGS84_FLATTENING + WGS84_FLATTENING - WGS84_FLATTENING * WGS84_FLATTENING),
                tempValue    = eccentricity * sin (middleLat),
                scaledEquRadius;

        scaleCoef       = cos (middleLat) / sqrt (1.0 - tempValue * tempValue);
        scaledEquRadius = scaleCoef * WGS84_EQUAT_RAD_M / METERS_IN_NM;
        deltaPhi        = deltaPhiInline (eccentricity, begLat, endLat);

        deltaEasting  = scaledEquRadius * (westLonDiff < eastLonDiff ? westLonDiff : eastLonDiff);
        deltaNorthing = scaledEquRadius * deltaPhi;
                                
        *range   = hypoLen (deltaEasting, deltaNorthing);
        *bearing = fmod (atan2 (westLonDiff < eastLonDiff ? westLonDiff : -eastLonDiff, deltaPhi), TWO_PI);

        // Distance too big - must be splitted...
        if (*range > stepDistance) {
            double midLat, midLon, milesLeft, midBrg;

            // One step by elementary step distance...
            if (!_calcRhumblinePos (true, begLat, begLon, stepDistance, *bearing, midLat, midLon)) return false;

            // Calc distance for rest of distance (this may be splitted too...)
            if (!_calcRhumblineDistAndBrg (true, midLat, midLon, endLat, endLon, & milesLeft, & midBrg)) return false;

            *range   = stepDistance + milesLeft;
            *bearing = (*bearing + midBrg) * 0.5;
        }
    }

    return true;
}

bool calcRhumblinePos (bool useWgs84, Pos *origin, double range, double bearing, Pos *dest) {
    if (!origin || !dest) return false;

    auto eccentricity = sqrt (WGS84_FLATTENING + WGS84_FLATTENING - WGS84_FLATTENING * WGS84_FLATTENING);
    auto begLat = origin->lat;
    auto begLon = origin->lon;
    auto endLat = begLat;
    auto endLon = begLon;

    dest->lat = dest->lon = 0.0;

    // Belousov 4.4.0.16
    if (invalidVal (begLat) || invalidVal (begLon) || invalidVal (range) || invalidVal (bearing) || fabs (bearing) > 10000. || range > 50000. || range < 0.0) return false;

    normalizeAngle (& bearing);
    normalizeLon (& begLon);

    if (fabs (bearing) < DEFPRECISION || fabs (bearing - PI) < DEFPRECISION) {
        // Special case for course of 0 or 180 degrees
	    endLat = findEndLat (begLat, bearing, range, eccentricity, WGS84_EQUAT_RAD_M);
    } else if (fabs (bearing - HALF_PI) < DEFPRECISION || fabs (bearing - HALF_OF_THREE_PI) < DEFPRECISION) {
    	// Special case for course of 90 or 270
        auto partial  = eccentricity * sin (begLat);
        auto longDist = (range * sqrt (1.0 - partial * partial) / WGS84_EQUAT_RAD_M / METERS_IN_NM * cos (begLat));
        
        endLon += fabs (bearing - HALF_PI) < DEFPRECISION ? longDist : - longDist;
    } else {
        // General case
        endLat = findEndLat (begLat, bearing, range, eccentricity, WGS84_EQUAT_RAD_M);

        auto endLatDeg     = endLat * DEG_IN_RAD;
        auto meridPartDiff = meridionalPart (endLat, eccentricity) - meridionalPart (begLat, eccentricity);
        auto longDist      = fabs (meridPartDiff * tan (bearing)) * RAD_IN_DEG;

        endLon += (bearing < PI)? longDist : - longDist;
    }
	
    // Cope with crossing the date line
    if (endLon > PI)
        endLon -= TWO_PI;
    else if(endLon < -PI)
        endLon += TWO_PI;

    dest->lat = endLat;
    dest->lon = endLon;

    return true;
}

#ifdef __cplusplus
}
#endif
