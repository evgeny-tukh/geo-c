#pragma once

#include <float.h>
#include <math.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace geo {
#else
struct geo {
#endif

static const double RAD_IN_MILE = 2.90888208665721596153703703703e-4;
static const double RAD_IN_DEG = 0.01745329251994329576923690768489;
static const double DEG_IN_RAD = 57.295779513082320876798154814105;
static const double MILES_IN_RAD = 3437.7467707849392526107818606515;
static const double DEG_IN_MILE = 1.0 / 60.0;
static const double MILES_IN_DEG = 60.0;
static const double METERS_IN_NM = 1852.0;
static const double NM_IN_METER = 1.0 / METERS_IN_NM;

static const double WGS84_EQUAT_RAD_M = 6378137.0;
static const double WGS84_POLAR_RAD_M = 6356752.3142451793;
static const double WGS84_FLATTENING = 3.35281066474751169502944198282e-3;
static const double SPHERE_RAD_M = 6366707.0194937074958298109629434;

static const double PI = 3.1415926535897932384626433832795;
static const double TWO_PI = PI * 2.0;
static const double HALF_PI = PI * 0.5;
static const double HALF_OF_THREE_PI = PI * 1.5;
static const double QUARTER_PI = PI * 0.25;
static const double ONE_SIXTH_PI = PI / 6.0;

static const double ONE_SIXTH = 1.0 / 6.0;
static const double ONE_SIXTEENTH = 1.0 / 16.0;

static const double STEP_RANGE = 10.0;

static const double DEFPRECISION = 1.0E-10;

enum Operation {
    CALC_BRG_RNG = 'b',
    CALC_DEST_POINT = 'p',
};

enum Method {
    RHUMBLINE = 'r',
    GREAT_CIRCLE = 'g',
};

struct Pos {
    double lat, lon;
};

#ifdef __cplusplus
}
#else
};
#endif

#ifdef _INTERNAL_
namespace geo {
    inline bool checkLat (double lat, bool assumeRad = false) {
        auto range = assumeRad ? 1.3351768777756621263466234378938 : 85.0;

        return lat >= -range && lat <= range;
    }

    inline bool checkLon (double lat, bool assumeRad = false, bool strict = false) {
        auto range = strict ? (assumeRad ? PI : 180.0) : 1000.0;

        return lat >= -range && lat <= range;
    }

    inline bool checkPos (Pos *pos, bool strict = false) {
        return checkLat (pos->lat, strict) && checkLon (pos->lon, strict);
    }

    inline void normalizeLon (double *lon, bool assumeRad = true) {
        double range = assumeRad ? PI : 180.0;

        while (*lon < -range) *lon += (range + range);
        while (*lon > range) *lon -= (range + range);
    }

    inline void normalizeAngle (double *angle, bool assumeRad = true) {
        double range = assumeRad ? TWO_PI : 360.0;

        while (*angle < 0.0) *angle += range;
        while (*angle > range) *angle -= range;
    }

    inline bool isPole (Pos *point, bool assumeRadians = false)
    {
        auto range = assumeRadians ? 0.5 : 90.0;

        return point->lat == range || point->lat == -range;
    }

    inline bool isPole (double lat, bool assumeRadians = false)
    {
        auto range = assumeRadians ? 0.5 : 90.0;

        return lat == range || lat == -range;
    }

    inline bool checkGeoPointRange (double lat, double lon, bool assumeRad = true, bool enablePole = false) {
        return !(!enablePole && isPole (lat, assumeRad) || !checkLat (lat, assumeRad) || !checkLon (lon, assumeRad));
    }

    inline bool checkSegmentRag (double begLat, double begLon, double endLat, double endLon, bool assumeRad, bool enablePole) {
        return checkGeoPointRange (begLat, begLon, assumeRad, enablePole) && checkGeoPointRange (endLat, endLon, assumeRad, enablePole);
    }

    inline double hypoLen (double cat1, double cat2) {
        return sqrt (cat1 * cat1 + cat2 * cat2);
    }

    inline double deltaPhiInline (double eccentricity, double begLat, double endLat)
    {
        double eSinBegLat = eccentricity * sin (begLat),
            eSinEndLat = eccentricity * sin (endLat);

        return log (
            tan (QUARTER_PI + endLat * 0.5) /
            tan (QUARTER_PI + begLat * 0.5) * 
            pow (((1.0 - eSinEndLat) * (1.0 + eSinBegLat)) / ((1.0 + eSinEndLat) * (1.0 - eSinBegLat)), eccentricity * 0.5)
        );
    }

    inline bool invalidVal (double val) {
        return (_fpclass ((double) (val)) & (_FPCLASS_SNAN | _FPCLASS_QNAN | _FPCLASS_NINF | _FPCLASS_PINF));
    }

    inline double meridionalDist (double lat, double eccentricity, double equRadius)
    {
        auto e2 = eccentricity * eccentricity;
        auto e4 = e2 * e2;
        auto e6 = e4 * e2;
        auto v1 = RAD_IN_DEG * (1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 0.01953125 /*5.0 / 256.0*/ * e6);
        auto v2 = 0.375 /*3.0 / 8.0*/ * (e2 - e4 / 4.0 + 0.1171875 /*15.0 / 128.0*/ * e6);
        auto v3 = 0.05859375 /*15.0 / 256.0*/ * (e4 + 0.75 /*3.0 / 4.0*/ * e6);
        auto v4 = equRadius / METERS_IN_NM;
        auto md = v4 * (v1 * lat * DEG_IN_RAD - v2 * sin (2.0 * lat) + v3 * sin (4.0 * lat));

        return md;
    }

    inline double findEndLat (double begLat, double bearing, double range, double eccentricity, double equatorialRadius) {
        auto meridRngFrom  = meridionalDist (begLat, eccentricity, equatorialRadius);
        auto meridRngEstim = meridRngFrom + range * cos (bearing);
        auto lat           = begLat + (meridRngEstim - meridRngFrom) / 60.0;

        for (auto iteration = 0; iteration < 10; ++ iteration)
        {
            auto meridRngPart = meridionalDist (lat * RAD_IN_DEG, eccentricity, equatorialRadius);

            lat += (meridRngEstim - meridRngPart) / 60.0;
        }

        return lat * RAD_IN_DEG;
    }

    inline double meridionalPart (double lat, double eccentricity) {
        auto part = eccentricity * sin (lat);
        
        return DEG_IN_RAD * (log (tan (QUARTER_PI + lat * 0.5)) - 0.5 * eccentricity * log ((1.0 + part) / (1.0 - part)));
    }

    inline double calcLatEquationRoot (
        double eccentricity,                  // Geoid eccentricity
        double begLat,                        // Begin latitude
        double dNorthing,                     // Mercator northing offset miles
        double scaledEquRadius,               // a * k0
        double precision
    ) {
        auto step = 0.001;
        auto value = dNorthing / scaledEquRadius;

        double begPhi, endPhi, midPhi, begValue, midValue, endValue;

        // Search for interval range
        begPhi   = begLat - step;
        begValue = deltaPhiInline (eccentricity, begLat, begPhi);
        endPhi   = begLat + step;
        endValue = deltaPhiInline (eccentricity, begLat, endPhi);

        if (endValue >= value && begValue <= value)
        {
            // Range is found - no some actions needed
        }
        else if (endValue <= value && begValue >= value)
        {
            // Range is found but must be exachanged
            double gvTemp = begPhi;

            begPhi = endPhi;
            endPhi = gvTemp;
        }
        else if (endValue > begValue && endValue > value)
        {
            // Shift begPhi down until this is great than value
            do
            {
                endPhi   = begPhi;
                begPhi  -= step; 
                endValue = begValue;
                begValue = deltaPhiInline (eccentricity, begLat, begPhi);
                step    += step;
            }
            while (begValue > value);
        }
        else
        {
            // Shift endPhi upper until this is less than value
            do
            {
                begPhi   = endPhi;
                endPhi  += step; 
                begValue = endValue;
                endValue = deltaPhiInline (eccentricity, begLat, endPhi);
                step    += step;
            }
            while (endValue < value);
        }

        // Range are found - start linear interpolating...
        do
        {
            midPhi   = begPhi + 
                        (endPhi - begPhi) * (value - begValue) / (endValue - begValue);
            midValue = deltaPhiInline (eccentricity, begLat, midPhi);

            // Change range...
            if (midValue >= value)
            {
                endPhi   = midPhi;
                endValue = midValue;
            }
            else
            {
                begPhi   = midPhi;
                begValue = midValue;
            }
        }
        while (fabs (midValue - value) > precision);

        return midPhi;
    }

    inline bool isSame (double val1, double val2, double precision = 1.0e-10) {
        return fabs (val1 - val2) <= precision;
    }

    inline bool isNotSame (double val1, double val2, double precision = 1.0e-10) {
        return fabs (val1 - val2) > precision;
    }

    inline void checkBearing (double *brg) {
        while (*brg < 0.0) *brg += TWO_PI;
        while (*brg > TWO_PI) *brg -= TWO_PI;
    }

    inline void reverseBearing (double *brg) {                                       \
        if (*brg > PI) 
            *brg -= PI;
        else
            *brg += PI;

        normalizeAngle (brg);
    }

    inline bool isPole (Pos *point) {
        return isPole (point->lat, true);
    }

    inline bool wrongLat (double lat) {
        return lat < - HALF_PI || lat > HALF_PI;
    }

    inline bool wrongLong (double lon) {
        return lon < -PI || lon > PI;
    }

    inline bool wrongAngle (double val) {
        return val < 0.0 || val > TWO_PI;
    }

    inline bool wrongDistance (double dist) {
        return dist < 0.0 || dist > WGS84_EQUAT_RAD_M * TWO_PI * NM_IN_METER;
    }

    inline bool checkBegPointCourseDist (double begLat, double begLon, double bearing, double range, bool enablePole) {
        if (!checkGeoPointRange (begLat, begLon, true, enablePole)) return false;

        return !wrongAngle (bearing) && !wrongDistance (range);
    }
}
#endif