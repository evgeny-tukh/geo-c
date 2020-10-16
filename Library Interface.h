/******************************************************************************
 *                                                                            *
 * MARIS spherical and flat geometry utilities                                *
 * Written by Eugeny Toukh                                                    *
 * Started at 18.05.2000                                                      *
 *                                                                            *
 * Flat functions interface                                                   *
 * Allows to call all functions indirect, without some objects                *
 *                                                                            *
 ******************************************************************************/

#ifndef _MARIS_GEOMETRY_INTERFACE_
#define _MARIS_GEOMETRY_INTERFACE_

#pragma pack(1)

// **************************************************************************************************
// Copy data from Global Definitions.h
//
#ifndef _MARIS_GEOMETRY_GLOBAL_DEFS_
#define _MARIS_GEOMETRY_GLOBAL_DEFS_

#include <StdLib.h>

#define ARR_SIZE(array)         (sizeof (array) / sizeof (*(array)))

// Undefined cooredinate
#define GD_UNDEFINED            1.0E99

// Default values
#define GD_DEFPRECISION         1.0E-10
#define GD_DEFSTEPDISTANCE      10.0

// Zero range
#define GD_IZERO                0
#define GD_ZERO                 0.0
#define GD_ZERORANGE            1.0E-14

// Main true rates
#define GD_DIV_2                0.5
#define GD_DIV_3                0.33333333333333333333333333333333
#define GD_DIV_4                0.25
#define GD_DIV_5                0.2
#define GD_DIV_6                0.16666666666666666666666666666667
#define GD_DIV_10               0.1
#define GD_DIV_16               0.0625

// Natural logarithm base
#define GD_E                    2.7182818284590452353602874713527

// Different values related PI value
#define GD_PI_RAD               3.1415926535897932384626433832795   // PI in radians
#define GD_PI_DEG               180.0                               // PI in degree
#define GD_2_PI_RAD             6.283185307179586476925286766559
#define GD_2_PI_DEG             360.0                             
#define GD_PI_DIV_2_RAD         1.5707963267948966192313216916398   // PI/2 in radians
#define GD_PI_DIV_2_DEG         90.0                                // PI/2 in degree
#define GD_PI_DIV_3_RAD         1.0471975511965977461542144610932   // PI/3 in radians
#define GD_PI_DIV_3_DEG         60.0                                // PI/3 in degree
#define GD_3_PI_DIV_2_RAD       4.7123889803846898576939650749193   // 3PI/2 in radians
#define GD_3_PI_DIV_2_DEG       270.0                               // 3PI/2 in degree
#define GD_PI_DIV_4_RAD         0.78539816339744830961566084581988  // PI/4 in radians
#define GD_PI_DIV_4_DEG         45.0                                // PI/4 in degree
#define GD_PI_DIV_10_RAD        0.31415926535897932384626433832795  // PI/10 in radians
#define GD_PI_DIV_10_DEG        18.0                                // PI/10 in degree
#define GD_PI_DIV_20_RAD        0.15707963267948966192313216916398  // PI/20 in radians
#define GD_PI_DIV_20_DEG        9.0                                 // PI/20 in degree
#define GD_PI_DIV_40_RAD        7.85398163397448309615660845819e-2  // PI/40 in radians
#define GD_PI_DIV_40_DEG        4.5                                 // PI/40 in degree

// Upper latitude ranges when distance steps should be scaled
#define GD_UPPERSCALE_LAT_RAD   1.5620696805349249713467032377973
#define GD_UPPERSCALE_LAT_DEG   89.5

// Another mathematic constants
#define GD_GOLDENSECTION_LO     0.38196601125010515179541316563436  // Lower golden section part
#define GD_GOLDENSECTION_HI     0.61803398874989484820458683436564  // Higher golden section part

// Degree<->Radians conversion
#define GD_DEG_IN_RAD           57.295779513082320876798154814105   // One radian in degree
#define GD_RAD_IN_DEG           1.74532925199432957692222222222e-2  // One degree in radians                                
#define StdAngle(a,flg,t,f)     ((flg) ? a##_##t : a##_##f)
#define Radian2Degree(Rad)      ((Rad) * GD_DEG_IN_RAD)
#define Degree2Radian(Deg)      ((Deg) * GD_RAD_IN_DEG)

// Radians<->Nautical miles<->Meters conversion
#define GD_RAD_IN_MILE          2.90888208665721596153703703703e-4
#define GD_MILES_IN_RAD         3437.7467707849392526107818606515
#define GD_DEG_IN_MILE          1.66666666666666666666666666666e-2
#define GD_MILES_IN_DEG         60.0
#define GD_METERS_IN_MILE       1852.0
#define GD_MILES_IN_METER       5.39956803455723542116630669546e-4
#define Miles2Radian(Miles)     ((Miles) * GD_RAD_IN_MILE)
#define Radian2Miles(Rad)       ((Rad) * GD_MILES_IN_RAD)
#define Miles2Degree(Miles)     ((Miles) * GD_DEG_IN_MILE)
#define Degree2Miles(Deg)       ((Deg) * GD_MILES_IN_RAD)
#define CvMiles2Radian(Miles)   ((Miles) *= GD_RAD_IN_MILE)
#define CvRadian2Miles(Rad)     ((Rad) *= GD_MILES_IN_RAD)
#define CvMiles2Degree(Miles)   ((Miles) *= GD_DEG_IN_MILE)
#define CvDegree2Miles(Deg)     ((Deg) *= GD_MILES_IN_RAD)
#define Miles2Meters(Miles)     ((Miles) * GD_METERS_IN_MILE)
#define Meters2Miles(Meters)    ((Meters) * GD_MILES_IN_METER)
#define CvMiles2Meters(Miles)   ((Miles) *= GD_METERS_IN_MILE)
#define CvMeters2Miles(Meters)  ((Meters) *= GD_MILES_IN_METER)

// Default longitude/latitude definitions, world sides
#define GD_LATITUDE_FORMAT      "dd_mm.mmmS"
#define GD_LONGITUDE_FORMAT     "ddd_mm.mmmS"
#define GD_LATITUDE_SIZE        ARR_SIZE (GD_LATITUDE_FORMAT)
#define GD_LONGITUDE_SIZE       ARR_SIZE (GD_LONGITUDE_FORMAT)
#define GD_NORTH_SIDE           'N'
#define GD_SOUTH_SIDE           'S'
#define GD_EAST_SIDE            'E'
#define GD_WEST_SIDE            'W'
#define GD_NORTHLAT_DEG         GD_PI_DIV_2_DEG
#define GD_SOUTHLAT_DEG         - GD_PI_DIV_2_DEG
#define GD_MAXLONG_DEG          GD_PI_DEG
#define GD_NORTHLAT_RAD         GD_PI_DIV_2_RAD
#define GD_SOUTHLAT_RAD         - GD_PI_DIV_2_RAD
#define GD_MAXLONG_RAD          GD_PI_RAD
#define GD_EQUATOR              GD_ZERO

// Config profile settings
#define GD_LAT_FORMAT_KEY       "Latitude Format"
#define GD_LONG_FORMAT_KEY      "Longitude Format"
#define GD_ERROR_MSG_KEY        "Error Messages"
#define GD_BREAK_ON_ERROR_KEY   "Break On Error"
#define GD_ASSUME_RAD_KEY       "Assume Radian Args"
#define GD_PRECISION_KEY        "Iterations Precision"
#define GD_STEPDISTANCE_KEY     "Step Distance"

// Earth properties (WGS84)
#define GD_WGS84_EQUAT_RAD_M    6378137.0
#define GD_WGS84_POLAR_RAD_M    6356752.3142451793
#define GD_WGS84_FLATTENING     3.35281066474751169502944198282e-3
#define GD_SPHERE_RAD_M         6366707.0194937074958298109629434

// fmod() redefinition (Microsoft-specific fmod() returns remainder with the same sign as the 
// divident; this is a mathematical error!
#ifdef  fmod
#undef  fmod
#endif
#define fmod(x,y)               ((x) >= GD_ZERO ? fmod (x, y) : fmod (x, y) + (y))

// Validation of the double arg
#define IsInvalidGeoVal(val)    (_fpclass ((double) (val)) & (_FPCLASS_SNAN | _FPCLASS_QNAN | _FPCLASS_NINF | _FPCLASS_PINF))

// Angle normalize macros
#define NormalizeLong(a,u)      { while (a < - GD_PI_##u) a += GD_2_PI_##u; while (a > GD_PI_##u) a -= GD_2_PI_##u; }
#define NormalizeLongNoPi(a,u)  { while (a < - GD_PI_##u) a += GD_2_PI_##u; while (a >= GD_PI_##u) a -= GD_2_PI_##u; }
#define NormalizeBearing(a,u)   { while (a < 0.0) a += GD_2_PI_##u; while (a >= GD_2_PI_##u) a -= GD_2_PI_##u; }
#define NormalizeArcAngle(a,u)  { while (a < 0.0) a += GD_2_PI_##u; while (a > GD_2_PI_##u) a -= GD_2_PI_##u; }

// Safe by-ref-vars assignment
#define SafeAssign(ptr,expr)    if (ptr) try         \
                                {                    \
                                    *(ptr) = (expr); \
                                }                    \
                                catch (...)          \
                                {                    \
                                    /* Do nothing */ \
                                }

// Safe heap memory block freeing
#define SafeFree(ptr)           if (ptr) try         \
                                {                    \
                                    free (ptr);      \
                                }                    \
                                catch (...)          \
                                {                    \
                                    /* Do nothing */ \
                                }

// Course turn-over (circulation)
#define CheckCourse(crs,unit)   while (crs < GD_ZERO) crs += GD_2_PI_##unit
#define TurnCourse(course,unit) {                                       \
                                    course = ((course) > GD_PI_##unit ? \
                                        (course) - GD_PI_##unit :       \
                                        (course) + GD_PI_##unit);       \
                                                                        \
                                    NormalizeBearing (course, unit);    \
                                }

// Longitude cyclic check & correction
#define CheckLongitude(l,unit)  if (l < - GD_MAXLONG_##unit)      \
                                    l += GD_2_PI_##unit;        \
                                else if (l > GD_MAXLONG_##unit) \
                                    l -= GD_2_PI_##unit

// Misc definitions and macros
#define Sgn(val)                (((val) > GD_ZERO) ? 1.0 : -1.0 )
#define ISgn(val)               (((val) > GD_IZERO) ? 1 : -1)
#define Str2GeoVal(str)         (atof (str))
#define Odd(val)                ((val) & 1)
#define Even(val)               (Odd(val) == FALSE)
#define InRange(val,beg,end)    ((val) < (end) && (val) > (beg))
#define InRange2(val,beg,end)   (InRange(val,beg,end) || InRange(val,end,beg))

// Indirect comparision two double values is a great proramming error but the runtime does not
// work properly without this action where values are near the zero - I don't know why...
#define IsZero(val)             (((val) == GD_ZERO) || (((val) > -GD_ZERORANGE) && ((val) < GD_ZERORANGE)))
#define IsNotZero(val)          (((val) != GD_ZERO) && (((val) < -GD_ZERORANGE) || ((val) > GD_ZERORANGE)))
#define CheckZero(val)          (IsZero (val) ? GD_ZERO : (val))
#define Abs(val)                ((val) == GD_ZERO ? GD_ZERO : ((val) > GD_ZERORANGE ? (val) : ((val) < - GD_ZERORANGE ? - (val) : 0.0)))
#define IsSame(val1,val2)       ((val1) == (val2))
#define IsNotSame(val1,val2)    ((val1) != (val2))
#define Difference(val1,val2)   ((val1) - (val2))
#define IsEqual(val1,val2)      (                                               \
                                 IsSame (val1, val2) ||                         \
                                 (                                              \
                                  (Difference (val1, val2) < GD_ZERORANGE) &&   \
                                  (Difference (val1, val2) > -GD_ZERORANGE)     \
                                 )                                              \
                                )
#define IsNotEqual(val1,val2)   (                                               \
                                 IsNotSame (val1, val2) &&                      \
                                 (                                              \
                                  (Difference (val1, val2) >= GD_ZERORANGE) ||  \
                                  (Difference (val1, val2) <= -GD_ZERORANGE)    \
                                 )                                              \
                                )
#define Half(val)               ((val) * 0.5)
#define Two(val)                ((val) + (val))
#define ChangeSign(var)         (var = - var)
#define Square(val)             ((val) * (val))
#define HypoLen(cat1,cat2)      sqrt (Square (cat1) + Square (cat2))
//#define CvRadians(val)        (val * GD_RAD_IN_DEG)
//#define CvDegree(val)         (val * GD_DEG_IN_RAD)
#define ToRadians(val)          (val) *= GD_RAD_IN_DEG
#define ToDegree(val)           (val) *= GD_DEG_IN_RAD
#define WrongLat(val,unit)      ((val) < GD_SOUTHLAT_##unit || (val) > GD_NORTHLAT_##unit) // unit 
#define WrongLong(val,unit)     ((val) < -GD_MAXLONG_##unit || (val) > GD_MAXLONG_##unit)  // should be
#define WrongCourse(val,unit)   ((val) < GD_ZERO || (val) > GD_2_PI_##unit)                // DEG or RAD
#define WrongDistance(dist)     ((dist) < GD_ZERO)

// App-related functions
HINSTANCE GetInstanceHandle (void);

#endif

// **************************************************************************************************
// Copy data from Geo Kernel Errors.h
//
#ifndef _MARIS_GEOMETRY_ERRORCODES_
#define _MARIS_GEOMETRY_ERRORCODES_

#define GKERR_WRONGLATFORMAT            1                   // Wrong latitude format
#define GKERR_WRONGLONGFORMAT           2                   // Wrong longitude format
#define GKERR_LATOUTOFRNG               3                   // Latitude value out of range
#define GKERR_LONGOUTOFRNG              4                   // Longitude value out of range
#define GKERR_INVALIDVALUE              5                   // Invalid value
#define GKERR_VALOUTOFRNG               6                   // Value out of range
#define GKERR_NOTAPPLICABLEFORPOLE      7                   // Function is not applicable for pole
#define GKERR_COURSEOUTOFRNG            8                   // Course out of range
#define GKERR_INVALIDDISTANCE           9                   // Invalid distance
#define GKERR_DATUMNOTSUPPORTED         10                  // Datum is not supported
#define GKERR_SEGCOLLAPSEDTOPOINT       11                  // Line/segment is collapsed to point
#define GKERR_TOOFEWPOINTS              12                  // Too few points in the contour
#define GKERR_MEMALLOCFAILED            13                  // Memory allocation failed
#define GKERR_TOOMANYPOINTS             14                  // Too many points in the contour
#define GKERR_POINTDOESNTEXIST          15                  // Searching point does not exist
#define GKERR_FUNCTIONISNOTAPPLICABLE   16                  // Function is not applicable for this type of line

#define GKERR_MAXERRORS                 16

#endif

// **************************************************************************************************
// Copy data from Geo Kernel Types.h
//
#ifndef _MARIS_GEOMETRY_GLOBAL_TYPEDEFS_
#define _MARIS_GEOMETRY_GLOBAL_TYPEDEFS_

#ifdef _MARIS_GEOMETRY_INTERNAL_
#define GK_API  __declspec (dllexport)
#else
#define GK_API  __declspec (dllimport)
#endif

#define _USE_ECDIS_TYPES_
#ifdef _USE_ECDIS_TYPES_

#ifndef _GD_GEOPOINT_DEFINED_
#define _GD_GEOPOINT_DEFINED_

#pragma pack(1)

typedef struct
{
    double lat,
           lon;
}
SGeoPoint;

#pragma pack()

#endif // _GD_GEOPOINT_DEFINED_

#pragma pack(1)

typedef SGeoPoint AVector [2];

typedef SGeoPoint Gpoint_t;

#pragma pack()

#define gvLatitude      lat
#define gvLongitude     lon
#define GEOVAL          double
#define PGEOVAL         double *

#else // _USE_ECDIS_TYPES_

typedef double GEOVAL, *PGEOVAL;

#pragma pack(1)

typedef struct tagGeoPoint
{
    union
    {
        GEOVAL gvLatitude,
               lat;
    };
    union
    {
        GEOVAL gvLongitude,
               lon;
    };
}
SGeoPoint;

#pragma pack()

#endif // _USE_ECDIS_TYPES_

#pragma pack(1)

typedef union tagGeoRect
{
    struct
    {
        SGeoPoint nw,
                  se;
    };

    struct
    {
        double nw_lat,
               nw_lon,
               se_lat,
               se_lon;
    };
}
UGeoRect;

typedef UGeoRect Garea_t;

typedef union tagGeoSegment
{
    struct tagGeoSegmentCoord
    {
        GEOVAL      gvBegLatitude,
                    gvBegLongitude,
                    gvEndLatitude,
                    gvEndLongitude;
    } 
    Coord;

    struct tagGeoSegmentEnds
    {
        SGeoPoint   gptBegPoint,
                    gptEndPoint;
    } 
    EndPoints;
} 
UGeoSegment;

typedef struct tagGeoSegmentS
{
    SGeoPoint gpBegin;
    GEOVAL    gvBearing,
              gvDistance;
}
SGeoSegment;

typedef struct tagGeoArc
{
    SGeoPoint gpBegin,
              gpCenter;
    GEOVAL    gvAngle;
}
SGeoArc;

typedef struct tagMercatorPoint
{
    GEOVAL gvNorthing,
           gvEasting;
}
SMercatorPoint;

typedef union tagMercatorSegment
{
    struct tagMercatorSegmentCoord
    {
        GEOVAL          gvBegLatitude,
                        gvBegLongitude,
                        gvEndLatitude,
                        gvEndLongitude;
    } 
    Coord;

    struct tagMercatorSegmentEnds
    {
        SMercatorPoint  gptBegPoint,
                        gptEndPoint;
    } 
    EndPoints;
} 
UMercatorSegment;

typedef struct tagGeoContour
{
    int         nVertexCount;
    SGeoPoint  *pgptVertexes;
    GEOVAL      gvPhaseShift;
	int         tag;
}
SGeoContour;

typedef struct tagGeoRegion
{
    int           nExtContourCount,
                  nIntContourCount;
    GEOVAL        gvPhaseShift;                     // The same for all the contours
    SGeoContour  *pgcContours;
}
SGeoRegion;

typedef enum tagPointPrimRelation
{
    pprInside  = 1,                                 // Point is inside of the primitive
    pprOutside = 2,                                 // Point is outside of the primitive
    pprOnEdge  = 4,                                 // Point touches the edge of the primitive
    pprUnknown = 8                                  // Relation is unknown or error is detected
}
EPointPrimRelation;

typedef enum tagContourContourRelation
{
    ccrNoCommonParts = 1,                           // Two contours have no common parts
    ccrContains      = 2,                           // This contour contains another one
    ccrIncluded      = 4,                           // This contour included into another one
    ccrCrosses       = 8                            // This contour crosses into another one
}
EContourContourRelation;

typedef enum tagPointEntryTestMethod                // Method for test is point included into
{                                                   // contour
    petmCountCrosses = 1,                           // By count of crosses with edges
    petmByGCAngles   = 2,                           // By summ of angles between dirs to vertexes
    petmByRLAngles   = 3                            // (according great circle and rhumb line)
}
EPointEntryTestMethod;

typedef enum tagLineType
{
    ltSimple      = 1,
    ltRhumbLine   = 2,
    ltGreatCircle = 3
}
ELineType;

typedef enum tagDistCalcMethod
{
    dcmByShortestDistance = 1,
    dcmByPerpendicular    = 2
}
EDistCalcMethod;

#pragma pack()

#endif

// **************************************************************************************************
// Copy data from Datums.h
//
#ifndef _MARIS_GEOMETRY_DATUMS_
#define _MARIS_GEOMETRY_DATUMS_

// Geoid definitions, Internal ids
#define GD_GEOID_WGS84              1
#define GD_GEOID_WGS72              2
#define GD_GEOID_WGS60              3
#define GD_GEOID_WGS66              4
#define GD_GEOID_NAD27              5
#define GD_GEOID_GRS67              6
#define GD_GEOID_GRS80              7
#define GD_GEOID_KRASSOVSKY         8
#define GD_GEOID_AIRY1830           9
#define GD_GEOID_MODIFIEDAIRY       10
#define GD_GEOID_CLARKE1866         11
#define GD_GEOID_CLARKE1880         12
#define GD_GEOID_BESSEL1841         13
#define GD_GEOID_BESSEL1841_NAM     14
#define GD_GEOID_EVEREST_INDIA1830  15
#define GD_GEOID_EVEREST_SABAH      16
#define GD_GEOID_EVEREST_INDIA1956  17
#define GD_GEOID_AUSTRALIA_NAT      18
#define GD_GEOID_MALAYSIA1969       19
#define GD_GEOID_EVEREST_MALAY_SING 20
#define GD_GEOID_EVEREST_PAKISTAN   21
#define GD_GEOID_FISCHER1960        22
#define GD_GEOID_FISCHER1968        23
#define GD_GEOID_HELMERT1906        24
#define GD_GEOID_HOUGH              25
#define GD_GEOID_INDONESIAN1974     26
#define GD_GEOID_INTERNATIONAL1924  27
#define GD_GEOID_MODIFIEDEVEREST    28
#define GD_GEOID_MODIFIEDFISHER     29
#define GD_GEOID_SOUTHAMERICAN1969  30
#define GD_GEOID_SPHERE             0

#define GD_DATUM_UNKNOWN                                   0
#define GD_DATUM_WGS84                                     1
#define GD_DATUM_ADINDAN_BURKINA_FASO                      2
#define GD_DATUM_ADINDAN_CAMEROON                          3
#define GD_DATUM_ADINDAN_ETHIOPIA                          4
#define GD_DATUM_ADINDAN_MALI                              5
#define GD_DATUM_ADINDAN_MEAN_EPTH_SUDAN                   6
#define GD_DATUM_ADINDAN_SENEGAL                           7
#define GD_DATUM_ADINDAN_SUDAN                             8
#define GD_DATUM_AFGOOYE                                   9
#define GD_DATUM_AIN_EL_ABD_1970_BAHRAIN                   10
#define GD_DATUM_AIN_EL_ABD_1970_SAUDARABIA                11
#define GD_DATUM_AMER_SAMOA_1962                           12
#define GD_DATUM_ANNA_1_ASTRO_1965                         13
#define GD_DATUM_ANTIQ_ASTRO_1943                          14
#define GD_DATUM_ARC_1950_BOTSWANA                         15
#define GD_DATUM_ARC_1950_BURUNDI                          16
#define GD_DATUM_ARC_1950_LESOTO                           17
#define GD_DATUM_ARC_1950_MALAWI                           18
#define GD_DATUM_ARC_1950_MEAN                             19
#define GD_DATUM_ARC_1950_SWAZILAND                        20
#define GD_DATUM_ARC_1950_ZAIRE                            21
#define GD_DATUM_ARC_1950_ZAMBIA                           22
#define GD_DATUM_ARC_1950_ZIMBABWE                         23
#define GD_DATUM_ARC_1960_MEAN                             24
#define GD_DATUM_ARC_1960_KENIA                            25
#define GD_DATUM_ARC_1960_TANZANIA                         26
#define GD_DATUM_ASCENSION_ISLAND_1958                     27
#define GD_DATUM_ASTRO_BEACON_E                            28
#define GD_DATUM_ASTRO_POS_71_74                           29
#define GD_DATUM_ASTRO_TERN_ISLAND_1961                    30
#define GD_DATUM_ASTRONOMIC_STATION_1952                   31
#define GD_DATUM_AUSTRALIAN_GEODETIC_1966                  32
#define GD_DATUM_AUSTRALIAN_GEODETIC_1984                  33
#define GD_DATUM_AYABELLE_LIGHTHOUSE                       34
#define GD_DATUM_BELLEVUE_IGN                              35
#define GD_DATUM_BERMUDA_1957                              36
#define GD_DATUM_BISSAU                                    37
#define GD_DATUM_BOGOTA_OBSERVATORY                        38
#define GD_DATUM_BUKIT_RIMPAH                              39
#define GD_DATUM_CAMP_AREA_ASTRO                           40
#define GD_DATUM_CAMPO_INCHAUSPE                           41
#define GD_DATUM_CANTON_ISLAND_1966                        42
#define GD_DATUM_CAPE                                      43
#define GD_DATUM_CAPE_CANAVERAL_MEAN                       44
#define GD_DATUM_CARTHAGE                                  45
#define GD_DATUM_CHATHAM_1971                              46
#define GD_DATUM_CHUA_ASTRO                                47
#define GD_DATUM_CORREGO_ALEGRE                            48
#define GD_DATUM_DABOLA                                    49
#define GD_DATUM_DECEPTIONISLAND                           50
#define GD_DATUM_DJAKARTA_BATAVIA                          51
#define GD_DATUM_DOS_1968                                  52
#define GD_DATUM_EASTER_LSLAND_1967                        53
#define GD_DATUM_DAI_ESTONIA_1937                          54
#define GD_DATUM_ED_50_CYPRUS                              55
#define GD_DATUM_ED_50_EGYPT                               56
#define GD_DATUM_ED_50_NORTH_ENGLAND                       57
#define GD_DATUM_ED_50_NW_ENGLAND                          58
#define GD_DATUM_ED_50_FINLAND_NORWAY                      59
#define GD_DATUM_ED_50_GREECE                              60
#define GD_DATUM_ED_50_IRAN                                61
#define GD_DATUM_ED_50_SARDINIA                            62
#define GD_DATUM_ED_50_SICILY                              63
#define GD_DATUM_ED_50_MALTA                               64
#define GD_DATUM_ED_50_MEAN                                65
#define GD_DATUM_ED_50_MEAN_FOR_CENTER                     66
#define GD_DATUM_ED_50_MEAN_FOR_EAST                       67
#define GD_DATUM_ED_50_PORTUGAL_SPAIN                      68
#define GD_DATUM_ED_50_TUNISIA                             69
#define GD_DATUM_EUROPEAN_1979_MEAN                        70
#define GD_DATUM_FORT_THOMAS_1955                          71
#define GD_DATUM_GAN_1970                                  72
#define GD_DATUM_GEODETIC_DATUM_1949                       73
#define GD_DATUM_GRACIOSA_BASE_SW_1948                     74
#define GD_DATUM_GUAM_1963                                 75
#define GD_DATUM_GUNUNG_SEGARA                             76
#define GD_DATUM_GUX_1_ASTRO                               77
#define GD_DATUM_HERAT_NORTH                               78
#define GD_DATUM_HERMANNSKOGEL                             79
#define GD_DATUM_HJORSEY_1955                              80
#define GD_DATUM_HONG_KONG_1963                            81
#define GD_DATUM_HU_TZU_SHAN                               82
#define GD_DATUM_INDIAN_BANGLADESH                         83
#define GD_DATUM_INDIAN_INDIA_NEPAL                        84
#define GD_DATUM_INDIAN_PAKISTAN                           85
#define GD_DATUM_INDIAN_1954_THAILAND                      86
#define GD_DATUM_INDIAN_1960_CON_SON_ISLAND                87
#define GD_DATUM_INDIAN_1960_16oN                          88
#define GD_DATUM_INDIAN_1975_THAILAND                      89
#define GD_DATUM_INDONESIAN_1974                           90
#define GD_DATUM_IRELAND_1965                              91
#define GD_DATUM_ISTS_061_ASTRO_1968                       92
#define GD_DATUM_ISTS_073_ASTRO_1969                       93
#define GD_DATUM_JOHNSTON_ISLAND_1961                      94
#define GD_DATUM_KANDAWALA                                 95
#define GD_DATUM_KERGUELEN_ISLAND                          96
#define GD_DATUM_KERTAU_1948                               97
#define GD_DATUM_KUSAIE_ASTRO_1951                         98
#define GD_DATUM_KOREAN_GEODETIC_SYSTEM                    99
#define GD_DATUM_LC_5_ASTRO                                100
#define GD_DATUM_LEIGON                                    101
#define GD_DATUM_LIBERIA_1964                              102
#define GD_DATUM_LUSON_PHILIPINES                          103
#define GD_DATUM_LUSON_MINDANAO                            104
#define GD_DATUM_M_PORALOKO                                105
#define GD_DATUM_MAHE_1971                                 106
#define GD_DATUM_MASSAWA                                   107
#define GD_DATUM_MERCHICH                                  108
#define GD_DATUM_MIDWAY_ASTRO_1961                         109
#define GD_DATUM_MINNA_CAMERIIN                            110
#define GD_DATUM_MINNA_NIGERIA                             111
#define GD_DATUM_MONTSERRAT_ISLAND_ASTRO_58                112
#define GD_DATUM_MASIRAH_ISLAND_NAHRWAN                    113
#define GD_DATUM_NAHRWAN                                   114
#define GD_DATUM_UNITES_ARAB_EMIRATES_NAHRWAN              115
#define GD_DATUM_NAPARIMA_BWI                              116
#define GD_DATUM_NAD_27_ALASKA                             117
#define GD_DATUM_NAD_27_EAST_ALEUTAN                       118
#define GD_DATUM_NAD_27_WEST_ALEUTAN                       119
#define GD_DATUM_NAD_27_BAHAMAS                            120
#define GD_DATUM_NAD_27_SAN_SALVADOR_ISLND                 121
#define GD_DATUM_NAD_27_CANADA_ALB_BRITCOL                 122
#define GD_DATUM_NAD_27_CANADA_MANIT_ONTARIO               123
#define GD_DATUM_NAD_27_CANADA_OTHER                       124
#define GD_DATUM_NAD_27_NW_CANADA                          125
#define GD_DATUM_NAD_27_CANADA_YUKON                       126
#define GD_DATUM_NAD_27_CANAL_ZONE                         127
#define GD_DATUM_NAD_27_CUBA                               128
#define GD_DATUM_NAD_27_GREENLAND                          129
#define GD_DATUM_NAD_27_CENTRAL_AMERICA                    130
#define GD_DATUM_NAD_27_MEAN_FOR_BELIZE                    131
#define GD_DATUM_NAD_27_CANADA_MEAN                        132
#define GD_DATUM_NAD_27_MEAN_FOR_CONUS                     133
#define GD_DATUM_NAD_27_MEAN_FOR_CONUS_E                   134
#define GD_DATUM_NAD_27_MEAN_FOR_CONUS_W                   135
#define GD_DATUM_NAD_27_MEXICO                             136
#define GD_DATUM_NAD_83_ALASKA                             137
#define GD_DATUM_NAD_83_ALEUTAN                            138
#define GD_DATUM_NAD_83_CANADA                             139
#define GD_DATUM_NAD_83_CONUS                              140
#define GD_DATUM_NAD_83_HAWAII                             141
#define GD_DATUM_NAD_83_MEXICO                             142
#define GD_DATUM_NORTH_SAHARA_1959                         143
#define GD_DATUM_OBSERVATORY_WET_1939                      144
#define GD_DATUM_OLD_EGYPTIAN                              145
#define GD_DATUM_OLD_HAWAIIAN                              146
#define GD_DATUM_OLD_HAWAIIAN_KAUAI                        147
#define GD_DATUM_OLD_HAWAIIAN_MAUI                         148
#define GD_DATUM_OLD_HAWAIIAN_MEAN                         149
#define GD_DATUM_OLD_HAWAIIAN_OAHU                         150
#define GD_DATUM_OMAN                                      151
#define GD_DATUM_ORDN_SURVEY_GB_1936                       152
#define GD_DATUM_ORDN_SURVEY_GB_1936_MAN_WALES             153
#define GD_DATUM_ORDN_SURVEY_GB_1936_MEAN                  154
#define GD_DATUM_ORDN_SURVEY_GB_1936_SCOTLAND              155
#define GD_DATUM_ORDN_SURVEY_GB_1936_WALES                 156
#define GD_DATUM_PICO_DE_LAS_NIEVES                        157
#define GD_DATUM_PITCAIRN_ASTRO_1967                       158
#define GD_DATUM_POINT_59                                  159
#define GD_DATUM_POINTE_NOIRE_1948                         160
#define GD_DATUM_PORTO_SANTO_1936                          161
#define GD_DATUM_PROVISIONAL_SAD_56_BOLIVIA                162
#define GD_DATUM_PROVISIONAL_SAD_56_NCHILE                 163
#define GD_DATUM_PROVISIONAL_SAD_56_SCHILE                 164
#define GD_DATUM_PROVISIONAL_SAD_56_COLUMBIA               165
#define GD_DATUM_PROVISIONAL_SAD_56_ECUADOR                166
#define GD_DATUM_PROVISIONAL_SAD_56_GUIANA                 167
#define GD_DATUM_PROVISIONAL_SAD_56_MEAN                   168
#define GD_DATUM_PROVISIONAL_SAD_56_PERU                   169
#define GD_DATUM_PROVISIONAL_SAD_56_VENEZUELA              170
#define GD_DATUM_PROVISIONAL_SOUTH_CHILEAN_1963            171
#define GD_DATUM_PUERTO_RICO                               172
#define GD_DATUM_PULKOVO_1942                              173
#define GD_DATUM_QUATAR_NATIONAL                           174
#define GD_DATUM_QORNOQ                                    175
#define GD_DATUM_REUNION                                   176
#define GD_DATUM_ROME_1940                                 177
#define GD_DATUM_S_42_HUNGARY                              178
#define GD_DATUM_S_42_POLAND                               179
#define GD_DATUM_S_42_CZECHOSLAVAKIA                       180
#define GD_DATUM_S_42_LATVIA                               181
#define GD_DATUM_S_42_KAZAKHSTAN                           182
#define GD_DATUM_S_42_ALBANIA                              183
#define GD_DATUM_S_42_ROMANIA                              184
#define GD_DATUM_S_JTSK                                    185
#define GD_DATUM_SANTO_DOS                                 186
#define GD_DATUM_SAO_BRAZ                                  187
#define GD_DATUM_SAPPER_HILL_1943                          188
#define GD_DATUM_SCHWARZECK                                189
#define GD_DATUM_SELVAGEM_GRANDE_1938                      190
#define GD_DATUM_SIERRA_LEONE_1960                         191
#define GD_DATUM_SAD_69_ARGENTINA                          192
#define GD_DATUM_SAD_69_BOLIVIA                            193
#define GD_DATUM_SAD_69_BRAZIL                             194
#define GD_DATUM_SAD_69_CHILE                              195
#define GD_DATUM_SAD_69_COLOMBIA                           196
#define GD_DATUM_SAD_69_ECUADOR                            197
#define GD_DATUM_SAD_69_ECUADOR_BALTRA_GALAP               198
#define GD_DATUM_SAD_69_GUYANA                             199
#define GD_DATUM_SAD_69_MEAN                               200
#define GD_DATUM_SAD_69_PARAGUAY                           201
#define GD_DATUM_SAD_69_PERU                               202
#define GD_DATUM_SAD_69_TRINIDAD_TOBAGO                    203
#define GD_DATUM_SAD_69_VENEZUELA                          204
#define GD_DATUM_SOUTH_ASIA                                205
#define GD_DATUM_TANANARIVE_OBSERVATORY_1925               206
#define GD_DATUM_TIMBALAI_1948                             207
#define GD_DATUM_TOKYO_JAPAN                               208
#define GD_DATUM_TOKYO_MEAN                                209
#define GD_DATUM_TOKYO_OKINAWA                             210
#define GD_DATUM_TOKYO_SKOREA                              211
#define GD_DATUM_TRISTAN_ASTRO_1968                        212
#define GD_DATUM_VITI_LEVU_1916                            213
#define GD_DATUM_VOIROL_1960                               214
#define GD_DATUM_WAKE_ISLAND_ASTRO_1952                    215
#define GD_DATUM_WAKE_ENIWETOK_1960                        216
#define GD_DATUM_WGS_72                                    217
#define GD_DATUM_YACARE                                    218
#define GD_DATUM_ZANDERIJ                                  219
#define GD_DATUM_RT90_RT70                                 220
#define GD_DATUM_NORSKGRADNET                              221

#pragma pack(1)

typedef struct tagGeoid
{
    int    nGeoidID;
    double dblEquatorialRadius,
           dblPolarRadius,
           dblFlattening;
}
SGeoid;

typedef struct tagDatumInfo
{
    int     nGeoidID,
            nDatumID;
    short   nXShift,
            nYShift,
            nZShift;
    TCHAR  *pszName;
}
SDatumInfo;

#pragma pack()

BOOL GetGeoid (int nDatumID, SGeoid *pdatDatum);

TCHAR GK_API *gkGetDatumName (int nDatumID);

BOOL GK_API gkdCalcDatumShifts (int              nDatumID, 
                                const SGeoPoint *pgpPoint, 
                                double&          dblDeltaLat, 
                                double&          dblDeltaLon);
BOOL GK_API gkrCalcDatumShifts (int              nDatumID, 
                                const SGeoPoint *pgpPoint, 
                                double&          dblDeltaLat, 
                                double&          dblDeltaLon);

BOOL GK_API gkdCalcDatumShifts (int              nDatumID, 
                                const SGeoPoint *pgpPoint, 
                                double&          dblDeltaLat, 
                                double&          dblDeltaLon);
BOOL GK_API gkrCalcDatumShifts (int              nDatumID, 
                                const SGeoPoint *pgpPoint, 
                                double&          dblDeltaLat, 
                                double&          dblDeltaLon);
#endif // _MARIS_GEOMETRY_DATUMS_

// ********************************************************
// Copy from GeoValueFormat.h
//
#ifndef _MARIS_GEOMETRY_GEOVALFMT_
#define _MARIS_GEOMETRY_GEOVALFMT_

#if !defined(AFX_GEOVALUEFORMAT_H__A0B280E0_2D86_11D4_9C13_0090271CF0EA__INCLUDED_)
#define AFX_GEOVALUEFORMAT_H__A0B280E0_2D86_11D4_9C13_0090271CF0EA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CGeoValueFormat  
{
    private:
        BOOL    m_bPrepared;
        UINT    m_uiWrongFmtErrCode,
                m_uiValOutOfRngErrCode;
        GEOVAL  m_gvMaxValue,
                m_gvMinValue;
             
        struct tagData
        {
            unsigned short  usDegreePos,
                            usDegreeSize,
                            usMinutesIntPartPos,
                            usMinutesIntPartSize,
                            usMinutesFactPartPos,
                            usMinutesFractPartSize,
                            usMinutesSize,
                            usWorldSidePos,
                            usGeneralSize;
            TCHAR           cPositiveWorldSide,
                            cNegativeWorldSide,
                            szDegreeFormat [10],
                            szMinuteFormat [10];
            LPTSTR          lpszFormat;

        } m_Data;

        void ClearData      ();
        void FillByAsterisk (LPTSTR lpszBuffer, int nMaxSize);

    public:
	    CGeoValueFormat (UINT   uiWrongFmtErrCode, 
                         UINT   uiValOutOfRngErrCode, 
                         GEOVAL gvMinValue, 
                         GEOVAL gvMaxValue);
	    virtual ~CGeoValueFormat ();

        BOOL   IsPrepared      () { return m_bPrepared; }
        void   Unprepare       () { m_bPrepared = FALSE; }
        BOOL   Prepare         (LPTSTR lpszFormat, TCHAR cPositiveWorldSide, TCHAR cNegativeWorldSide);
        BOOL   FormatValue     (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize, BOOL bAssumeRadians = FALSE);
        BOOL   FormatValue     (GEOVAL gvValue, LPTSTR  lpszDegBuffer, LPTSTR lpszMinBuffer, TCHAR *pcWorldSideChar, int nMaxDegSize, int nMaxMinSize, BOOL bAssumeRadians = FALSE);
        void   GetData         (LPTSTR *ppszFormat, TCHAR *pcPositiveWorldSide, TCHAR *pcNegativeWorldSide);
        LPTSTR GetMinuteFormat () { return m_Data.szMinuteFormat; }
};

#endif // !defined(AFX_GEOVALUEFORMAT_H__A0B280E0_2D86_11D4_9C13_0090271CF0EA__INCLUDED_)
#endif

// ********************************************************
// Copy from Maris Geometry Utilities.h
//
#ifndef _MARIS_GEOMETRY_CLASS_
#define _MARIS_GEOMETRY_CLASS_

#define _MARIS_GEOMETRY_INTERNAL_

class CMarisGeometryUtilities
{
    private:
        HINSTANCE        m_hInstance;
        BOOL             m_bShowErrorMessages,      // Shalls the system show error boxes?
                         m_bBreakOnError,           // Shalls the system break when error detected?
                         m_bAssumeRadianArgs;       // Assume radian arg values
        GEOVAL           m_gvPrecision,             // Interative calculations precision
                         m_gvStepDistance;          // Length of single interval for rhumb-line
                                                    // forward interpolation
        CGeoValueFormat *m_pgvfLatitude,            // ECDIS format parameters
                        *m_pgvfLongitude;

    public:
	    CMarisGeometryUtilities  ();

        // Access to m_hInstance property field
        HINSTANCE GetLibInstance ()                    { return m_hInstance; }
        void SetLibInstance      (HINSTANCE hInstance) { m_hInstance = hInstance; }

        // Load/save all kernel settings from INI-file
        void LoadSettings        (LPCTSTR lpszProfileName, LPCTSTR lpszSection);
        void SaveSettings        (LPCTSTR lpszProfileName, LPCTSTR lpszSection);

        // Configure longitude-latitude formatting
        BOOL SetLatitudeFormat   (LPTSTR lpszFormat, TCHAR cNorthChar, TCHAR cSouthChar);
        BOOL SetLongitudeFormat  (LPTSTR lpszFormat, TCHAR cEastChar, TCHAR cWestChar);
        BOOL GetLatitudeFormat   (LPTSTR *ppszFormat, TCHAR *pcNorthChar, TCHAR *pcSouthChar);
        BOOL GetLongitudeFormat  (LPTSTR *ppszFormat, TCHAR *pcEastChar, TCHAR *pcWestChar);
        BOOL GetUnpreparedFormat (CGeoValueFormat& gvfFormat, LPTSTR *ppszData);

        // Separate minute formatting
        LPTSTR GetLatMinuteFormat  () { return m_pgvfLatitude->GetMinuteFormat  (); }
        LPTSTR GetLongMinuteFormat () { return m_pgvfLongitude->GetMinuteFormat (); }

        // Error handling interface
        void SetUnhandledError   (UINT uiNewErrorCode);
        void ClearError          ();
        void ErrorMessage        ();
        void ErrorMessage        (LPTSTR lpszErrorMsg);
        void ErrorMessage        (UINT uiErrorMsg);

        // User interface settings
        BOOL GetErrorMsgFlag     ()                  { return m_bShowErrorMessages; }
        BOOL SetErrorMsgFlag     (BOOL bFlag = TRUE) { return m_bShowErrorMessages = bFlag; }
        void SetDegreeMode       ()                  { m_bAssumeRadianArgs = FALSE; }
        void SetRadianMode       ()                  { m_bAssumeRadianArgs = TRUE; }
        BOOL IsDegreeAssumed     ()                  { return !m_bAssumeRadianArgs; }
        BOOL IsRadiansAssumed    ()                  { return m_bAssumeRadianArgs; }

        // Debugging support
        BOOL GetBreakOnErrorFlag ()                  { return m_bBreakOnError; }
        BOOL SetBreakOnErrorFlag (BOOL bFlag = TRUE) { return m_bBreakOnError = bFlag; }

        // Calculation control
        GEOVAL GetPrecision      ()                  { return m_gvPrecision; }
        void   SetPrecision      (GEOVAL gvPrecis)   { m_gvPrecision = gvPrecis; }
        GEOVAL GetStepDistance   ()                  { return m_gvStepDistance; }
        void   SetStepDistance   (GEOVAL gvDistance) { m_gvStepDistance = gvDistance; }

        // Geographical cooredinates formatting
        BOOL FormatLatitude      (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize) { return m_pgvfLatitude->FormatValue (gvValue, lpszBuffer, nMaxSize, m_bAssumeRadianArgs); }
        BOOL FormatLatitude      (GEOVAL gvValue, LPTSTR lpszDegBuffer, LPTSTR lpszMinBuffer, TCHAR *pcWorldSideChar, int nMaxDegSize, int nMaxMinSize) { return m_pgvfLatitude->FormatValue (gvValue, lpszDegBuffer, lpszMinBuffer, pcWorldSideChar, nMaxDegSize, nMaxMinSize, m_bAssumeRadianArgs); }
        BOOL FormatLongitide     (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize) { return m_pgvfLongitude->FormatValue (gvValue, lpszBuffer, nMaxSize, m_bAssumeRadianArgs); }
        BOOL FormatLongitide     (GEOVAL gvValue, LPTSTR lpszDegBuffer, LPTSTR lpszMinBuffer, TCHAR *pcWorldSideChar, int nMaxDegSize, int nMaxMinSize) { return m_pgvfLatitude->FormatValue (gvValue, lpszDegBuffer, lpszMinBuffer, pcWorldSideChar, nMaxDegSize, nMaxMinSize, m_bAssumeRadianArgs); }
};

// Global functions
void GK_API gkDeInitialize ();

CMarisGeometryUtilities& GetLibKernel ();

#endif

// Radian/degree args interpretation mode switching
void GK_API gkSetDegreeMode       (void);
void GK_API gkSetRadianMode       (void);
BOOL GK_API gkIsDegreeAssumed     (void);
BOOL GK_API gkIsRadiansAssumed    (void);
void GK_API gkSetRadianAssumeFlag (BOOL bFlag = TRUE);

// To read description of function purpose see the file "Maris Geometry Utilities.h"

// Common purpose calls
HINSTANCE GK_API gkGetLibInstance (void);

// Error handling
void   GK_API gkClearError    (void);
void   GK_API gkSetError      (UINT uiErrorCode);
UINT   GK_API gkGetErrorCode  (void);
UINT   GK_API gkGetErrorId    (void);
LPTSTR GK_API gkGetError      (LPTSTR lpszBuffer, int nMaxSize);

// Calculation control
GEOVAL GK_API gkGetPrecision    (void);
void   GK_API gkSetPrecision    (GEOVAL gvPrecis);
GEOVAL GK_API gkGetStepDistance (void);
void   GK_API gkSetStepDistance (GEOVAL gvDistance);

// Load/save config from/to INI file
void GK_API gkLoadSettings (LPCTSTR lpszProfileName, LPCTSTR lpszSection);
void GK_API gkSaveSettings (LPCTSTR lpszProfileName, LPCTSTR lpszSection);

// Configure longitude-latitude formatting
BOOL GK_API gkSetLatitudeFormat  (LPTSTR lpszFormat, TCHAR cNorthChar, TCHAR cSouthChar);
BOOL GK_API gkSetLongitudeFormat (LPTSTR lpszFormat, TCHAR cEastChar, TCHAR cWestChar);
BOOL GK_API gkGetLatitudeFormat  (LPTSTR *ppszFormat, TCHAR *pcNorthChar, TCHAR *cSouthChar);
BOOL GK_API gkGetLongitudeFormat (LPTSTR *ppszFormat, TCHAR *pcEastChar, TCHAR *pcWestChar);

// User interface parameters
BOOL GK_API gkGetErrorMsgFlag (void);                  // Return TRUE if error msg will be shown
BOOL GK_API gkSetErrorMsgFlag (BOOL bFlag = TRUE);     // Set/clear this flag, return new value

// Debugging support
BOOL GK_API gkGetBreakOnErrorFlag (void);              // Return TRUE if execution will be broken on error
BOOL GK_API gkSetBreakOnErrorFlag (BOOL bFlag = TRUE); // Set/clear this flag, return new value

#ifndef _DLGFUNC_DEFINED_
#define _DLGFUNC_DEFINED_

// Dialog functions
BOOL GK_API gkQueryGeoPoint (GEOVAL    *pgvLatitude,
                             GEOVAL    *pgvLongitude,
                             HWND       hwndParent = 0 /*HWND_DESKTOP*/,
                             LPCTSTR    lpszTitle  = NULL,
                             LPCTSTR    lpszPrompt = NULL);
BOOL GK_API gkQueryGeoPoint (SGeoPoint *pgptPoint,
                             HWND       hwndParent = 0 /*HWND_DESKTOP*/,
                             LPCTSTR    lpszTitle  = NULL,
                             LPCTSTR    lpszPrompt = NULL);

#endif // _DLGFUNC_DEFINED_

// Formatting/deformatting functions
#ifndef _GK_FORMATTING_
#define _GK_FORMATTING_

BOOL GK_API gkFormatLatitude (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize);
BOOL GK_API gkFormatLongitude (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize); 

BOOL GK_API gkAutoFormatLatitude (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize);
BOOL GK_API gkAutoFormatLongitude (GEOVAL gvValue, LPTSTR lpszBuffer, int nMaxSize); 

BOOL GK_API gkDeFormatLatitudeFree (LPTSTR lpszLatitude, GEOVAL *pgvLatitude);
BOOL GK_API gkDeFormatLongitudeFree (LPTSTR lpszLongitude, GEOVAL *pgvLongitude); 

BOOL GK_API gkDeFormatGeoPoint (LPTSTR lpszLatitude, LPTSTR lpszLongitude, GEOVAL *pgvLatitude, GEOVAL *pgvLongitude);
BOOL GK_API gkDeFormatGeoPoint (LPTSTR lpszLatitude, LPTSTR lpszLongitude, SGeoPoint *pgpPoint);

#endif // _GK_FORMATTING_

#ifndef _UNIVERSAL_FUNC_DEFINED_
#define _UNIVERSAL_FUNC_DEFINED_

// Universal direct/backward calculations
BOOL GK_API gkrCalcDistAndBrg (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgptBegin, 
                               SGeoPoint *pgptEnd,
                               GEOVAL    *pgvDistance,
                               GEOVAL    *pgvBegCourse,
                               GEOVAL    *pgvEndCourse,
                               BOOL       bUseLongestArc = FALSE);
BOOL GK_API gkrCalcDistAndBrg (int        nDatum,
                               ELineType  ltLineType,
                               GEOVAL     gvBegLat,
                               GEOVAL     gvBegLon,
                               GEOVAL     gvEndLat,
                               GEOVAL     gvEndLon,
                               GEOVAL    *pgvDistance,
                               GEOVAL    *pgvBegCourse,
                               GEOVAL    *pgvEndCourse,
                               BOOL       bUseLongestArc = FALSE);
BOOL GK_API gkrCalcPos        (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgptBegin, 
                               GEOVAL     gvDistance,
                               GEOVAL     gvBegCourse,
                               SGeoPoint *pgptEnd,
                               GEOVAL    *pgvEndCourse);
BOOL GK_API gkrCalcPos        (int        nDatum,
                               ELineType  ltLineType,
                               GEOVAL     gvBegLat,
                               GEOVAL     gvBegLon,
                               GEOVAL     gvDistance,
                               GEOVAL     gvBegCourse,
                               GEOVAL    *pgvEndLat,
                               GEOVAL    *pgvEndLon,
                               GEOVAL    *pgvEndCourse);
BOOL GK_API gkrCalcPosInv     (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgpBegin, 
                               GEOVAL     gvDistance,
                               GEOVAL     gvEndCourse,
                               SGeoPoint *pgpEnd,
                               GEOVAL    *pgvBegCourse);

BOOL GK_API gkdCalcDistAndBrg (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgptBegin, 
                               SGeoPoint *pgptEnd,
                               GEOVAL    *pgvDistance,
                               GEOVAL    *pgvBegCourse,
                               GEOVAL    *pgvEndCourse,
                               BOOL       bUseLongestArc = FALSE);
BOOL GK_API gkdCalcDistAndBrg (int        nDatum,
                               ELineType  ltLineType,
                               GEOVAL     gvBegLat,
                               GEOVAL     gvBegLon,
                               GEOVAL     gvEndLat,
                               GEOVAL     gvEndLon,
                               GEOVAL    *pgvDistance,
                               GEOVAL    *pgvBegCourse,
                               GEOVAL    *pgvEndCourse,
                               BOOL       bUseLongestArc = FALSE);
BOOL GK_API gkdCalcPos        (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgptBegin, 
                               GEOVAL     gvDistance,
                               GEOVAL     gvBegCourse,
                               SGeoPoint *pgptEnd,
                               GEOVAL    *pgvEndCourse);
BOOL GK_API gkdCalcPos        (int        nDatum,
                               ELineType  ltLineType,
                               GEOVAL     gvBegLat,
                               GEOVAL     gvBegLon,
                               GEOVAL     gvDistance,
                               GEOVAL     gvBegCourse,
                               GEOVAL    *pgvEndLat,
                               GEOVAL    *pgvEndLon,
                               GEOVAL    *pgvEndCourse);
BOOL GK_API gkdCalcPosInv     (int        nDatum,
                               ELineType  ltLineType,
                               SGeoPoint *pgpBegin, 
                               GEOVAL     gvDistance,
                               GEOVAL     gvEndCourse,
                               SGeoPoint *pgpEnd,
                               GEOVAL    *pgvBegCourse);
#endif // _UNIVERSAL_DEFINED_

#ifndef _RHUMBLINE_FUNC_DEFINED_
#define _RHUMBLINE_FUNC_DEFINED_

// Rhumb line direct/backward calculations
BOOL GK_API gkrCalcRhumblineDistAndBrg (int        nDatumID,
                                        GEOVAL     gvBegLat, 
                                        GEOVAL     gvBegLong,
                                        GEOVAL     gvEndLat,
                                        GEOVAL     gvEndLong,
                                        GEOVAL    *pgvDistance,
                                        GEOVAL    *pgvCourse = NULL);
BOOL GK_API gkrCalcRhumblinePos        (int        nDatumID,
                                        GEOVAL     gvBegLat, 
                                        GEOVAL     gvBegLong,
                                        GEOVAL     gvDistance,
                                        GEOVAL     gvCourse,
                                        GEOVAL    *pgvEndLat,
                                        GEOVAL    *pgvEndLong);

BOOL GK_API gkdCalcRhumblineDistAndBrg (int        nDatumID,
                                        GEOVAL     gvBegLat, 
                                        GEOVAL     gvBegLong,
                                        GEOVAL     gvEndLat,
                                        GEOVAL     gvEndLong,
                                        GEOVAL    *pgvDistance,
                                        GEOVAL    *pgvCourse = NULL);
BOOL GK_API gkdCalcRhumblinePos        (int        nDatumID,
                                        GEOVAL     gvBegLat, 
                                        GEOVAL     gvBegLong,
                                        GEOVAL     gvDistance,
                                        GEOVAL     gvCourse,
                                        GEOVAL    *pgvEndLat,
                                        GEOVAL    *pgvEndLong);

// Precisious rhumb-line calculations
BOOL GK_API gkrCalcRhumblineDistAndBrg2 (int        nDatumID,
                                         GEOVAL     gvBegLat, 
                                         GEOVAL     gvBegLong,
                                         GEOVAL     gvEndLat,
                                         GEOVAL     gvEndLong,
                                         GEOVAL    *pgvDistance,
                                         GEOVAL    *pgvCourse = NULL);
BOOL GK_API gkrCalcRhumblinePos2        (int        nDatumID,
                                         GEOVAL     gvBegLat, 
                                         GEOVAL     gvBegLong,
                                         GEOVAL     gvDistance,
                                         GEOVAL     gvCourse,
                                         GEOVAL    *pgvEndLat,
                                         GEOVAL    *pgvEndLong);

BOOL GK_API gkdCalcRhumblineDistAndBrg2 (int        nDatumID,
                                         GEOVAL     gvBegLat, 
                                         GEOVAL     gvBegLong,
                                         GEOVAL     gvEndLat,
                                         GEOVAL     gvEndLong,
                                         GEOVAL    *pgvDistance,
                                         GEOVAL    *pgvCourse = NULL);
BOOL GK_API gkdCalcRhumblinePos2        (int        nDatumID,
                                         GEOVAL     gvBegLat, 
                                         GEOVAL     gvBegLong,
                                         GEOVAL     gvDistance,
                                         GEOVAL     gvCourse,
                                         GEOVAL    *pgvEndLat,
                                         GEOVAL    *pgvEndLong);

__forceinline BOOL gkrCalcRhumblinePos2 (int        nDatumID,
                                         SGeoPoint *pgpBegin, 
                                         GEOVAL     gvDistance,
                                         GEOVAL     gvCourse,
                                         SGeoPoint *pgpEnd)
{
    return gkrCalcRhumblinePos2 (nDatumID, pgpBegin->lat, pgpBegin->lon, gvDistance, gvCourse, & pgpEnd->lat, & pgpEnd->lon);
}

__forceinline BOOL gkrCalcRhumblineDistAndBrg2 (int        nDatumID,
                                                SGeoPoint *pgptBegin, 
                                                SGeoPoint *pgptEnd,
                                                GEOVAL    *pgvDistance,
                                                GEOVAL    *pgvCourse)
{
    return gkrCalcRhumblineDistAndBrg2 (nDatumID,
                                        pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                        pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                        pgvDistance, pgvCourse);
}

__forceinline BOOL gkdCalcRhumblinePos2 (int        nDatumID,
                                         SGeoPoint *pgpBegin, 
                                         GEOVAL     gvDistance,
                                         GEOVAL     gvCourse,
                                         SGeoPoint *pgpEnd)
{
    return gkdCalcRhumblinePos2 (nDatumID, pgpBegin->lat, pgpBegin->lon, gvDistance, gvCourse, & pgpEnd->lat, & pgpEnd->lon);
}

__forceinline BOOL gkdCalcRhumblineDistAndBrg2 (int        nDatumID,
                                                SGeoPoint *pgptBegin, 
                                                SGeoPoint *pgptEnd,
                                                GEOVAL    *pgvDistance,
                                                GEOVAL    *pgvCourse)
{
    return gkdCalcRhumblineDistAndBrg2 (nDatumID,
                                        pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                        pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                        pgvDistance, pgvCourse);
}

/////////////
__forceinline BOOL gkrCalcRhumblinePos (int        nDatumID,
                                        SGeoPoint *pgpBegin, 
                                        GEOVAL     gvDistance,
                                        GEOVAL     gvCourse,
                                        SGeoPoint *pgpEnd)
{
    return gkrCalcRhumblinePos (nDatumID, pgpBegin->lat, pgpBegin->lon, gvDistance, gvCourse, & pgpEnd->lat, & pgpEnd->lon);
}

__forceinline BOOL gkrCalcRhumblineDistAndBrg (int        nDatumID,
                                               SGeoPoint *pgptBegin, 
                                               SGeoPoint *pgptEnd,
                                               GEOVAL    *pgvDistance,
                                               GEOVAL    *pgvCourse)
{
    return gkrCalcRhumblineDistAndBrg (nDatumID,
                                       pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                       pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                       pgvDistance, pgvCourse);
}

__forceinline BOOL gkdCalcRhumblinePos (int        nDatumID,
                                        SGeoPoint *pgpBegin, 
                                        GEOVAL     gvDistance,
                                        GEOVAL     gvCourse,
                                        SGeoPoint *pgpEnd)
{
    return gkdCalcRhumblinePos (nDatumID, pgpBegin->lat, pgpBegin->lon, gvDistance, gvCourse, & pgpEnd->lat, & pgpEnd->lon);
}

__forceinline BOOL gkdCalcRhumblineDistAndBrg (int        nDatumID,
                                               SGeoPoint *pgptBegin, 
                                               SGeoPoint *pgptEnd,
                                               GEOVAL    *pgvDistance,
                                               GEOVAL    *pgvCourse)
{
    if (pgvDistance) *pgvDistance = 0.0;
    if (pgvCourse) *pgvCourse = 0.0;
    return gkdCalcRhumblineDistAndBrg (nDatumID,
                                       pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                       pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                       pgvDistance, pgvCourse);
}

__forceinline BOOL gkCalcRhumblineDistAndBrg (int        nDatumID,
                                               GEOVAL     gvBegLat, 
                                               GEOVAL     gvBegLong,
                                               GEOVAL     gvEndLat,
                                               GEOVAL     gvEndLong,
                                               GEOVAL    *pgvDistance,
                                               GEOVAL    *pgvCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcRhumblineDistAndBrg (nDatumID, gvBegLat, gvBegLong, gvEndLat, gvEndLong, pgvDistance, pgvCourse);
    else
        return gkrCalcRhumblineDistAndBrg (nDatumID, gvBegLat, gvBegLong, gvEndLat, gvEndLong, pgvDistance, pgvCourse);
}

__forceinline BOOL gkCalcRhumblineDistAndBrg (int        nDatumID,
                                              SGeoPoint *pgptBegin, 
                                              SGeoPoint *pgptEnd,
                                              GEOVAL    *pgvDistance,
                                              GEOVAL    *pgvCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcRhumblineDistAndBrg (nDatumID, pgptBegin, pgptEnd, pgvDistance, pgvCourse);
    else
        return gkrCalcRhumblineDistAndBrg (nDatumID, pgptBegin, pgptEnd, pgvDistance, pgvCourse);
}

__forceinline BOOL gkCalcRhumblinePos (int        nDatumID,
                                       GEOVAL     gvBegLat, 
                                       GEOVAL     gvBegLong,
                                       GEOVAL     gvDistance,
                                       GEOVAL     gvCourse,
                                       GEOVAL    *pgvEndLat,
                                       GEOVAL    *pgvEndLong)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcRhumblinePos (nDatumID, gvBegLat, gvBegLong, gvDistance, gvCourse, pgvEndLat, pgvEndLong);
    else
        return gkrCalcRhumblinePos (nDatumID, gvBegLat, gvBegLong, gvDistance, gvCourse, pgvEndLat, pgvEndLong);
}

__forceinline BOOL gkCalcRhumblinePos (int        nDatumID,
                                       SGeoPoint *pgptBegin, 
                                       GEOVAL     gvDistance,
                                       GEOVAL     gvCourse,
                                       SGeoPoint *pgptEnd)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcRhumblinePos (nDatumID, pgptBegin, gvDistance, gvCourse, pgptEnd);
    else
        return gkrCalcRhumblinePos (nDatumID, pgptBegin, gvDistance, gvCourse, pgptEnd);
}

#endif // _RHUMBLINE_FUNC_DEFINED_

#ifndef _GREATCIRCLE_FUNC_DEFINED_
#define _GREATCIRCLE_FUNC_DEFINED_

// Vincenty formulae realisation
BOOL GK_API gkrCalcGreatCircleDistAndBrg (int        nDatum,
                                          GEOVAL     gvBegLat, 
                                          GEOVAL     gvBegLong,
                                          GEOVAL     gvEndLat,
                                          GEOVAL     gvEndLong,
                                          GEOVAL    *pgvDistance,
                                          GEOVAL    *pgvBegCourse = NULL,
                                          GEOVAL    *pgvEndCourse = NULL);
BOOL GK_API gkrCalcGreatCirclePos        (int        nDatum,
                                          GEOVAL     gvBegLat, 
                                          GEOVAL     gvBegLong,
                                          GEOVAL     gvDistance,
                                          GEOVAL     gvBegCourse,
                                          GEOVAL    *pgvEndLat,
                                          GEOVAL    *gvEndLong,
                                          GEOVAL    *pgvEndCourse = NULL);

BOOL GK_API gkdCalcGreatCircleDistAndBrg (int        nDatum,
                                          GEOVAL     gvBegLat, 
                                          GEOVAL     gvBegLong,
                                          GEOVAL     gvEndLat,
                                          GEOVAL     gvEndLong,
                                          GEOVAL    *pgvDistance,
                                          GEOVAL    *pgvBegCourse = NULL,
                                          GEOVAL    *pgvEndCourse = NULL);
BOOL GK_API gkdCalcGreatCirclePos        (int        nDatum,
                                          GEOVAL     gvBegLat, 
                                          GEOVAL     gvBegLong,
                                          GEOVAL     gvDistance,
                                          GEOVAL     gvBegCourse,
                                          GEOVAL    *pgvEndLat,
                                          GEOVAL    *gvEndLong,
                                          GEOVAL    *pgvEndCourse = NULL);

// Sodano formulae realisation
BOOL GK_API gkrCalcGreatCircleDistAndBrg2 (int        nDatum,
                                           GEOVAL     gvBegLat, 
                                           GEOVAL     gvBegLong,
                                           GEOVAL     gvEndLat,
                                           GEOVAL     gvEndLong,
                                           GEOVAL    *pgvDistance,
                                           GEOVAL    *pgvBegCourse = NULL,
                                           GEOVAL    *pgvEndCourse = NULL);
BOOL GK_API gkrCalcGreatCirclePos2        (int        nDatum,
                                           GEOVAL     gvBegLat, 
                                           GEOVAL     gvBegLong,
                                           GEOVAL     gvDistance,
                                           GEOVAL     gvBegCourse,
                                           GEOVAL    *pgvEndLat,
                                           GEOVAL    *gvEndLong,
                                           GEOVAL    *pgvEndCourse = NULL);

BOOL GK_API gkdCalcGreatCircleDistAndBrg2 (int        nDatum,
                                           GEOVAL     gvBegLat, 
                                           GEOVAL     gvBegLong,
                                           GEOVAL     gvEndLat,
                                           GEOVAL     gvEndLong,
                                           GEOVAL    *pgvDistance,
                                           GEOVAL    *pgvBegCourse = NULL,
                                           GEOVAL    *pgvEndCourse = NULL);
BOOL GK_API gkdCalcGreatCirclePos2        (int        nDatum,
                                           GEOVAL     gvBegLat, 
                                           GEOVAL     gvBegLong,
                                           GEOVAL     gvDistance,
                                           GEOVAL     gvBegCourse,
                                           GEOVAL    *pgvEndLat,
                                           GEOVAL    *gvEndLong,
                                           GEOVAL    *pgvEndCourse = NULL);

__forceinline BOOL gkCalcGreatCircleDistAndBrg (int     nDatum,
                                                GEOVAL  gvBegLat, 
                                                GEOVAL  gvBegLong,
                                                GEOVAL  gvEndLat,
                                                GEOVAL  gvEndLong,
                                                GEOVAL *pgvDistance,
                                                GEOVAL *pgvBegCourse = NULL,
                                                GEOVAL *pgvEndCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcGreatCircleDistAndBrg (nDatum, gvBegLat, gvBegLong, gvEndLat, gvEndLong, pgvDistance, 
                                             pgvBegCourse, pgvEndCourse);
    else
        return gkrCalcGreatCircleDistAndBrg (nDatum, gvBegLat, gvBegLong, gvEndLat, gvEndLong, pgvDistance, 
                                             pgvBegCourse, pgvEndCourse);
}

// Calculate great circle pgvDistance and bearing by two geographical points
__forceinline BOOL gkrCalcGreatCircleDistAndBrg (int        nDatum,
                                                 SGeoPoint *pgptBegin, 
                                                 SGeoPoint *pgptEnd,
                                                 GEOVAL    *pgvDistance,
                                                 GEOVAL    *pgvBegCourse,
                                                 GEOVAL    *pgvEndCourse)
{
    return gkrCalcGreatCircleDistAndBrg (nDatum, 
                                         pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                         pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                         pgvDistance, pgvBegCourse, pgvEndCourse);
}

__forceinline BOOL gkdCalcGreatCircleDistAndBrg (int        nDatum,
                                                 SGeoPoint *pgptBegin, 
                                                 SGeoPoint *pgptEnd,
                                                 GEOVAL    *pgvDistance,
                                                 GEOVAL    *pgvBegCourse,
                                                 GEOVAL    *pgvEndCourse)
{
    if (pgvDistance) *pgvDistance = 0.0;
    if (pgvBegCourse) *pgvBegCourse = 0.0;
    if (pgvEndCourse) *pgvEndCourse = 0.0;
    return gkdCalcGreatCircleDistAndBrg (nDatum, 
                                         pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                         pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                         pgvDistance, pgvBegCourse, pgvEndCourse);
}

// Calculate great circle end position by begin point, pgvDistance and bearing
__forceinline BOOL gkrCalcGreatCirclePos (int        nDatum,
                                          SGeoPoint *pgptBegin, 
                                          GEOVAL     gvpgvDistance,
                                          GEOVAL     gvBegCourse,
                                          SGeoPoint *pgptEnd,
                                          GEOVAL    *pgvEndCourse)
{
    
    return gkrCalcGreatCirclePos (nDatum,
                                  pgptBegin->gvLatitude, pgptBegin->gvLongitude, gvpgvDistance, gvBegCourse,
                                  & (pgptEnd->gvLatitude), & (pgptEnd->gvLongitude), pgvEndCourse);
}

__forceinline BOOL gkdCalcGreatCirclePos (int        nDatum,
                                          SGeoPoint *pgptBegin, 
                                          GEOVAL     gvpgvDistance,
                                          GEOVAL     gvBegCourse,
                                          SGeoPoint *pgptEnd,
                                          GEOVAL    *pgvEndCourse)
{
    
    return gkdCalcGreatCirclePos (nDatum,
                                  pgptBegin->gvLatitude, pgptBegin->gvLongitude, gvpgvDistance, gvBegCourse,
                                  & (pgptEnd->gvLatitude), & (pgptEnd->gvLongitude), pgvEndCourse);
}
/////////////////////
// Calculate great circle pgvDistance and bearing by two geographical points
__forceinline BOOL gkrCalcGreatCircleDistAndBrg2 (int        nDatum,
                                                  SGeoPoint *pgptBegin, 
                                                  SGeoPoint *pgptEnd,
                                                  GEOVAL    *pgvDistance,
                                                  GEOVAL    *pgvBegCourse,
                                                  GEOVAL    *pgvEndCourse)
{
    return gkrCalcGreatCircleDistAndBrg2 (nDatum, 
                                          pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                          pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                          pgvDistance, pgvBegCourse, pgvEndCourse);
}

__forceinline BOOL gkdCalcGreatCircleDistAndBrg2 (int        nDatum,
                                                  SGeoPoint *pgptBegin, 
                                                  SGeoPoint *pgptEnd,
                                                  GEOVAL    *pgvDistance,
                                                  GEOVAL    *pgvBegCourse,
                                                  GEOVAL    *pgvEndCourse)
{
    return gkdCalcGreatCircleDistAndBrg2 (nDatum, 
                                          pgptBegin->gvLatitude, pgptBegin->gvLongitude,
                                          pgptEnd->gvLatitude, pgptEnd->gvLongitude,
                                          pgvDistance, pgvBegCourse, pgvEndCourse);
}

// Calculate great circle end position by begin point, pgvDistance and bearing
__forceinline BOOL gkrCalcGreatCirclePos2 (int        nDatum,
                                           SGeoPoint *pgptBegin, 
                                           GEOVAL     gvpgvDistance,
                                           GEOVAL     gvBegCourse,
                                           SGeoPoint *pgptEnd,
                                           GEOVAL    *pgvEndCourse)
{
    
    return gkrCalcGreatCirclePos2 (nDatum,
                                   pgptBegin->gvLatitude, pgptBegin->gvLongitude, gvpgvDistance, gvBegCourse,
                                   & (pgptEnd->gvLatitude), & (pgptEnd->gvLongitude), pgvEndCourse);
}

__forceinline BOOL gkdCalcGreatCirclePos2 (int        nDatum,
                                           SGeoPoint *pgptBegin, 
                                           GEOVAL     gvpgvDistance,
                                           GEOVAL     gvBegCourse,
                                           SGeoPoint *pgptEnd,
                                           GEOVAL    *pgvEndCourse)
{
    
    return gkdCalcGreatCirclePos2 (nDatum,
                                   pgptBegin->gvLatitude, pgptBegin->gvLongitude, gvpgvDistance, gvBegCourse,
                                   & (pgptEnd->gvLatitude), & (pgptEnd->gvLongitude), pgvEndCourse);
}
////////////////////
__forceinline BOOL gkCalcGreatCircleDistAndBrg (int        nDatum,
                                                SGeoPoint *pgptBegin, 
                                                SGeoPoint *pgptEnd,
                                                GEOVAL    *pgvDistance,
                                                GEOVAL    *pgvBegCourse = NULL,
                                                GEOVAL    *pgvEndCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcGreatCircleDistAndBrg (nDatum, pgptBegin, pgptEnd, pgvDistance, pgvBegCourse, pgvEndCourse);
    else
        return gkrCalcGreatCircleDistAndBrg (nDatum, pgptBegin, pgptEnd, pgvDistance, pgvBegCourse, pgvEndCourse);
}

__forceinline BOOL gkCalcGreatCirclePos (int     nDatum,
                                         GEOVAL  gvBegLat, 
                                         GEOVAL  gvBegLong,
                                         GEOVAL  gvDistance,
                                         GEOVAL  gvBegCourse,
                                         GEOVAL *pgvEndLat,
                                         GEOVAL *gvEndLong,
                                         GEOVAL *pgvEndCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcGreatCirclePos (nDatum, gvBegLat, gvBegLong, gvDistance, gvBegCourse, pgvEndLat, gvEndLong, pgvEndCourse);
    else
        return gkrCalcGreatCirclePos (nDatum, gvBegLat, gvBegLong, gvDistance, gvBegCourse, pgvEndLat, gvEndLong, pgvEndCourse);
}

__forceinline BOOL gkCalcGreatCirclePos (int        nDatum,
                                         SGeoPoint *pgptBegin, 
                                         GEOVAL     gvDistance,
                                         GEOVAL     gvBegCourse,
                                         SGeoPoint *pgptEnd,
                                         GEOVAL    *pgvEndCourse = NULL)
{
    if (gkIsDegreeAssumed ())
        return gkdCalcGreatCirclePos (nDatum, pgptBegin, gvDistance, gvBegCourse, pgptEnd, pgvEndCourse);
    else
        return gkrCalcGreatCirclePos (nDatum, pgptBegin, gvDistance, gvBegCourse, pgptEnd, pgvEndCourse);
}


#endif // _GREATCIRCLE_FUNC_DEFINED_

#ifndef _TOPOLOGY_FUNC_DEFINED_
#define _TOPOLOGY_FUNC_DEFINED_

// Are two lines parallel?
BOOL GK_API gkAreParallel (UGeoSegment *pgsgmLine1, UGeoSegment *pgsgmLine2);
BOOL GK_API gkAreParallel (SGeoPoint   *pgptBegPoint1, 
                           SGeoPoint   *pgptEndPoint1,
                           SGeoPoint   *pgptBegPoint2, 
                           SGeoPoint   *pgptEndPoint2);
BOOL GK_API gkAreParallel (GEOVAL       gvBegLat1,
                           GEOVAL       gvBegLong1,
                           GEOVAL       gvEndLat1,
                           GEOVAL       gvEndLong1,
                           GEOVAL       gvBegLat2,
                           GEOVAL       gvBegLong2,
                           GEOVAL       gvEndLat2,
                           GEOVAL       gvEndLong2);

// FOR RADIANS:
// Crossing line and another line
BOOL GK_API gkrCrossLineByLine (UGeoSegment *pgsgmLine1, 
                                UGeoSegment *pgsgmLine2, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossLineByLine (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossLineByLine (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);
// Crossing line and segment
BOOL GK_API gkrCrossLineBySegm (UGeoSegment *pgsgmLine, 
                                UGeoSegment *pgsgmSegm, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossLineBySegm (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossLineBySegm (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Crossing segment and line
BOOL GK_API gkrCrossSegmByLine (UGeoSegment *pgsgmSegm, 
                                UGeoSegment *pgsgmLine, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossSegmByLine (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossSegmByLine (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Crossing two different segments
BOOL GK_API gkrCrossSegmBySegm (UGeoSegment *pgsgmSegm1, 
                                UGeoSegment *pgsgmSegm2, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossSegmBySegm (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkrCrossSegmBySegm (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Get relation of the point and the line
BOOL GK_API gkrIsPointOnLine (GEOVAL       gvPointLat,
                              GEOVAL       gvPointLong,
                              GEOVAL       gvBegLat,
                              GEOVAL       gvBegLong,
                              GEOVAL       gvEndLat,
                              GEOVAL       gvEndLong,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkrIsPointOnLine (SGeoPoint   *pgptPoint,
                              SGeoPoint   *pgptBegPoint,
                              SGeoPoint   *pgptEndPoint,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkrIsPointOnLine (SGeoPoint   *pgptPoint,
                              UGeoSegment *pgsLine,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);

// Get relation of the point and the segment
BOOL GK_API gkrIsPointOnSegm (GEOVAL       gvPointLat,
                              GEOVAL       gvPointLong,
                              GEOVAL       gvBegLat,
                              GEOVAL       gvBegLong,
                              GEOVAL       gvEndLat,
                              GEOVAL       gvEndLong,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkrIsPointOnSegm (SGeoPoint   *pgptPoint,
                              SGeoPoint   *pgptBegPoint,
                              SGeoPoint   *pgptEndPoint,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkrIsPointOnSegm (SGeoPoint   *pgptPoint,
                              UGeoSegment *pgsSegm,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
 
// Get relation of the point to the contour
EPointPrimRelation GK_API gkrGetPointAndContourRel (SGeoPoint            *pgptPoint, 
                                                    SGeoContour          *pgcContour,
                                                    BOOL                  bThroughDCL = FALSE,
                                                    BOOL                  bExactCheck = TRUE,
                                                    EPointEntryTestMethod petmMethod  = petmByGCAngles,
                                                    int                   nDatum      = GD_GEOID_WGS84);
EPointPrimRelation GK_API gkrGetPointAndContourRel (SGeoPoint            *pgptPoint, 
                                                    SGeoPoint            *pgptVertexes, 
                                                    int                   nVertexCount,
                                                    BOOL                  bThroughDCL = FALSE,
                                                    BOOL                  bExactCheck = TRUE,
                                                    EPointEntryTestMethod petmMethod  = petmByGCAngles,
                                                    int                   nDatum      = GD_GEOID_WGS84);
 
EContourContourRelation GK_API gkrGetContourAndContourRel (SGeoPoint   *pgptVertexes1,
                                                           int          nVertexCount1,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoPoint   *pgptVertexes2,
                                                           int          nVertexCount2,
                                                           GEOVAL       gvPhaseShift2,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);
EContourContourRelation GK_API gkrGetContourAndContourRel (SGeoContour *pgcContour1,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoContour *pgcContour2,
                                                           GEOVAL       gvPhaseShift2,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);

BOOL GK_API gkrCrossSegmByContour (SGeoPoint *pgptPoint1,
                                   SGeoPoint *pgptPoint2,
                                   SGeoPoint *pgptVertexes,
                                   int        nVertexCount,
                                   int       *pnCrossCount);
BOOL GK_API gkrCrossSegmByPolygon (SGeoPoint *pgptPoint1,
                                   SGeoPoint *pgptPoint2,
                                   SGeoPoint *pgptVertexes,
                                   int        nVertexCount,
                                   int       *pnCrossCount);

BOOL GK_API gkrCheckSegmCrossing (int          nDatum,
                                  GEOVAL       gvBegLat1,
                                  GEOVAL       gvBegLon1,
                                  GEOVAL       gvBearing1,
                                  GEOVAL       gvDistance1,
                                  GEOVAL       gvBegLat2,
                                  GEOVAL       gvBegLon2,
                                  GEOVAL       gvBearing2,
                                  GEOVAL       gvDistance2,
                                  ELineType    ltLineType);
BOOL GK_API gkrCheckSegmCrossing (int          nDatum,
                                  SGeoPoint   *pgpBegPoint1,
                                  GEOVAL       gvBearing1,
                                  GEOVAL       gvDistance1,
                                  SGeoPoint   *pgpBegPoint2,
                                  GEOVAL       gvBearing2,
                                  GEOVAL       gvDistance2,
                                  ELineType    ltLineType);
BOOL GK_API gkrCheckSegmCrossing (int          nDatum,
                                  SGeoSegment *pgsSegemnt1,
                                  SGeoSegment *pgsSegemnt2,
                                  ELineType    ltLineType);

BOOL GK_API gkrCheckSegmCrossing2 (int       nDatum,
                                   GEOVAL    gvBegLat1,        // P1
                                   GEOVAL    gvBegLon1,
                                   GEOVAL    gvEndLat1,        // P1
                                   GEOVAL    gvEndLon1,
                                   GEOVAL    gvBegLat2,        // P3
                                   GEOVAL    gvBegLon2,
                                   GEOVAL    gvEndLat2,        // P4
                                   GEOVAL    gvEndLon2,
                                   ELineType ltLineType,
                                   BOOL      bUseLongestArc1 = FALSE,
                                   BOOL      bUseLongestArc2 = FALSE);

// FOR DEGREES:
// Crossing line and another line
BOOL GK_API gkdCrossLineByLine (UGeoSegment *pgsgmLine1, 
                                UGeoSegment *pgsgmLine2, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossLineByLine (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossLineByLine (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);
// Crossing line and segment
BOOL GK_API gkdCrossLineBySegm (UGeoSegment *pgsgmLine, 
                                UGeoSegment *pgsgmSegm, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossLineBySegm (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossLineBySegm (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Crossing segment and line
BOOL GK_API gkdCrossSegmByLine (UGeoSegment *pgsgmSegm, 
                                UGeoSegment *pgsgmLine, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossSegmByLine (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossSegmByLine (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Crossing two different segments
BOOL GK_API gkdCrossSegmBySegm (UGeoSegment *pgsgmSegm1, 
                                UGeoSegment *pgsgmSegm2, 
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossSegmBySegm (SGeoPoint   *pgptBegPoint1, 
                                SGeoPoint   *pgptEndPoint1,
                                SGeoPoint   *pgptBegPoint2, 
                                SGeoPoint   *pgptEndPoint2,
                                SGeoPoint   *pgptCrossPoint);
BOOL GK_API gkdCrossSegmBySegm (GEOVAL       gvBegLat1,
                                GEOVAL       gvBegLong1,
                                GEOVAL       gvEndLat1,
                                GEOVAL       gvEndLong1,
                                GEOVAL       gvBegLat2,
                                GEOVAL       gvBegLong2,
                                GEOVAL       gvEndLat2,
                                GEOVAL       gvEndLong2,
                                GEOVAL      *pgvCrossPointLat,
                                GEOVAL      *pgvCrossPointLong);

// Get relation of the point and the line
BOOL GK_API gkdIsPointOnLine (GEOVAL       gvPointLat,
                              GEOVAL       gvPointLong,
                              GEOVAL       gvBegLat,
                              GEOVAL       gvBegLong,
                              GEOVAL       gvEndLat,
                              GEOVAL       gvEndLong,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkdIsPointOnLine (SGeoPoint   *pgptPoint,
                              SGeoPoint   *pgptBegPoint,
                              SGeoPoint   *pgptEndPoint,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkdIsPointOnLine (SGeoPoint   *pgptPoint,
                              UGeoSegment *pgsLine,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);

// Get relation of the point and the segment
BOOL GK_API gkdIsPointOnSegm (GEOVAL       gvPointLat,
                              GEOVAL       gvPointLong,
                              GEOVAL       gvBegLat,
                              GEOVAL       gvBegLong,
                              GEOVAL       gvEndLat,
                              GEOVAL       gvEndLong,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkdIsPointOnSegm (SGeoPoint   *pgptPoint,
                              SGeoPoint   *pgptBegPoint,
                              SGeoPoint   *pgptEndPoint,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
BOOL GK_API gkdIsPointOnSegm (SGeoPoint   *pgptPoint,
                              UGeoSegment *pgsSegm,
                              BOOL         bExactCheck = TRUE,
                              ELineType    ltLineType  = ltGreatCircle,
                              int          nDatum      = GD_GEOID_WGS84);
 
// Get relation of the point to the contour
EPointPrimRelation GK_API gkdGetPointAndContourRel (SGeoPoint            *pgptPoint, 
                                                    SGeoContour          *pgcContour,
                                                    BOOL                  bThroughDCL = FALSE,
                                                    BOOL                  bExactCheck = TRUE,
                                                    EPointEntryTestMethod petmMethod  = petmByGCAngles,
                                                    int                   nDatum      = GD_GEOID_WGS84);
EPointPrimRelation GK_API gkdGetPointAndContourRel (SGeoPoint            *pgptPoint, 
                                                    SGeoPoint            *pgptVertexes, 
                                                    int                   nVertexCount,
                                                    BOOL                  bThroughDCL = FALSE,
                                                    BOOL                  bExactCheck = TRUE,
                                                    EPointEntryTestMethod petmMethod  = petmByGCAngles,
                                                    int                   nDatum      = GD_GEOID_WGS84);
 
EContourContourRelation GK_API gkdGetContourAndContourRel (SGeoPoint   *pgptVertexes1,
                                                           int          nVertexCount1,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoPoint   *pgptVertexes2,
                                                           int          nVertexCount2,
                                                           GEOVAL       gvPhaseShift2,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);
EContourContourRelation GK_API gkdGetContourAndContourRel (SGeoContour *pgcContour1,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoContour *pgcContour2,
                                                           GEOVAL       gvPhaseShift2,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);

EContourContourRelation GK_API gkdGetSegmentAndContourRel (SGeoPoint *pgpBegin,
                                                           SGeoPoint *pgpEnd,
                                                           GEOVAL     gvPhaseShift1,
                                                           SGeoPoint *pgptVertexes,
                                                           int        nVertexCount,
                                                           GEOVAL     gvPhaseShift2,
                                                           ELineType  ltLineType = ltGreatCircle,
                                                           int        nDatum     = GD_GEOID_WGS84);
EContourContourRelation GK_API gkdGetSegmentAndContourRel (SGeoPoint   *pgpBegin,
                                                           SGeoPoint   *pgpEnd,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoContour *pgcContour,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);
EContourContourRelation GK_API gkrGetSegmentAndContourRel (SGeoPoint *pgpBegin,
                                                           SGeoPoint *pgpEnd,
                                                           GEOVAL     gvPhaseShift1,
                                                           SGeoPoint *pgptVertexes,
                                                           int        nVertexCount,
                                                           GEOVAL     gvPhaseShift2,
                                                           ELineType  ltLineType = ltGreatCircle,
                                                           int        nDatum     = GD_GEOID_WGS84);
EContourContourRelation GK_API gkrGetSegmentAndContourRel (SGeoPoint   *pgpBegin,
                                                           SGeoPoint   *pgpEnd,
                                                           GEOVAL       gvPhaseShift1,
                                                           SGeoContour *pgcContour,
                                                           ELineType    ltLineType = ltGreatCircle,
                                                           int          nDatum     = GD_GEOID_WGS84);

BOOL GK_API gkdCrossSegmByContour (SGeoPoint *pgptPoint1,
                                   SGeoPoint *pgptPoint2,
                                   SGeoPoint *pgptVertexes,
                                   int        nVertexCount,
                                   int       *pnCrossCount);
BOOL GK_API gkdCrossSegmByPolygon (SGeoPoint *pgptPoint1,
                                   SGeoPoint *pgptPoint2,
                                   SGeoPoint *pgptVertexes,
                                   int        nVertexCount,
                                   int       *pnCrossCount);

BOOL GK_API gkdCheckSegmCrossing (int          nDatum,
                                  GEOVAL       gvBegLat1,
                                  GEOVAL       gvBegLon1,
                                  GEOVAL       gvBearing1,
                                  GEOVAL       gvDistance1,
                                  GEOVAL       gvBegLat2,
                                  GEOVAL       gvBegLon2,
                                  GEOVAL       gvBearing2,
                                  GEOVAL       gvDistance2,
                                  ELineType    ltLineType);
BOOL GK_API gkdCheckSegmCrossing (int          nDatum,
                                  SGeoPoint   *pgpBegPoint1,
                                  GEOVAL       gvBearing1,
                                  GEOVAL       gvDistance1,
                                  SGeoPoint   *pgpBegPoint2,
                                  GEOVAL       gvBearing2,
                                  GEOVAL       gvDistance2,
                                  ELineType    ltLineType);
BOOL GK_API gkdCheckSegmCrossing (int          nDatum,
                                  SGeoSegment *pgsSegemnt1,
                                  SGeoSegment *pgsSegemnt2,
                                  ELineType    ltLineType);

BOOL GK_API gkdCheckSegmCrossing2 (int       nDatum,
                                   GEOVAL    gvBegLat1,        // P1
                                   GEOVAL    gvBegLon1,
                                   GEOVAL    gvEndLat1,        // P1
                                   GEOVAL    gvEndLon1,
                                   GEOVAL    gvBegLat2,        // P3
                                   GEOVAL    gvBegLon2,
                                   GEOVAL    gvEndLat2,        // P4
                                   GEOVAL    gvEndLon2,
                                   ELineType ltLineType,
                                   BOOL      bUseLongestArc1 = FALSE,
                                   BOOL      bUseLongestArc2 = FALSE);

BOOL GK_API gkdQuickAre2VectCross (AVector vctVector1, AVector vctVector2, int nGeoid = GD_GEOID_WGS84, ELineType ltLineType = ltGreatCircle);
BOOL GK_API gkrQuickAre2VectCross (AVector vctVector1, AVector vctVector2, int nGeoid = GD_GEOID_WGS84, ELineType ltLineType = ltGreatCircle);

BOOL GK_API gkdQuickIsVectCrossPoly (AVector    vctVector, 
                                     SGeoPoint *pgpVertexes,
                                     int        nVertexCount,
                                     int        nGeoid     = GD_GEOID_WGS84,
                                     ELineType  ltLineType = ltGreatCircle);
BOOL GK_API gkrQuickIsVectCrossPoly (AVector    vctVector, 
                                     SGeoPoint *pgpVertexes,
                                     int        nVertexCount,
                                     int        nGeoid     = GD_GEOID_WGS84,
                                     ELineType  ltLineType = ltGreatCircle);

BOOL GK_API gkdIsAreaCrossArea (SGeoContour *pgcContour1,
                                SGeoContour *pgcContour2,
                                int          nGeoid     = GD_GEOID_WGS84,
                                ELineType    ltLineType = ltGreatCircle);
BOOL GK_API gkrIsAreaCrossArea (SGeoContour *pgcContour1,
                                SGeoContour *pgcContour2,
                                int          nGeoid     = GD_GEOID_WGS84,
                                ELineType    ltLineType = ltGreatCircle);

#endif // _TOPOLOGY_FUNC_DEFINED_

#ifndef _BUILDERS_DEFINED_
#define _BUILDERS_DEFINED_

// Build polygon by start and end point and begin step (will be scaled depends to latitude)
BOOL GK_API gkrBuildGeoPolygon (int         nDatumID,
                                SGeoPoint  *pgpBegin, 
                                SGeoPoint  *pgpEnd, 
                                ELineType   ltLineType,
                                GEOVAL      gvBeginStep,
                                SGeoPoint **ppgpPolygon, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);
BOOL GK_API gkrBuildGeoArc     (int         nDatumID,
                                SGeoPoint  *pgpCenter, 
                                SGeoPoint  *pgpBegin, 
                                GEOVAL      gvAngle, 
                                ELineType   ltLineType,
                                GEOVAL      gvBeginStep,
                                SGeoPoint **ppgpArc, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);
BOOL GK_API gkrBuildGeoArc2    (int         nDatumID,
                                SGeoPoint  *pgpCenter, 
                                SGeoPoint  *pgpBegin, 
                                GEOVAL      gvAngle, 
                                ELineType   ltLineType,
                                GEOVAL      gvAngleStep,
                                SGeoPoint **ppgpArc, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);

BOOL GK_API gkdBuildGeoPolygon (int         nDatumID,
                                SGeoPoint  *pgpBegin, 
                                SGeoPoint  *pgpEnd, 
                                ELineType   ltLineType,
                                GEOVAL      gvBeginStep,
                                SGeoPoint **ppgpPolygon, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);
BOOL GK_API gkdBuildGeoArc     (int         nDatumID,
                                SGeoPoint  *pgpCenter, 
                                SGeoPoint  *pgpBegin, 
                                GEOVAL      gvAngle, 
                                ELineType   ltLineType,
                                GEOVAL      gvBeginStep,
                                SGeoPoint **ppgpArc, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);
BOOL GK_API gkdBuildGeoArc2    (int         nDatumID,
                                SGeoPoint  *pgpCenter, 
                                SGeoPoint  *pgpBegin, 
                                GEOVAL      gvAngle, 
                                ELineType   ltLineType,
                                GEOVAL      gvAngleStep,
                                SGeoPoint **ppgpArc, 
                                int        *pnPointNumber,
                                BOOL        bIncludeEndPoint = TRUE);

void GK_API gkReleaseMem (void *pPointer);

#endif // _BUILDERS_DEFINED_

#ifndef _ADDIT_GEO_DEFINED_
#define _ADDIT_GEO_DEFINED_

// Additional geometry functions
BOOL GK_API gkrDistFromPointToLine (int             nDatumID,
                                    SGeoPoint      *pgpPoint, 
                                    SGeoPoint      *pgpBegin, 
                                    SGeoPoint      *pgpEnd, 
                                    ELineType       ltLineType,
                                    EDistCalcMethod dcmMethod,
                                    GEOVAL         *pgvDistance,
                                    SGeoPoint      *pgpNearestPoint,
                                    BOOL            bOnlyBetween = FALSE);
BOOL GK_API gkrDistFromPointToArc  (int             nDatumID,
                                    SGeoPoint      *pgpPoint, 
                                    SGeoPoint      *pgpCenter, 
                                    SGeoPoint      *pgpBegin, 
                                    GEOVAL          gvAngle, 
                                    ELineType       ltLineType,
                                    GEOVAL         *pgvDistance,
                                    SGeoPoint      *pgpNearestPoint);
BOOL GK_API gkdDistFromPointToLine (int             nDatumID,
                                    SGeoPoint      *pgpPoint, 
                                    SGeoPoint      *pgpBegin, 
                                    SGeoPoint      *pgpEnd, 
                                    ELineType       ltLineType,
                                    EDistCalcMethod dcmMethod,
                                    GEOVAL         *pgvDistance,
                                    SGeoPoint      *pgpNearestPoint,
                                    BOOL            bOnlyBetween = FALSE);
BOOL GK_API gkdDistFromPointToArc  (int             nDatumID,
                                    SGeoPoint      *pgpPoint, 
                                    SGeoPoint      *pgpCenter, 
                                    SGeoPoint      *pgpBegin, 
                                    GEOVAL          gvAngle, 
                                    ELineType       ltLineType,
                                    GEOVAL         *pgvDistance,
                                    SGeoPoint      *pgpNearestPoint);
void GK_API gkCvPointToDeg         (SGeoPoint      *pgpPoint);
void GK_API gkCvPointToRad         (SGeoPoint      *pgpPoint);

#endif // _ADDIT_GEO_DEFINED_

#ifndef _GEOLINESCROSS_DEFINED_
#define _GEOLINESCROSS_DEFINED_

BOOL GK_API gkrCrossGeoSegments      (int          nDatumID,
                                      SGeoSegment *pgsSegment1,
                                      SGeoSegment *pgsSegment2,
                                      ELineType    ltLineType,
                                      SGeoPoint   *pgpCrossPoint  = NULL,
                                      BOOL         bSkipAllChecks = FALSE, // For internal use only! Does not set TRUE here!
                                      SGeoPoint   *pgpEnd1        = NULL,  // To speed up the calculation
                                      SGeoPoint   *pgpEnd2        = NULL);
BOOL GK_API gkrCrossGeoSegmentAndArc (int          nDatumID,
                                      SGeoSegment *pgsSegment,
                                      SGeoArc     *pgaArc,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkrCrossGeoSegmentAndArc (int          nDatumID,
                                      SGeoSegment *pgsSegment,
                                      SGeoPoint   *pgpCenter,
                                      SGeoPoint   *pgpBegin,
                                      GEOVAL       gvAngle,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkrCrossGeoArcs          (int          nDatumID,
                                      SGeoArc     *pgaArc1,
                                      SGeoArc     *pgaArc2,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkrCrossGeoArcs          (int          nDatumID,
                                      SGeoPoint   *pgpCenter1,
                                      SGeoPoint   *pgpBegin1,
                                      GEOVAL       gvAngle1,
                                      SGeoPoint   *pgpCenter2,
                                      SGeoPoint   *pgpBegin2,
                                      GEOVAL       gvAngle2,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);

BOOL GK_API gkdCrossGeoSegments      (int          nDatumID,
                                      SGeoSegment *pgsSegment1,
                                      SGeoSegment *pgsSegment2,
                                      ELineType    ltLineType,
                                      SGeoPoint   *pgpCrossPoint  = NULL,
                                      BOOL         bSkipAllChecks = FALSE, // For internal use only! Does not set TRUE here!
                                      SGeoPoint   *pgpEnd1        = NULL,  // To speed up the calculation
                                      SGeoPoint   *pgpEnd2        = NULL);
BOOL GK_API gkdCrossGeoSegmentAndArc (int          nDatumID,
                                      SGeoSegment *pgsSegment,
                                      SGeoArc     *pgaArc,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkdCrossGeoSegmentAndArc (int          nDatumID,
                                      SGeoSegment *pgsSegment,
                                      SGeoPoint   *pgpCenter,
                                      SGeoPoint   *pgpBegin,
                                      GEOVAL       gvAngle,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkdCrossGeoArcs          (int          nDatumID,
                                      SGeoArc     *pgaArc1,
                                      SGeoArc     *pgaArc2,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
BOOL GK_API gkdCrossGeoArcs          (int          nDatumID,
                                      SGeoPoint   *pgpCenter1,
                                      SGeoPoint   *pgpBegin1,
                                      GEOVAL       gvAngle1,
                                      SGeoPoint   *pgpCenter2,
                                      SGeoPoint   *pgpBegin2,
                                      GEOVAL       gvAngle2,
                                      ELineType    ltLineType,
                                      int         *pnCrossCount,
                                      SGeoPoint   *pgpCrossPoint1,
                                      SGeoPoint   *pgpCrossPoint2);
                               
#endif // _GEOLINESCROSS_DEFINED_

// Tracing
void GK_API gkEnableTrace (BOOL bEnable);

#ifndef _MARIS_GEOMETRY_REGIONS_
#define _MARIS_GEOMETRY_REGIONS_

typedef enum tagRegionOperation
{
    roAnd,
    roOr,
    roSubtract,
    roNegSubtract
}
ERegionOperation;

// Contour serice
SGeoContour GK_API *gkCreateContour (void);
void GK_API gkDeleteContour (SGeoContour *pgcContour);
void GK_API gkAppendPointToContour (SGeoContour *pgcContour, SGeoPoint *pgpPoint);

void GK_API gkTurnOverContour (SGeoContour *pgcContour);

BOOL GK_API gkrIsContourClockwise     (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);
__forceinline BOOL gkrIsContourUnclockwise   (SGeoContour *pgcContour, 
                                              ELineType    ltLineType    = ltGreatCircle, 
                                              BOOL         bQuickCheck   = FALSE,
                                              GEOVAL       gvEccentrcity = -1.0,
                                              SGeoid      *pdatDatum     = NULL) { return !gkrIsContourClockwise (pgcContour, ltLineType, bQuickCheck, gvEccentrcity, pdatDatum); }
void GK_API gkrMakeContourClockwise   (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);
void GK_API gkrMakeContourUnclockwise (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);

BOOL GK_API gkdIsContourClockwise     (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);
__forceinline BOOL gkdIsContourUnclockwise   (SGeoContour *pgcContour, 
                                              ELineType    ltLineType    = ltGreatCircle, 
                                              BOOL         bQuickCheck   = FALSE,
                                              GEOVAL       gvEccentrcity = -1.0,
                                              SGeoid      *pdatDatum     = NULL) { return !gkdIsContourClockwise (pgcContour, ltLineType, bQuickCheck, gvEccentrcity, pdatDatum); }
void GK_API gkdMakeContourClockwise   (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);
void GK_API gkdMakeContourUnclockwise (SGeoContour *pgcContour, 
                                       ELineType    ltLineType    = ltGreatCircle, 
                                       BOOL         bQuickCheck   = FALSE,
                                       GEOVAL       gvEccentrcity = -1.0,
                                       SGeoid      *pdatDatum     = NULL);

void GK_API gkSaveContour (HANDLE hFile, SGeoContour *pgcContour, int nIndentLevel = 0);
SGeoContour GK_API *gkLoadContour (HANDLE hFile, TCHAR *pszStoppingString = NULL);

// Region service
SGeoRegion GK_API *gkCreateRegion (void);
void GK_API gkDeleteRegion (SGeoRegion *pgrRegion);
SGeoRegion GK_API *gkCreateSingleBufferedRegion (int        nExtContourCount,
                                                 int        nIntContourCount,
                                                 int       *pnContoursSizeArray,
                                                 SGeoPoint *pgptVertexArray,
                                                 GEOVAL     gvPhaseShift = 0.0);
void GK_API gkFreeSingleBufferedRegion (SGeoRegion *pgrRegion);
void GK_API gkAddContourToRegion (SGeoRegion *pgrRegion, SGeoContour *pgcContour, BOOL bExternal);
void GK_API gkRemoveContourFromRegion (SGeoRegion *pgrRegion, int nContourNo, BOOL bExternal);

void GK_API gkSaveRegion (HANDLE hFile, SGeoRegion *pgrRegion, int nIndentLevel = 0);
SGeoRegion GK_API *gkLoadRegion (HANDLE hFile, TCHAR *pszStoppingString = NULL);

void GK_API gkSaveRegionList (const TCHAR *pszPath, SGeoRegion **ppgrRegions, int nRegionCount, int nIndentLevel = 0);
void GK_API gkLoadRegionList (const TCHAR *pszPath, SGeoRegion **ppgrRegions, int nMaxRegions, int *pnRegionCount, BOOL *pbUseRadians);

SGeoRegion  GK_API *gkrCombineTwoRegions (SGeoRegion      *pgrRegion1, 
                                          SGeoRegion      *pgrRegion2, 
                                          ERegionOperation roOperation,
                                          int             *pnResultContoursNumber,
                                          ELineType        ltLineType  = ltGreatCircle,
                                          int              nDatumID    = GD_GEOID_WGS84,
                                          BOOL             bQuickCheck = FALSE);
BOOL GK_API gkrQuickIsPointInsideContour (const SGeoPoint *pgpCheckPoint,
                                          const SGeoPoint *pgpVertexes,
                                          int              nVertexCount,
                                          GEOVAL           gvPhaseShift,
                                          ELineType        ltLineType  = ltRhumbLine,
                                          int              nDatum      = GD_GEOID_WGS84,
                                          BOOL             bQuickCheck = FALSE);
BOOL GK_API gkrQuickCheckCrossGeoSegments (int        nDatumID,
                                           SGeoPoint *pgpBegin1,
                                           SGeoPoint *pgpEnd1,
                                           SGeoPoint *pgpBegin2,
                                           SGeoPoint *pgpEnd2,
                                           ELineType  ltLineType,
                                           BOOL       bQuickCheck,
                                           SGeoPoint *pgpCrossPoint = NULL);

SGeoRegion  GK_API *gkdCombineTwoRegions (SGeoRegion      *pgrRegion1, 
                                          SGeoRegion      *pgrRegion2, 
                                          ERegionOperation roOperation,
                                          int             *pnResultContoursNumber,
                                          ELineType        ltLineType  = ltGreatCircle,
                                          int              nDatumID    = GD_GEOID_WGS84,
                                          BOOL             bQuickCheck = FALSE);
BOOL GK_API gkdQuickIsPointInsideContour (const SGeoPoint *pgpCheckPoint,
                                          const SGeoPoint *pgpVertexes,
                                          int              nVertexCount,
                                          GEOVAL           gvPhaseShift,
                                          ELineType        ltLineType  = ltRhumbLine,
                                          int              nDatum      = GD_GEOID_WGS84,
                                          BOOL             bQuickCheck = FALSE);
BOOL GK_API gkdQuickCheckCrossGeoSegments (int        nDatumID,
                                           SGeoPoint *pgpBegin1,
                                           SGeoPoint *pgpEnd1,
                                           SGeoPoint *pgpBegin2,
                                           SGeoPoint *pgpEnd2,
                                           ELineType  ltLineType,
                                           BOOL       bQuickCheck,
                                           SGeoPoint *pgpCrossPoint = NULL);

void GK_API gkContourRad2Deg (SGeoContour *pgcContour);
void GK_API gkContourDeg2Rad (SGeoContour *pgcContour);

void GK_API gkRegionDeg2Rad (SGeoRegion *pgrRegion);
void GK_API gkRegionRad2Deg (SGeoRegion *pgrRegion);

#endif // _MARIS_GEOMETRY_REGIONS_

#pragma pack()

#endif // _MARIS_GEOMETRY_INTERFACE_

