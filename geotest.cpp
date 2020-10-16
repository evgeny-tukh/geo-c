#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "geo.h"

#include <Windows.h>
#include "Library Interface.h"

void die (char *msg = 0, int code = 1) {
    if (msg) printf (msg);

    exit (code);
}

void showHelp () {
    printf (
        "\n\nUsage:\n\n"
        "gt command options\n\n"
        "commands are:\n"
        "\tb\tcalculate bearing and range from one point to another one\n"
        "\tp\tcalculate point by origin point, bearing and range\n"
        "\th|?\thelp\n\n"
        "options are:\n"
        "\t-o:lat,lon\torigin position\n"
        "\t-d:lat,lon\tdestination position\n"
        "\t-r:value\trange from origin (nm)\n"
        "\t-b:value\tbearing from origin (deg)\n"
        "\t-m:o[rthodromy]|l[oxodromy]\n\n"
    );

    die ();
}

bool parsePosition (geo::Pos& position, char *source) {
    char *comma = strchr (source, ',');
    bool result = false;

    if (comma) {
        position.lat = atof (source);
        position.lon = atof (comma + 1);

        result = position.lat >= -85.0 && position.lat <= 85.0 && position.lon >= -180.0 && position.lon <= 180.0;
    }

    return result;
}

int main (int argCount, char *args []) {
    geo::Operation operation;
    geo::Pos origin { 1.0e3, 1.0e3 }, dest { 1.0e3, 1.0e3 };
    double range = -1.0, bearing = -1.0;
    bool rhumbline = true;

    printf ("Neptune Geo test tool\n");

    if (argCount < 2) showHelp ();

    char operChar = tolower (args [1][0]);

    switch (operChar) {
        case geo::Operation::CALC_BRG_RNG:
        case geo::Operation::CALC_DEST_POINT:
            operation = (geo::Operation) operChar; break;
        default:
            die ("Invalid command");
    }

    for (auto i = 2; i < argCount; ++ i) {
        auto arg = args [i];

        if (arg [0] != '-' && arg [0] != '/') die ("Invalid option");

        switch (tolower (arg [1])) {
            case 'm':
                if (arg [2] != ':') die ("Invalid origin option");

                switch (tolower (arg [3])) {
                    case 'o': rhumbline = true; break;
                    case 'l': rhumbline = false; break;
                    default: die ("Invalid origin option");
                }

                break;
            
            case 'o':
                if (arg [2] != ':' || !parsePosition (origin, arg + 3)) die ("Invalid origin option");
                break;

            case 'd':
                if (arg [2] != ':' || !parsePosition (dest, arg + 3)) die ("Invalid destination option");
                break;

            case 'r':
                if (arg [2] == ':') range  = (float) atof (arg + 3);
                if (range < 0.0f) die ("Invalid range option");                
                break;

            case 'b':
                if (arg [2] == ':') bearing  = (float) atof (arg + 3);
                if (bearing < 0.0f || bearing >= 360.0f) die ("Invalid bearing option");                
                break;

            default:
                die ("Unknown option");
        }
    }

    switch (operation) {
        case geo::Operation::CALC_BRG_RNG: {
            if (origin.lat > 1.0e3 || origin.lon > 1.0e3) die ("Origin point not specified");
            if (dest.lat > 1.0e3 || dest.lon > 1.0e3) die ("Destination point not specified");

            SGeoPoint org { origin.lat, origin.lon }, dst { dest.lat, dest.lon };
            double rng, brg;

            gkdCalcRhumblineDistAndBrg (1, & org, & dst, & rng, & brg);
            gkdCalcRhumblineDistAndBrg2 (1, & org, & dst, & rng, & brg);

            geo::Pos _origin { geo::valToRad (origin.lat), geo::valToRad (origin.lon )};
            geo::Pos _dest { geo::valToRad (dest.lat), geo::valToRad (dest.lon )};

            bool result = rhumbline ? geo::calcRhumblineDistAndBrg (true, & _origin, & _dest, & range, & bearing) : false;

            if (result) {
                geo::radToDeg (& bearing);
                printf ("Leg from [%.6f; %.6f] to [%.6f; %.6f] is %.2fnm length; bearing is %.1fdeg\n", origin.lat, origin.lon, dest.lat, dest.lon, range, bearing);
                printf ("MARIS GU gives %.2fnm length; bearing is %.1fdeg\n", rng, brg);
            } else {
                printf ("Invalid data\n");
            }

            break;
        }
        case geo::Operation::CALC_DEST_POINT: {
            if (range < 0.0f) die ("Range not specified");
            if (bearing < 0.0f) die ("Bearing not specified");

            geo::Pos _origin { geo::valToRad (origin.lat), geo::valToRad (origin.lon )};            
            geo::Pos _dest;
            SGeoPoint org { origin.lat, origin.lon }, dst;            
            double brg = geo::valToRad (bearing);
            
            gkdCalcRhumblinePos (1, & org, range, bearing, & dst);

            if (geo::calcRhumblinePos (true, & _origin, range, brg, & _dest)) {
                geo::ptToDeg (& _dest);
                
                printf ("Leg from [%.6f; %.6f] %.2f nm length to %.1fdeg ends at [%.6f; %.6f]\n", origin.lat, origin.lon, range, bearing, _dest.lat, _dest.lon);
                printf ("MARIS GU gives [%.6f; %.6f]\n", dst.lat, dst.lon);
            } else  {
                printf ("Invalid data\n");
            }

            break;
        }
    }
}