using Test

import Oiler

lax_lon = 118 + 24 / 60.;
lax_lat = 33 + 57 / 60.;

jfk_lat = 40 + 38 / 60.;
jfk_lon = 73 + 47 / 60.;


quien_sabe_coords = [-121.37575 36.93768;
                     -121.33534 36.89357;
                     -121.29735 36.85596;
                     -121.25795 36.81047;
                     -121.2073 36.75499];

@test Oiler.gc_distance(lax_lon, lax_lat, jfk_lon, jfk_lat) == 3972.857776250374

@test Oiler.azimuth(lax_lon, lax_lat, jfk_lon, jfk_lat) == 77.11787999624288

@test Oiler.average_azimuth(quien_sabe_coords[:,1], quien_sabe_coords[:,2]) == 324.1541623629684