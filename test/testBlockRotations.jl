import Oiler

using Test

@test Oiler.BlockRotations.build_PvGb_deg(0., 0.) == [ 0.0  -6.371e9  0.0;
                                  0.0   0.0      6.371e9;
                                  0.0   0.0      0.0]

@test Oiler.BlockRotations.build_Pv_deg(0., 0.) == [ -0.0  -0.0   1.0;
                                -0.0   1.0   0.0;
                                -1.0  -0.0  -0.0]

@test Oiler.BlockRotations.build_Gb_deg(0., 0.) == [  0.0   0.0      -0.0;
                                -0.0   0.0       6.371e9;
                                 0.0  -6.371e9   0.0]

