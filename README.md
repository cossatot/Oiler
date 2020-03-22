# Oiler

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cossatot.gitlab.io/Oiler.jl/dev)
[![Build Status](https://gitlab.com/cossatot/Oiler.jl/badges/master/build.svg)](https://gitlab.com/cossatot/Oiler.jl/pipelines)
[![Coverage](https://gitlab.com/cossatot/Oiler.jl/badges/master/coverage.svg)](https://gitlab.com/cossatot/Oiler.jl/commits/master)


`Oiler` is a program for modeling tectonic block motions on a sphere, in
particular for solving for best-fit, kinematically-consistent fault slip
rates given both geologic fault slip rate estimates and GNSS-type geodetic
velocities. The primary motivation is in estimating a complete and internally-
consistent set of fault slip rates for all faults in a region for seismic
hazard assessment, though of course more traditional earth-science research
may find use here as well. Because of this motivation, `Oiler` is fault-first
in that it directly uses fault data (traces and other geometry, slip rates, kinematics)
produced by geologic investigation and required to build fault-source models
for seismic hazard analysis.

Foundationally, `Oiler` solves for the poles and rates of rotation between
adjacent blocks separated by faults, constrained by relative velocities
between blocks, through the relation

V = Omega x RP

where V represents the velocity of a point P relative to some block, Omega is
the rotation vector (pole and rate) describing relative block motion, and RP
is the Cartesian representation of the location P.  If three or more blocks
are present in a model, velocity closures are imposed so that all velocities
in the solution are internally consistent.

`Oiler` treats both fault slip rate estimates and GNSS data as velocities V,
differing only in two aspects: Fault slip rates may only constrain motions
between adjacent blocks separated by the fault (though any number of faults may
separate two blocks), while GNSS velocities may constrain the motions between
any block and the geodetic reference frame (which may or may not be an actual
block in the model). GNSS data are also subject to earthquake-cycle velocity
perturbations (i.e., the effects of interseismic fault locking) that geologic
fault velocities are not.

The natural incorporation of both fault slip rates and GNSS velocities is common
to (typically global) plate motion models such as MORVEL (DeMets et al., 2010),
but not to regional-scale elastic block models such as Blocks (Meade and
Loveless, 2009) and DEFNODE/TDEFNODE (McCaffrey 2002, 200x). In the case of
`Oiler`, the use of both velocity data types allows relative velocities between
blocks to be estimated regardless of whether any or all blocks contain GNSS
stations, which allows the modeler to incorporate many more faults and honor
fault network geometry, particularly in remote and poorly-studied areas where
GNSS data are sparse to nonexistent.