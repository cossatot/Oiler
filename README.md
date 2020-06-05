# Oiler


[![Build
Status](https://travis-ci.com/cossatot/Oiler.svg?branch=master)](https://travis-ci.com/cossatot/Oiler)

[![Coverage Status](https://coveralls.io/repos/github/cossatot/Oiler/badge.svg?branch=master)](https://coveralls.io/github/cossatot/Oiler?branch=master)

*Note: This project is still in rapid development and not sufficiently tested
for use by those other than me.*

`Oiler` is a program for modeling tectonic block motions on a sphere, in
particular for solving for best-fit, kinematically-consistent fault slip
rates given both geologic fault slip rate estimates and GNSS-type geodetic
velocities. The primary motivation is in estimating a complete and
internally- consistent set of fault slip rates for all faults in a region for
seismic hazard assessment, though of course more traditional earth-science
research may find use here as well. Because of this motivation, `Oiler` is
fault-first in that it directly uses fault data (traces and other geometry,
slip rates, kinematics) produced by geologic investigation and required to
build fault-source models for seismic hazard analysis.

`Oiler` aims to be fast and (sort of) interactive, so that a scientist can
progressively build a fault and block network in a GIS program, quantifying
fault slip rate estimates given data and network details, and adjusting weights,
slip rates, fault connectivity and block geometry iteratively based on rapid
modeling results from `Oiler`. Once the workflow and IO are streamlined, more
details on this will be available.

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

The use of both velocity data types allows relative velocities between
blocks to be estimated regardless of whether any or all blocks contain GNSS
stations, which allows the modeler to incorporate many more faults and honor
fault network geometry, particularly in remote and poorly-studied areas where
GNSS data are sparse to nonexistent.

Oiler draws heavily from the previous research into plate/block theory and to
a lesser extent, previous implementations. The mathematical formulations of
Meade and Loveless (2009), Chase (1972), and Cox and Hart (1986) were extremely
helpful, and I am very grateful to those authors for their careful and thorough
explanations. Additionally, implementation details of certain parts of the
code here, particularly fault locking or earthquake cycle effects, were modeled
after the Blocks code by Meade and Loveless. I double thank Meade and Loveless 
for making their code public and with a permissive license.