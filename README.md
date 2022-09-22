# Simulation-of-Kirkwood-Gaps
Repository which contains my code for simulating the Kirkwood Gaps in the Asteroid belt as part of IITB's Krittika Summer Projects 2021

Kirkwood gaps are gaps in the orbital periods of Main-Belt asteroids corresponding to orbital resonances with Jupiter

This code, written in FORTRAN uses the 'Euler Richardson' differential equation solver to calculate positions of randomly generated asteroids after a set number of years under the influence of the Jupiter-Sun system.

To reduce computational load, the problem is modelled as a reduced three body problem. The equations have been nondimensionalized as well.

Refer the attached project report for more details.

A few plots are also attached

## How to run
1. Compile Asteroids.f90 file using any fortran compiler
2. Plot the position data from the output file 'Asteroids.txt' as desired

Notes:
Uncomment the OMP tags in the code to enable parallel processing
*Open MP does not work on Windows*
*Use linux to achieve parallel processing*
