# numba-jit-example
Small project illustrating the use of Anaconda's numba jit compiler to speed up calculations

There are two python files:  

* without.py -- the basic calculation.  This does the numerical integration, solves for the motions, and creates a plot showing the trajectories. It takes about 94.54 seconds to run on a 3.4 GHz Core i7.

* with.py -- the same calculation modified for Anaconda's Numba JIT compiler.  The changes are literally trivial -- just a matter of adding a single line in front o f most of the functions in program. The calculation is done, the motions are solved for, and a plot is created.  This only takes 46.05 seconds to run on the very same 3.4 GHz Core i7.

Wow!  So just by adding 8 trivial descriptive lines that provided hints to the 'numba' JIT compiler, we were able to more than DOUBLE the speed of the code!

The project doesn't really do anything terribly deep -- it's a basic physics ballistics calculation -- cannon shells traveling through a plane-parallel atmosphere in the presence of gravity and atmospheric drag.  It sets up and performs the calculation for 12 different cannon shells with the same size and shape but drastically different masses (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, and 5000 kgm) and launches them all from an identical cannon with the same initial velocity and then calculates the altitude and range they travel.  All of the shells are fired at a 45 degree angle to the horizon with a starting speed of 330 m/s.Notice that at first, increasing the mass of the projectile leads to longer ranges, but eventually there is vanishing benefit (that's because the projectile ends up with so much momentum that drag becomes irrelevant).

