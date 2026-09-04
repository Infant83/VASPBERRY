#!/bin/sh

mpif90 -DMPI_USE -fpp -assume byterecl -o berry berry.f

mpirun -np 8 ./berry -z2 1 -kx 12 -ky 12 -s 2 -ii 1 -if 10

gnuplot gnuZ2.gp

eps2eps Z2.eps a.eps;mv a.eps Z2.eps
