1. run VASP calculations to get proper WAVECAR file (do one-shot SCF calculations with given KPOINTS)
2. Note that KPOINTs should include Full BZ (ISYM = -1)
3. run vaspberry to get BERRYCURV.dat:
   ./vaspberry -kx 12 -ky 12 -ii 1 -if 18 -s 2
5. run contour.py : python contour.py   (to get contour plot for the Berry curvature map)

