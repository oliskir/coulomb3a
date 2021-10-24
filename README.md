# coulomb3a
Numerical solution of Newton's classical equations of motion for the three-alpha system implemented in Fortran.

[coul3a.f](coul3a.f) considers a fixed 8Be excitation energy, while [coul3a_bw.f](coul3a_bw.f) samples the 8Be excitation energy from a Breit-Wigner distribution.

These scripts can be used to estimate the effects of final-state Coulomb interactions on the 3a Dalitz distribution.

## Notes

Requires `cernlib`. Use [cernf](cernf) for compiling the Fortran scripts.
