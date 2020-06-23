# Kirkkosaari

This repository contains MATLAB® code for the use of the Kirkkosaari project.

## Benchmark

Use the following parameters:

```
params.d_borehole = 140e-3;
params.L_borehole = 300;

params.T_surface = 4.6;
params.q_geothermal = 0.042;

params.k_rock = 3.2;
params.Cp_rock = 728;
params.rho_rock = 2750;

params.H_soil = 20;
params.k_soil = 3.2;
params.Cp_soil = 728;
params.rho_soil = 2750;

params.H_clay = 50;
params.k_clay = 3.2;
params.Cp_clay = 728;
params.rho_clay = 2750;
```

Minimization of **the average borehole wall temperature (`T_ave`)**
should give around 5928 W of constant heat extraction for 50 years
which corresponds to 5928 W × 8760 h = 51.9 MWh of energy each year.

Similarly, minimization of **the minimum borehole wall temperature (`T_min`)**
should give around 4687 W of constant heat extraction for 50 years
which corresponds to 4687 W × 8760 h = 41.1 MWh of energy each year.

Setting up an EED model with the above parameters and solving the maximal
extractable energy while not allowing the borehole wall temperature drop
below 0 centigrades will give 51.2 MWh/a.
