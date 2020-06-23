# Kirkkosaari

This repository contains MATLAB® code for the use of the Kirkkosaari project.

## COMSOL vs. EED

Optimizing the simple and complex COMSOL models with the following parameters

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

gives the result of 4687 W of constant heat extraction for 50 years
which corresponds to 41 MWh/a of energy being extracted.
However, for the same parameters, EED gives 51.2 MWh/a of energy for 50 years.
