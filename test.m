clc

clear all
close all

params.R_model = 500;
params.H_model = 2500;

params.T_surface = 7.0;
params.q_geothermal = 0.040;

params.L_borehole = 2000;
params.d_borehole = 200e-3;

params.d_inner = 32e-3;
params.d_outer = 50e-3;

params.k_rock = 3.2;
params.Cp_rock = 728;
params.rho_rock = 2850;

params.H_soil = 20;
params.k_soil = 1.7;
params.Cp_soil = 1500;
params.rho_soil = 1600;

params.H_clay = 50;
params.k_clay = 1.86;
params.Cp_clay = 820;
params.rho_clay = 2560;

params.k_pipe = 0.1;
params.Cp_pipe = 1900;
params.rho_pipe = 950;

params.k_fluid = 0.6;
params.Cp_fluid = 4186;
params.rho_fluid = 1000;

params.Q_fluid = 0.6;
params.Q_extraction = 70000;

params.rtol = 1e-3;
params.istep = 1e-6;

params.h_buffer = 1;
params.t_simulation = 50;

com.comsol.model.util.ModelUtil.showProgress(true);

model = init_coaxial(params);

mphsave(model, 'coaxial.mph');
