clc

clear all
close all

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

params.Q_extraction = 4687.07;

params.t_simulation = 50;

com.comsol.model.util.ModelUtil.showProgress(true);

opts = optimset('display', 'iter', 'tolfun', 0.1, 'tolx', 1);

% model = init_simple(params);
% 
% Q_extraction1 = fminbnd(@(x)cost_function(model,x), 4000, 8000, opts);

model = init_complex(params);

Q_extraction2 = fminbnd(@(x)cost_function(model,x), 4000, 5000, opts);
