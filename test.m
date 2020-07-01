clc

clear all
close all

params.R_model = 500;
params.H_model = 2500;

params.T_surface = 7.0;
params.q_geothermal = 0.042;

params.L_borehole = 2000;
params.d_borehole = 200e-3;

params.d_inner = 76e-3;
params.d_outer = 114e-3;

params.k_rock = 3.2;
params.Cp_rock = 728;
params.rho_rock = 2850;

params.H_soil = 20;
params.k_soil = params.k_rock;%1.7;
params.Cp_soil = params.Cp_rock;%1500;
params.rho_soil = params.rho_rock;%1600;

params.H_clay = 50;
params.k_clay = params.k_rock;%1.86;
params.Cp_clay = params.Cp_rock;%820;
params.rho_clay = params.rho_rock;%2560;

params.k_pipe = 0.02;
params.Cp_pipe = 1900;
params.rho_pipe = 950;

params.k_fluid = 0.6;
params.Cp_fluid = 4186;
params.rho_fluid = 1000;

params.Q_fluid = 3;
params.Q_extraction = 70000;

params.rtol = 1e-3;
params.h_buffer = 5;
params.t_simulation = 50;

com.comsol.model.util.ModelUtil.showProgress(true);

model = init_coaxial(params);

% mphsave(model, 'coaxial.mph');

% model.sol('sol1').runAll;

%options = optimset('display', 'iter', 'tolfun', 0.1, 'tolx', 10, 'plotfcn', {@plot_x, @plot_fval});
options = optimset('display', 'iter', 'tolfun', 0.1, 'tolx', 10, 'plotfcn', @plot_xy);

x0 = [100e3, 1.5];

% x = fminbnd(@(x)cost_function(model,x), 0.1, 10.0, options)
x = fminsearch(@(x)cost_function(model,x), x0, options)

% si = mphsolinfo(model);
% t = si.solvals / (365.2425 * 24 * 3600);
% T_min = mphglobal(model,'T_min','unit','degC');
% Q_wall = mphglobal(model, 'Q_wall', 'unit', 'W');
% figure;
% semilogx(t, T_min);
% figure;
% semilogx(t, Q_wall);

% a = linspace(0, 500000, 100);
% b = linspace(0, 25, 100);
% [a, b] = meshgrid(a, b);
% surf(a,b,rost_function(a,b))
% x = fminsearchbnd(@(x)rost_function(x), x0, lb, ub, options)
% x = fminunc(@(x)rost_function(x), x0, options)
