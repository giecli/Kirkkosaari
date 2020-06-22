clc

clear all
close all

T_surface = [4.6 5.5];
q_geothermal = [0.042 0.042];



global model

for i = 1:2

    param.T_surface = T_surface(i);
    param.q_geothermal = q_geothermal(i);
    param.k_rock = 3.2;
    param.Cp_rock = 728;
    param.rho_rock = 2750;
    param.Q_extraction = 8000;
    
    model = init_simple(param);

    opts = optimset('display', 'iter');

    Q = fminbnd(@cost_function, 1, 100000, opts);

end
