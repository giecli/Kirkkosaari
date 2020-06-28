clc

clear all
close all

params.d_borehole = 140e-3;
params.L_borehole = 1000;
params.T_surface = 4.6;
params.q_geothermal = 0.048;
params.k_rock = 3.2;
params.Cp_rock = 728;
params.rho_rock = 2750;
params.Cp_soil = 1500;
params.rho_soil = 1600;
params.k_clay = 1.86;
params.Cp_clay = 820;
params.rho_clay = 2560;
params.Q_extraction = 4687.07;
params.t_simulation = 50;
params.hmax_vertical_edge = 20e-3;
params.hmax_horizontal_edge = 10e-3;

H_soil = [20, 40, 60, 80, 100, 120];
H_clay = [50, 100, 200, 400, 600, 800, 1000];
k_soil = [1.7, 2.4];

com.comsol.model.util.ModelUtil.showProgress(true);

counter = 1;

options = optimset('display', 'iter', 'tolfun', 0.1, 'tolx', 1, 'plotfcn', {@plot_x, @plot_fval});

fid = fopen('results_bhe_1000m.csv', 'w');

fprintf(fid, '#;H_clay;H_soil;k_soil;Q_extraction\n');

for i = 1:length(H_clay)
    for j = 1:length(H_soil)
        for k = 1:length(k_soil)
            
            tic;
            
            params.H_clay = H_clay(i);
            params.H_soil = H_soil(j);
            params.k_soil = k_soil(k);
            
            model = init_complex(params);
            
            Q_extraction = fminbnd(@(x)cost_function(model,x), 7000, 13000, options);
            
            fprintf(fid, '%d;%.0f;%.0f;%.2f;%.3f\n', counter, H_clay(i), H_soil(j), k_soil(k), Q_extraction);
            
            t_elapsed = toc / 60;
            
            fprintf(1, 'i=%d H_clay=%.0f H_soil=%.0f k_soil=%.2f Q_extraction=%.3f t_elapsed=%.1f\n', counter, H_clay(i), H_soil(j), k_soil(k), Q_extraction, t_elapsed);
            
            counter = counter + 1;
            
        end
    end
end

fclose(fid);
