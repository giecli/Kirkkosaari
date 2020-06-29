% =========================================================================
% Creates a new coaxial model.
% =========================================================================

function model = init_coaxial(params)

% -------------------------------------------------------------------------
% Imports Java packages.
% -------------------------------------------------------------------------

import com.comsol.model.*
import com.comsol.model.util.*

% -------------------------------------------------------------------------
% Creates a new model and initializes it with a component and a mesh.
% -------------------------------------------------------------------------

model = ModelUtil.create('Model');

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').axisymmetric(true);
model.component('comp1').mesh.create('mesh1');

% -------------------------------------------------------------------------
% Sets up model parameters.
% -------------------------------------------------------------------------

model.param.set('R_model', sprintf('%f[m]', params.R_model));
model.param.set('H_model', sprintf('%f[m]', params.H_model));

model.param.set('T_surface', sprintf('%f[degC]', params.T_surface));
model.param.set('q_geothermal', sprintf('%f[W/m^2]', params.q_geothermal));

model.param.set('L_borehole', sprintf('%f[m]', params.L_borehole));
model.param.set('r_borehole', sprintf('%f[m]', 0.5*params.d_borehole));

model.param.set('r_inner', sprintf('%f[m]', 0.5*params.d_inner));
model.param.set('r_outer', sprintf('%f[m]', 0.5*params.d_outer));

model.param.set('k_rock', sprintf('%f[W/(m*K)]', params.k_rock));
model.param.set('Cp_rock', sprintf('%f[J/(kg*K)]', params.Cp_rock));
model.param.set('rho_rock', sprintf('%f[kg/m^3]', params.rho_rock));

model.param.set('H_soil', sprintf('%f[m]', params.H_soil));
model.param.set('k_soil', sprintf('%f[W/(m*K)]', params.k_soil));
model.param.set('Cp_soil', sprintf('%f[J/(kg*K)]', params.Cp_soil));
model.param.set('rho_soil', sprintf('%f[kg/m^3]', params.rho_soil));

model.param.set('H_clay', sprintf('%f[m]', params.H_clay));
model.param.set('k_clay', sprintf('%f[W/(m*K)]', params.k_clay));
model.param.set('Cp_clay', sprintf('%f[J/(kg*K)]', params.Cp_clay));
model.param.set('rho_clay', sprintf('%f[kg/m^3]', params.rho_clay));

model.param.set('k_pipe', sprintf('%f[W/(m*K)]', params.k_pipe));
model.param.set('Cp_pipe', sprintf('%f[J/(kg*K)]', params.Cp_pipe));
model.param.set('rho_pipe', sprintf('%f[kg/m^3]', params.rho_pipe));

model.param.set('k_fluid', sprintf('%f[W/(m*K)]', params.k_fluid));
model.param.set('Cp_fluid', sprintf('%f[J/(kg*K)]', params.Cp_fluid));
model.param.set('rho_fluid', sprintf('%f[kg/m^3]', params.rho_fluid));

model.param.set('Q_fluid', sprintf('%f[L/s]', params.Q_fluid));
model.param.set('Q_extraction', sprintf('%f[W]', params.Q_extraction));

model.param.set('A_inner', 'pi*r_inner^2');
model.param.set('A_outer', 'pi*r_borehole^2-pi*r_outer^2');

model.param.set('v_inner', 'Q_fluid/A_inner');
model.param.set('v_outer', 'Q_fluid/A_outer');

model.param.set('delta_T', 'Q_extraction/(rho_fluid*Cp_fluid*Q_fluid)');

model.param.set('h_buffer', sprintf('%f[m]', params.h_buffer));

% -------------------------------------------------------------------------
% Creates the initial temperature function.
% -------------------------------------------------------------------------

model.func.create('an1', 'Analytic');
model.func('an1').set('funcname', 'T_initial');
model.func('an1').set('expr', 'T_surface-q_geothermal/k_rock*z');
model.func('an1').set('args', {'z'});
model.func('an1').set('argunit', 'm');
model.func('an1').set('fununit', 'K');

model.func.create('gp1', 'GaussianPulse');
model.func('gp1').set('location', '-0.5*L_borehole');
model.func('gp1').set('sigma', '0.2*L_borehole');
model.func('gp1').set('normalization', 'peak');

model.func.create('an2', 'Analytic');
model.func('an2').set('expr', 'max(min(gp1(z[1/m]),0.5),0.001)');
model.func('an2').set('args', {'z'});
model.func('an2').set('argunit', 'm');
model.func('an2').set('fununit', 'K');
model.func('an2').set('plotargs', {'z' '-L_borehole' '0'});

model.func.create('step1', 'Step');
model.func('step1').set('location', '1/365');
model.func('step1').set('smooth', '2/365');

% -------------------------------------------------------------------------
% Creates the geometry and runs it.
% -------------------------------------------------------------------------

model.component('comp1').geom('geom1').label('Geometry');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'0' '-H_soil'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'R_model' 'H_soil'});

model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0' '-H_soil-H_clay'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'R_model' 'H_clay'});

model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'0' '-H_model'});
model.component('comp1').geom('geom1').feature('r3').set('size', {'R_model' 'H_model-H_soil-H_clay'});

model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', {'0' '-L_borehole'; 'r_borehole' '-L_borehole'; 'r_borehole' '0'; 'r_borehole+h_buffer' '0'; 'r_borehole+h_buffer' '-L_borehole-h_buffer'; '0' '-L_borehole-h_buffer'});

model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'0' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'r_inner' 'L_borehole'});

model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').set('pos', {'r_inner' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r5').set('size', {'r_outer-r_inner' 'L_borehole'});

model.component('comp1').geom('geom1').create('r6', 'Rectangle');
model.component('comp1').geom('geom1').feature('r6').set('pos', {'r_outer' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r6').set('size', {'r_borehole-r_outer' 'L_borehole'});

model.component('comp1').geom('geom1').run('fin');

% -------------------------------------------------------------------------
% Creates component couplings and variables.
% -------------------------------------------------------------------------

i = get_boundaries([0.5*params.d_borehole 0], [0.5*params.d_borehole -params.L_borehole]);

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').label('Borehole Wall Minimum Operator');
model.component('comp1').cpl('minop1').selection.geom('geom1', 1);
model.component('comp1').cpl('minop1').selection.set(i);

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').label('Borehole Wall Integration Operator');
model.component('comp1').cpl('intop1').selection.geom('geom1', 1);
model.component('comp1').cpl('intop1').selection.set(i);
model.component('comp1').cpl('intop1').set('axisym', true);

i = get_boundaries([0.5*params.d_outer 0], [0.5*params.d_borehole 0]);

model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').label('Top Inlet Average Operator');
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop1').selection.set(i);
model.component('comp1').cpl('aveop1').set('axisym', true);

i = get_boundaries([0.5*params.d_outer -params.L_borehole], [0.5*params.d_borehole -params.L_borehole]);

model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl('aveop2').label('Bottom Outlet Average Operator');
model.component('comp1').cpl('aveop2').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop2').selection.set(i);
model.component('comp1').cpl('aveop2').set('axisym', true);

i = get_boundaries([0 0], [0.5*params.d_inner 0]);

model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('aveop3').label('Top Outlet Average Operator');
model.component('comp1').cpl('aveop3').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop3').selection.set(i);
model.component('comp1').cpl('aveop3').set('axisym', true);

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('T_min', 'minop1(T)');
model.component('comp1').variable('var1').set('Q_wall', 'intop1(ht.ndflux)');
model.component('comp1').variable('var1').set('T_inlet', 'aveop1(T)');
model.component('comp1').variable('var1').set('T_bottom', 'aveop2(T)');
model.component('comp1').variable('var1').set('T_outlet', 'aveop3(T)');

% -------------------------------------------------------------------------
% Creates the physics.
% -------------------------------------------------------------------------

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');

model.component('comp1').physics('ht').prop('ShapeProperty').set('order_temperature', 1);

model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock'; '0'; '0'; '0'; 'k_rock'; '0'; '0'; '0'; 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');

model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');

i = get_boundaries([0.5*params.d_borehole 0], [params.R_model 0]);

model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').label('Surface Heat Flux');
model.component('comp1').physics('ht').feature('hf1').selection.set(i);
model.component('comp1').physics('ht').feature('hf1').set('q0', '-q_geothermal');

i = get_boundaries([0 -params.H_model], [params.R_model -params.H_model]);

model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').label('Bottom Heat Flux');
model.component('comp1').physics('ht').feature('hf2').selection.set(i);
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');

i = get_domains([0.5*params.d_borehole 0], [params.R_model -params.H_soil]);

model.component('comp1').physics('ht').create('solid2', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid2').label('Soil Solid');
model.component('comp1').physics('ht').feature('solid2').selection.set(i);
model.component('comp1').physics('ht').feature('solid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('k', {'k_soil'; '0'; '0'; '0'; 'k_soil'; '0'; '0'; '0'; 'k_soil'});
model.component('comp1').physics('ht').feature('solid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('rho', 'rho_soil');
model.component('comp1').physics('ht').feature('solid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('Cp', 'Cp_soil');

i = get_domains([0.5*params.d_borehole -params.H_soil], [params.R_model -params.H_soil-params.H_clay]);

model.component('comp1').physics('ht').create('solid3', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid3').label('Claystone Solid');
model.component('comp1').physics('ht').feature('solid3').selection.set(i);
model.component('comp1').physics('ht').feature('solid3').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('k', {'k_clay'; '0'; '0'; '0'; 'k_clay'; '0'; '0'; '0'; 'k_clay'});
model.component('comp1').physics('ht').feature('solid3').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('rho', 'rho_clay');
model.component('comp1').physics('ht').feature('solid3').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('Cp', 'Cp_clay');

i = get_domains([0 0], [0.5*params.d_inner -params.L_borehole]);

model.component('comp1').physics('ht').create('fluid1', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid1').label('Inner Fluid');
model.component('comp1').physics('ht').feature('fluid1').selection.set(i);
model.component('comp1').physics('ht').feature('fluid1').set('u', {'0'; '0'; 'step1(t[1/a])*v_inner'});
model.component('comp1').physics('ht').feature('fluid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('k', {'1000' '0' '0' '0' '1000' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('gamma', '1');

i = get_domains([0.5*params.d_inner 0], [0.5*params.d_outer -params.L_borehole]);

model.component('comp1').physics('ht').create('solid4', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid4').label('Pipe Solid');
model.component('comp1').physics('ht').feature('solid4').selection.set(i);
model.component('comp1').physics('ht').feature('solid4').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('k', {'k_pipe'; '0'; '0'; '0'; 'k_pipe'; '0'; '0'; '0'; 'k_pipe'});
model.component('comp1').physics('ht').feature('solid4').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('rho', 'rho_pipe');
model.component('comp1').physics('ht').feature('solid4').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('Cp', 'Cp_pipe');

i = get_domains([0.5*params.d_outer 0], [0.5*params.d_borehole -params.L_borehole]);

model.component('comp1').physics('ht').create('fluid2', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid2').label('Outer Fluid');
model.component('comp1').physics('ht').feature('fluid2').selection.set(i);
model.component('comp1').physics('ht').feature('fluid2').set('u', {'0'; '0'; '-step1(t[1/a])*v_outer'});
model.component('comp1').physics('ht').feature('fluid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('k', {'1000' '0' '0' '0' '1000' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('gamma', '1');

i = get_boundaries([0.5*params.d_outer 0], [0.5*params.d_borehole 0]);

model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp2').label('Top Inlet Temperature');
model.component('comp1').physics('ht').feature('temp2').selection.set(i);
model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-step1(t[1/a])*delta_T');

i = get_boundaries([0 -params.L_borehole], [0.5*params.d_inner -params.L_borehole]);

model.component('comp1').physics('ht').create('temp3', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp3').label('Bottom Inlet Temperature');
model.component('comp1').physics('ht').feature('temp3').selection.set(i);
model.component('comp1').physics('ht').feature('temp3').set('T0', 'T_bottom');

% -------------------------------------------------------------------------
% Creates the mesh and runs it.
% -------------------------------------------------------------------------

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.001);

i = get_boundaries([0 0], [0.5*params.d_borehole 0]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set(i);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '2.5[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

i = get_boundaries([0 0], [0 -params.L_borehole]);
j = get_boundaries([0.5*params.d_inner 0], [0.5*params.d_inner -params.L_borehole]);
k = get_boundaries([0.5*params.d_outer 0], [0.5*params.d_outer -params.L_borehole]);
l = get_boundaries([0.5*params.d_borehole 0], [0.5*params.d_borehole -params.L_borehole]);

model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([i j k l]);
model.component('comp1').mesh('mesh1').feature('edg2').create('se1', 'SizeExpression');
model.component('comp1').mesh('mesh1').feature('edg2').feature('se1').set('sizeexpr', 'an2(z)');

i = get_domains([0 0], [0.5*params.d_borehole -params.L_borehole]);

model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set(i);

i = get_domains([0.5*params.d_borehole 0], [0.5*params.d_borehole+params.h_buffer -params.H_soil]);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set(i);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');

i = get_domains([0.5*params.d_borehole -params.H_soil], [0.5*params.d_borehole+params.h_buffer -params.H_soil-params.H_clay]);

model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set(i);
model.component('comp1').mesh('mesh1').feature('ftri2').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 1);

i = get_domains([0 -params.H_soil-params.H_clay], [0.5*params.d_borehole+params.h_buffer -params.L_borehole-params.h_buffer]);
j = get_domains([0 -params.H_soil-params.H_clay], [0.5*params.d_borehole -params.L_borehole]);
k = setdiff(i, j);

model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set(k);
model.component('comp1').mesh('mesh1').feature('ftri3').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 1);

i = get_domains([0.5*params.d_borehole+params.h_buffer 0], [params.R_model -params.H_soil]);

model.component('comp1').mesh('mesh1').create('ftri4', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri4').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri4').selection.set(i);
model.component('comp1').mesh('mesh1').feature('ftri4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri4').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmax', 5);

i = get_domains([0.5*params.d_borehole+params.h_buffer -params.H_soil], [params.R_model -params.H_soil-params.H_clay]);

model.component('comp1').mesh('mesh1').create('ftri5', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri5').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri5').selection.set(i);
model.component('comp1').mesh('mesh1').feature('ftri5').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri5').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmax', 25);

model.component('comp1').mesh('mesh1').create('ftri6', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri6').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri6').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hauto', 1);

model.component('comp1').mesh('mesh1').run;

% -------------------------------------------------------------------------
% Creates study and solution.
% -------------------------------------------------------------------------

tlist = sprintf('0 %d', params.t_simulation);

model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', tlist);
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', params.rtol);

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('clist', {tlist '1e-6[a]'});
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').set('rtol', params.rtol);
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('dDef').set('pivotperturb', 1.0e-13);

%model.sol('sol1').runAll;

    function i = get_boundaries(pt1, pt2)
        i = mphselectbox(model, 'geom1', [pt1(1) pt2(1); pt1(2) pt2(2)], 'boundary');
        assert(~isempty(i));
    end

    function i = get_domains(pt1, pt2)
        i = mphselectbox(model, 'geom1', [pt1(1)-0.001 pt2(1)+0.001; pt1(2)+0.001 pt2(2)-0.001], 'domain');
        assert(~isempty(i));
    end

end
