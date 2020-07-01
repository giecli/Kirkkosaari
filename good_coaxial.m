function out = model
%
% good_coaxial.m
%
% Model exported on Jul 1 2020, 14:38 by COMSOL 5.5.0.292.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Work\Kirkkosaari');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.param.set('R_model', '500.000000[m]');
model.param.set('H_model', '2500.000000[m]');
model.param.set('T_surface', '7.000000[degC]');
model.param.set('q_geothermal', '0.042000[W/m^2]');
model.param.set('L_borehole', '2000.000000[m]');
model.param.set('r_borehole', '0.100000[m]');
model.param.set('r_inner', '0.038000[m]');
model.param.set('r_outer', '0.057000[m]');
model.param.set('k_rock', '3.200000[W/(m*K)]');
model.param.set('Cp_rock', '728.000000[J/(kg*K)]');
model.param.set('rho_rock', '2850.000000[kg/m^3]');
model.param.set('H_soil', '20.000000[m]');
model.param.set('k_soil', '3.200000[W/(m*K)]');
model.param.set('Cp_soil', '728.000000[J/(kg*K)]');
model.param.set('rho_soil', '2850.000000[kg/m^3]');
model.param.set('H_clay', '50.000000[m]');
model.param.set('k_clay', '3.200000[W/(m*K)]');
model.param.set('Cp_clay', '728.000000[J/(kg*K)]');
model.param.set('rho_clay', '2850.000000[kg/m^3]');
model.param.set('k_pipe', '0.020000[W/(m*K)]');
model.param.set('Cp_pipe', '1900.000000[J/(kg*K)]');
model.param.set('rho_pipe', '950.000000[kg/m^3]');
model.param.set('k_fluid', '0.600000[W/(m*K)]');
model.param.set('Cp_fluid', '4186.000000[J/(kg*K)]');
model.param.set('rho_fluid', '1000.000000[kg/m^3]');
model.param.set('Q_fluid', '3.000000[L/s]');
model.param.set('Q_extraction', '70000.000000[W]');
model.param.set('A_inner', 'pi*r_inner^2');
model.param.set('A_outer', 'pi*r_borehole^2-pi*r_outer^2');
model.param.set('v_inner', 'Q_fluid/A_inner');
model.param.set('v_outer', 'Q_fluid/A_outer');
model.param.set('delta_T', 'Q_extraction/(rho_fluid*Cp_fluid*Q_fluid)');
model.param.set('h_buffer', '1.000000[m]');

model.func.create('pw1', 'Piecewise');
model.func('pw1').set('funcname', 'T_initial');
model.func('pw1').set('arg', 'z');
model.func('pw1').set('extrap', 'interior');
model.func('pw1').set('pieces', {'-H_soil' '0' 'T_surface-q_geothermal/k_soil*z'; '-H_soil-H_clay' '-H_soil' 'T_surface+q_geothermal/k_soil*H_soil-q_geothermal/k_clay*(z+H_soil)'; '-H_model' '-H_soil-H_clay' 'T_surface+q_geothermal/k_soil*H_soil+q_geothermal/k_clay*H_clay-q_geothermal/k_rock*(z+H_soil+H_clay)'});
model.func('pw1').set('argunit', 'm');
model.func('pw1').set('fununit', 'K');
model.func.create('gp1', 'GaussianPulse');
model.func('gp1').set('funcname', 'gaussian_pulse');
model.func('gp1').set('location', '-0.5*L_borehole');
model.func('gp1').set('sigma', '0.2*L_borehole');
model.func('gp1').set('normalization', 'peak');
model.func.create('an2', 'Analytic');
model.func('an2').set('expr', 'max(min(gaussian_pulse(z[1/m]),0.5),0.001)');
model.func('an2').set('funcname', 'mesh_size');
model.func('an2').set('args', 'z');
model.func('an2').set('argunit', 'm');
model.func('an2').set('fununit', 'K');
model.func('an2').set('plotargs', {'z' '-L_borehole' '0'});
model.func.create('step1', 'Step');
model.func('step1').set('funcname', 'ramp_function');
model.func('step1').set('location', '1/365.2425');
model.func('step1').set('smooth', '1/182.62125');

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
model.component('comp1').geom('geom1').feature('pol1').set('table', {'0' '-L_borehole';  ...
'r_borehole' '-L_borehole';  ...
'r_borehole' '0';  ...
'r_borehole+h_buffer' '0';  ...
'r_borehole+h_buffer' '-L_borehole-h_buffer';  ...
'0' '-L_borehole-h_buffer'});
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
model.component('comp1').geom('geom1').run;

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').label('Borehole Wall Minimum Operator');
model.component('comp1').cpl('minop1').selection.geom('geom1', 1);
model.component('comp1').cpl('minop1').selection.set([26 27 29]);
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').label('Borehole Wall Integration Operator');
model.component('comp1').cpl('intop1').selection.geom('geom1', 1);
model.component('comp1').cpl('intop1').selection.set([26 27 29]);
model.component('comp1').cpl('intop1').set('axisym', true);

model.component('comp1').geom('geom1').run;

model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').label('Top Inlet Average Operator');
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop1').selection.set([25]);
model.component('comp1').cpl('aveop1').set('axisym', true);

model.component('comp1').geom('geom1').run;

model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl('aveop2').label('Bottom Outlet Average Operator');
model.component('comp1').cpl('aveop2').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop2').selection.set([20]);
model.component('comp1').cpl('aveop2').set('axisym', true);

model.component('comp1').geom('geom1').run;

model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('aveop3').label('Top Outlet Average Operator');
model.component('comp1').cpl('aveop3').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop3').selection.set([11]);
model.component('comp1').cpl('aveop3').set('axisym', true);

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('T_min', 'minop1(T)');
model.component('comp1').variable('var1').set('Q_wall', 'intop1(ht.ndflux)');
model.component('comp1').variable('var1').set('T_inlet', 'aveop1(T)');
model.component('comp1').variable('var1').set('T_bottom', 'aveop2(T)');
model.component('comp1').variable('var1').set('T_outlet', 'aveop3(T)');

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').prop('ShapeProperty').set('order_temperature', 1);
model.component('comp1').physics('ht').feature('solid1').label('Bedrock Solid');
model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock' '0' '0' '0' 'k_rock' '0' '0' '0' 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');
model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').label('Surface Heat Flux BC');
model.component('comp1').physics('ht').feature('hf1').selection.set([31 37]);
model.component('comp1').physics('ht').feature('hf1').set('q0', '-q_geothermal');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').label('Bottom Heat Flux BC');
model.component('comp1').physics('ht').feature('hf2').selection.set([2]);
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('solid2', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid2').label('Soil Solid');
model.component('comp1').physics('ht').feature('solid2').selection.set([13 15]);
model.component('comp1').physics('ht').feature('solid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('k', {'k_soil' '0' '0' '0' 'k_soil' '0' '0' '0' 'k_soil'});
model.component('comp1').physics('ht').feature('solid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('rho', 'rho_soil');
model.component('comp1').physics('ht').feature('solid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('Cp', 'Cp_soil');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('solid3', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid3').label('Claystone Solid');
model.component('comp1').physics('ht').feature('solid3').selection.set([12 14]);
model.component('comp1').physics('ht').feature('solid3').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('k', {'k_clay' '0' '0' '0' 'k_clay' '0' '0' '0' 'k_clay'});
model.component('comp1').physics('ht').feature('solid3').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('rho', 'rho_clay');
model.component('comp1').physics('ht').feature('solid3').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('Cp', 'Cp_clay');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('fluid1', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid1').label('Inner Fluid');
model.component('comp1').physics('ht').feature('fluid1').selection.set([3 4 5]);
model.component('comp1').physics('ht').feature('fluid1').set('u', {'0' '0' 'ramp_function(t[1/a])*v_inner'});
model.component('comp1').physics('ht').feature('fluid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('k', {'1000' '0' '0' '0' '1000' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('gamma', '1');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('solid4', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid4').label('Pipe Solid');
model.component('comp1').physics('ht').feature('solid4').selection.set([6 7 8]);
model.component('comp1').physics('ht').feature('solid4').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('k', {'k_pipe' '0' '0' '0' 'k_pipe' '0' '0' '0' 'k_pipe'});
model.component('comp1').physics('ht').feature('solid4').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('rho', 'rho_pipe');
model.component('comp1').physics('ht').feature('solid4').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('Cp', 'Cp_pipe');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('fluid2', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid2').label('Outer Fluid');
model.component('comp1').physics('ht').feature('fluid2').selection.set([9 10 11]);
model.component('comp1').physics('ht').feature('fluid2').set('u', {'0' '0' '-ramp_function(t[1/a])*v_outer'});
model.component('comp1').physics('ht').feature('fluid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('k', {'1000' '0' '0' '0' '1000' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('gamma', '1');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp2').label('Top Inlet Temperature BC');
model.component('comp1').physics('ht').feature('temp2').selection.set([25]);
model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-ramp_function(t[1/a])*delta_T');

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ht').create('temp3', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp3').label('Bottom Inlet Temperature BC');
model.component('comp1').physics('ht').feature('temp3').selection.set([6]);
model.component('comp1').physics('ht').feature('temp3').set('T0', 'T_bottom');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.001);

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([11 18 25]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '2.5[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([5 7 9 12 14 16 19 26 27 29]);
model.component('comp1').mesh('mesh1').feature('edg2').create('se1', 'SizeExpression');
model.component('comp1').mesh('mesh1').feature('edg2').feature('se1').set('sizeexpr', 'mesh_size(z)');

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([3 4 5 6 7 8 9 10 11]);

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([13]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([12]);
model.component('comp1').mesh('mesh1').feature('ftri2').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 1);

model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri3').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 1);

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('ftri4', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri4').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri4').selection.set([15]);
model.component('comp1').mesh('mesh1').feature('ftri4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri4').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmax', 5);

model.component('comp1').geom('geom1').run;

model.component('comp1').mesh('mesh1').create('ftri5', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri5').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri5').selection.set([14]);
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

model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', '0 50');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', 0.001);

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('clist', {'0 50' '1e-6[a]'});
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('dDef').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.param.set('Q_extraction', '100000.000000[W]');
model.param.set('Q_fluid', '2.000000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 78);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '105000.000000[W]');
model.param.set('Q_fluid', '2.000000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '100000.000000[W]');
model.param.set('Q_fluid', '2.100000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '105000.000000[W]');
model.param.set('Q_fluid', '1.900000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 85);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '85');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '101250.000000[W]');
model.param.set('Q_fluid', '2.050000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '96250.000000[W]');
model.param.set('Q_fluid', '2.050000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '91875.000000[W]');
model.param.set('Q_fluid', '2.075000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 76);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '76');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95000.000000[W]');
model.param.set('Q_fluid', '2.000000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '91875.000000[W]');
model.param.set('Q_fluid', '1.975000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 74);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '74');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '91250.000000[W]');
model.param.set('Q_fluid', '2.050000[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '93437.500000[W]');
model.param.set('Q_fluid', '2.037500[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '97812.500000[W]');
model.param.set('Q_fluid', '2.012500[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 82);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '82');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '94531.250000[W]');
model.param.set('Q_fluid', '2.031250[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 76);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '76');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '93281.250000[W]');
model.param.set('Q_fluid', '1.981250[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95507.812500[W]');
model.param.set('Q_fluid', '2.032812[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '94023.437500[W]');
model.param.set('Q_fluid', '1.998438[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 83);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '83');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95136.718750[W]');
model.param.set('Q_fluid', '2.024219[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95605.468750[W]');
model.param.set('Q_fluid', '1.992969[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '94799.804688[W]');
model.param.set('Q_fluid', '2.021680[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95336.914063[W]');
model.param.set('Q_fluid', '2.002539[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '94934.082031[W]');
model.param.set('Q_fluid', '2.016895[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95070.800781[W]');
model.param.set('Q_fluid', '2.041113[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95017.700195[W]');
model.param.set('Q_fluid', '2.010278[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95053.100586[W]');
model.param.set('Q_fluid', '2.030835[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 75);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '75');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95044.250488[W]');
model.param.set('Q_fluid', '2.025696[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95035.400391[W]');
model.param.set('Q_fluid', '2.020557[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95077.209473[W]');
model.param.set('Q_fluid', '2.017249[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95094.909668[W]');
model.param.set('Q_fluid', '2.027527[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95081.634521[W]');
model.param.set('Q_fluid', '2.019818[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95090.484619[W]');
model.param.set('Q_fluid', '2.024957[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95083.847046[W]');
model.param.set('Q_fluid', '2.021103[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95185.165405[W]');
model.param.set('Q_fluid', '2.024765[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 82);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '82');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95072.841644[W]');
model.param.set('Q_fluid', '2.021609[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95019.969940[W]');
model.param.set('Q_fluid', '2.018493[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95107.531548[W]');
model.param.set('Q_fluid', '2.022787[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 82);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '82');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95078.344345[W]');
model.param.set('Q_fluid', '2.021356[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95110.282898[W]');
model.param.set('Q_fluid', '2.022661[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 81);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '81');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95051.908493[W]');
model.param.set('Q_fluid', '2.019798[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95095.689297[W]');
model.param.set('Q_fluid', '2.021945[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 79);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '79');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95101.191998[W]');
model.param.set('Q_fluid', '2.021692[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95084.056258[W]');
model.param.set('Q_fluid', '2.021440[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 78);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95072.214007[W]');
model.param.set('Q_fluid', '2.020598[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 80);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '80');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95072.004795[W]');
model.param.set('Q_fluid', '2.020261[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95081.043392[W]');
model.param.set('Q_fluid', '2.021145[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 76);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '76');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95092.676431[W]');
model.param.set('Q_fluid', '2.021650[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 78);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95077.329613[W]');
model.param.set('Q_fluid', '2.020861[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 77);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '77');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95082.445219[W]');
model.param.set('Q_fluid', '2.021124[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 76);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '76');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.param.set('Q_extraction', '95078.030527[W]');
model.param.set('Q_fluid', '2.020850[L/s]');

model.sol('sol1').runAll;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);

model.sol('sol1').getSizeMulti;

model.result.numerical('global_internal').set('solnum', 78);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.sol('sol1').getSizeMulti;
model.sol('sol1').getSize;
model.sol('sol1').getPVals;
model.sol('sol1').getPValsImag;

model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'T_min');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);
model.result.numerical('global_internal').set('unit', {'degC'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('solnum', {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15' '16' '17' '18' '19' '20'  ...
'21' '22' '23' '24' '25' '26' '27' '28' '29' '30'  ...
'31' '32' '33' '34' '35' '36' '37' '38' '39' '40'  ...
'41' '42' '43' '44' '45' '46' '47' '48' '49' '50'  ...
'51' '52' '53' '54' '55' '56' '57' '58' '59' '60'  ...
'61' '62' '63' '64' '65' '66' '67' '68' '69' '70'  ...
'71' '72' '73' '74' '75' '76' '77' '78'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'T_min'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');
model.result.numerical.remove('global_internal');
model.result.numerical.create('global_internal', 'Global');
model.result.numerical('global_internal').set('expr', 'Q_wall');
model.result.numerical('global_internal').set('matherr', 'off');
model.result.numerical('global_internal').set('phase', 0);
model.result.numerical('global_internal').set('unit', {'W'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('solnum', {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15' '16' '17' '18' '19' '20'  ...
'21' '22' '23' '24' '25' '26' '27' '28' '29' '30'  ...
'31' '32' '33' '34' '35' '36' '37' '38' '39' '40'  ...
'41' '42' '43' '44' '45' '46' '47' '48' '49' '50'  ...
'51' '52' '53' '54' '55' '56' '57' '58' '59' '60'  ...
'61' '62' '63' '64' '65' '66' '67' '68' '69' '70'  ...
'71' '72' '73' '74' '75' '76' '77' '78'});
model.result.numerical('global_internal').set('outersolnum', 1);
model.result.numerical('global_internal').set('expr', {'Q_wall'});
model.result.numerical('global_internal').set('evalmethodactive', false);
model.result.numerical('global_internal').set('outertype', 'none');
model.result.numerical('global_internal').set('solvertype', 'time');
model.result.numerical('global_internal').set('outersolnum', '1');
model.result.numerical('global_internal').set('solnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78');
model.result.numerical('global_internal').set('timeinterp', 'off');
model.result.numerical('global_internal').getData;
model.result.numerical.remove('global_internal');

model.label('coaxial.mph');

model.func('step1').set('location', '2/12');
model.func('step1').set('smooth', '1/12');
model.func('step1').set('location', '1/12');
model.func('step1').set('smooth', '1/6');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset('rev1').label('Revolution 2D');
model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result.dataset('rev1').set('data', 'dset1');
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label('Temperature, 3D (ht)');
model.result('pg1').set('showlooplevel', {'off' 'off' 'off'});
model.result('pg1').set('data', 'rev1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('colortable', 'ThermalLight');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label('Isothermal Contours (ht)');
model.result('pg2').set('showlooplevel', {'off' 'off' 'off'});
model.result('pg2').set('dataisaxisym', 'off');
model.result('pg2').set('data', 'dset1');
model.result('pg2').feature.create('con1', 'Contour');
model.result('pg2').feature('con1').label('Contour');
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('colortable', 'ThermalLight');
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('data', 'parent');
model.result.remove('pg2');
model.result.remove('pg1');
model.result.dataset.remove('rev1');

model.func('step1').set('location', 1);
model.func('step1').set('smooth', 1);
model.func('step1').set('location', 0.5);

model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset('rev1').label('Revolution 2D');
model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result.dataset('rev1').set('data', 'dset1');
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label('Temperature, 3D (ht)');
model.result('pg1').set('showlooplevel', {'off' 'off' 'off'});
model.result('pg1').set('data', 'rev1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('colortable', 'ThermalLight');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label('Isothermal Contours (ht)');
model.result('pg2').set('showlooplevel', {'off' 'off' 'off'});
model.result('pg2').set('dataisaxisym', 'off');
model.result('pg2').set('data', 'dset1');
model.result('pg2').feature.create('con1', 'Contour');
model.result('pg2').feature('con1').label('Contour');
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('colortable', 'ThermalLight');
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('data', 'parent');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result.create('pg3', 'PlotGroup1D');
model.result('pg3').run;
model.result('pg3').create('glob1', 'Global');
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', false);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', true);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').run;
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', false);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', true);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').feature('glob1').set('linemarker', 'point');
model.result('pg3').feature('glob1').set('markerpos', 'datapoints');

model.component('comp1').physics('ht').prop('InconsistentStabilization').set('HeatIsotropicDiffusion', true);

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;

model.func('an2').set('fununit', 'm');

model.result('pg3').run;

model.func('step1').set('location', '1/12');
model.func('step1').set('smooth', '1/12');
model.func('step1').set('location', '1/24');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;

model.component('comp1').physics('ht').feature('fluid1').set('u', {'0' '0' 'v_inner'});
model.component('comp1').physics('ht').feature('fluid2').set('u', {'0' '0' '-v_outer'});

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-ramp_function(t[1/a])*delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;

model.study('std1').feature('time').set('rtol', '0.0001');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;

model.study('std1').feature('time').set('rtol', 0.01);

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').set('xlog', false);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', true);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', false);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', true);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);

model.study('std1').feature('time').set('rtol', 0.001);

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-ramp_function(t[1/a])*delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-ramp_function(t[1/a])*delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;

model.func('step1').set('location', '2/24');
model.func('step1').set('smooth', '2/12');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;

model.study('std1').feature('time').set('rtol', '1e-2');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;

model.param.set('Q_fluid', '13[L/s]');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'T_min', 0);
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg3').run;

model.param.set('Q_fluid', '2[L/s]');
model.param.set('Q_extraction', '95[kW]');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg3').run;
model.result.create('pg4', 'PlotGroup1D');
model.result('pg4').run;
model.result('pg4').create('glob1', 'Global');
model.result('pg4').feature('glob1').setIndex('expr', 'T_inlet', 0);
model.result('pg4').feature('glob1').setIndex('unit', 'degC', 0);
model.result('pg4').feature('glob1').setIndex('expr', 'T_outlet', 1);
model.result('pg4').feature('glob1').setIndex('unit', 'degC', 1);
model.result('pg4').run;
model.result('pg4').set('xlog', true);
model.result('pg4').set('showmanualgrid', true);
model.result('pg4').set('showxspacing', false);
model.result('pg4').set('showyspacing', true);
model.result('pg4').set('showsecyspacing', false);
model.result('pg4').set('showsecyextra', false);
model.result('pg4').set('xlog', false);
model.result('pg4').set('showmanualgrid', true);
model.result('pg4').set('showxspacing', true);
model.result('pg4').set('showyspacing', true);
model.result('pg4').set('showsecyspacing', false);
model.result('pg4').set('showsecyextra', false);

model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-delta_T');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg4').run;
model.result('pg3').run;
model.result('pg3').feature('glob1').setIndex('expr', 'Q_wall', 0);
model.result('pg3').run;
model.result('pg4').run;
model.result('pg3').run;
model.result('pg3').set('xlog', false);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', true);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);
model.result('pg3').set('xlog', true);
model.result('pg3').set('showmanualgrid', true);
model.result('pg3').set('showxspacing', false);
model.result('pg3').set('showyspacing', true);
model.result('pg3').set('showsecyspacing', false);
model.result('pg3').set('showsecyextra', false);

model.func.create('rm1', 'Ramp');
model.func.remove('rm1');

model.result('pg3').run;
model.result.create('pg5', 'PlotGroup1D');
model.result('pg5').run;
model.result('pg5').create('lngr1', 'LineGraph');
model.result('pg5').feature('lngr1').selection.set([29]);
model.result('pg5').feature('lngr1').set('expr', 'z');
model.result('pg5').feature('lngr1').set('xdata', 'expr');
model.result('pg5').feature('lngr1').set('xdataunit', 'degC');
model.result('pg5').run;
model.result.dataset.create('cln1', 'CutLine2D');
model.result.dataset('cln1').setIndex('genpoints', '0.5*r_inner', 0, 0);
model.result.dataset('cln1').setIndex('genpoints', '0.5*r_inner', 1, 0);
model.result.dataset('cln1').setIndex('genpoints', '-L_borehole', 1, 1);
model.result.dataset.duplicate('cln2', 'cln1');
model.result.dataset('cln2').setIndex('genpoints', '0.5*(r_outer+r_borehole)', 0, 0);
model.result.dataset('cln2').setIndex('genpoints', '0.5*(r_outer+r_borehole)', 1, 0);
model.result('pg5').run;
model.result('pg5').feature('lngr1').set('data', 'cln1');
model.result('pg5').run;
model.result('pg5').feature.duplicate('lngr2', 'lngr1');
model.result('pg5').run;
model.result('pg5').feature('lngr2').set('data', 'cln2');
model.result('pg5').run;

model.param.set('Q_fluid', '13[L/s]');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg5').run;

model.label('coaxial.mph');

model.study('std1').setGenPlots(false);

out = model;
