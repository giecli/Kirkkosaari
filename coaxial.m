function out = model
%
% coaxial.m
%
% Model exported on Jun 30 2020, 11:37 by COMSOL 5.5.0.292.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Work\Kirkkosaari');

model.label('coaxial.mph');

model.param.set('R_model', '500.000000[m]');
model.param.set('H_model', '2500.000000[m]');
model.param.set('T_surface', '7.000000[degC]');
model.param.set('q_geothermal', '0.042000[W/m^2]');
model.param.set('L_borehole', '2000.000000[m]');
model.param.set('r_borehole', '0.038000[m]');
model.param.set('r_inner', '0.016000[m]');
model.param.set('r_outer', '0.025000[m]');
model.param.set('k_rock', '3[W/(m*K)]');
model.param.set('Cp_rock', '728.000000[J/(kg*K)]');
model.param.set('rho_rock', '2850.000000[kg/m^3]');
model.param.set('H_soil', '20.000000[m]');
model.param.set('k_soil', '1[W/(m*K)]');
model.param.set('Cp_soil', '728.000000[J/(kg*K)]');
model.param.set('rho_soil', '2850.000000[kg/m^3]');
model.param.set('H_clay', '50.000000[m]');
model.param.set('k_clay', '2[W/(m*K)]');
model.param.set('Cp_clay', '728.000000[J/(kg*K)]');
model.param.set('rho_clay', '2850.000000[kg/m^3]');
model.param.set('k_pipe', '0.100000[W/(m*K)]');
model.param.set('Cp_pipe', '1900.000000[J/(kg*K)]');
model.param.set('rho_pipe', '950.000000[kg/m^3]');
model.param.set('k_fluid', '0.600000[W/(m*K)]');
model.param.set('Cp_fluid', '4186.000000[J/(kg*K)]');
model.param.set('rho_fluid', '1000.000000[kg/m^3]');
model.param.set('Q_fluid', '0.600000[L/s]');
model.param.set('Q_extraction', '70000.000000[W]');
model.param.set('A_inner', 'pi*r_inner^2');
model.param.set('A_outer', 'pi*r_borehole^2-pi*r_outer^2');
model.param.set('v_inner', 'Q_fluid/A_inner');
model.param.set('v_outer', 'Q_fluid/A_outer');
model.param.set('delta_T', 'Q_extraction/(rho_fluid*Cp_fluid*Q_fluid)');
model.param.set('h_buffer', '1.000000[m]');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.func.create('an1', 'Analytic');
model.func.create('gp1', 'GaussianPulse');
model.func.create('an2', 'Analytic');
model.func.create('step1', 'Step');
model.func.create('pw1', 'Piecewise');
model.func('an1').set('funcname', 'T_initial');
model.func('an1').set('expr', 'T_surface-q_geothermal/k_rock*z');
model.func('an1').set('args', {'z'});
model.func('an1').set('argunit', 'm');
model.func('an1').set('fununit', 'K');
model.func('an1').set('plotargs', {'z' '0' '1'});
model.func('gp1').set('funcname', 'gaussian_pulse');
model.func('gp1').set('location', '-0.5*L_borehole');
model.func('gp1').set('sigma', '0.2*L_borehole');
model.func('gp1').set('normalization', 'peak');
model.func('an2').set('funcname', 'mesh_size');
model.func('an2').set('expr', 'max(min(gaussian_pulse(z[1/m]),0.5),0.001)');
model.func('an2').set('args', {'z'});
model.func('an2').set('argunit', 'm');
model.func('an2').set('fununit', 'K');
model.func('an2').set('plotargs', {'z' '-L_borehole' '0'});
model.func('step1').set('funcname', 'ramp_function');
model.func('step1').set('location', '1/365.2425');
model.func('step1').set('smooth', '1/182.62125');
model.func('pw1').set('arg', 'z');
model.func('pw1').set('extrap', 'interior');
model.func('pw1').set('pieces', {'-H_soil' '0' 'T_surface-q_geothermal/k_soil*z'; '-H_soil-H_clay' '-H_soil' 'T_surface+q_geothermal/k_soil*H_soil-q_geothermal/k_clay*(z+H_soil)'; '-H_model' '-H_soil-H_clay' 'T_surface+q_geothermal/k_soil*H_soil+q_geothermal/k_clay*H_clay-q_geothermal/k_rock*(z+H_soil+H_clay)'});
model.func('pw1').set('argunit', 'm');
model.func('pw1').set('fununit', 'K');

model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

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
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('T_min', 'minop1(T)');
model.component('comp1').variable('var1').set('Q_wall', 'intop1(ht.ndflux)');
model.component('comp1').variable('var1').set('T_inlet', 'aveop1(T)');
model.component('comp1').variable('var1').set('T_bottom', 'aveop2(T)');
model.component('comp1').variable('var1').set('T_outlet', 'aveop3(T)');

model.view.create('view2', 3);

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('minop1').selection.geom('geom1', 1);
model.component('comp1').cpl('minop1').selection.set([26 27 29]);
model.component('comp1').cpl('intop1').selection.geom('geom1', 1);
model.component('comp1').cpl('intop1').selection.set([26 27 29]);
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop1').selection.set([25]);
model.component('comp1').cpl('aveop2').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop2').selection.set([20]);
model.component('comp1').cpl('aveop3').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop3').selection.set([11]);

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').selection.set([31 37]);
model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').selection.set([2]);
model.component('comp1').physics('ht').create('solid2', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid2').selection.set([13 15]);
model.component('comp1').physics('ht').create('solid3', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid3').selection.set([12 14]);
model.component('comp1').physics('ht').create('fluid1', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid1').selection.set([3 4 5]);
model.component('comp1').physics('ht').create('solid4', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid4').selection.set([6 7 8]);
model.component('comp1').physics('ht').create('fluid2', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid2').selection.set([9 10 11]);
model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp2').selection.set([25]);
model.component('comp1').physics('ht').create('temp3', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp3').selection.set([6]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri4', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri5', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri6', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([11 18 25]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([5 7 9 12 14 16 19 21 23 26 27 29]);
model.component('comp1').mesh('mesh1').feature('edg2').create('se1', 'SizeExpression');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([3 4 5 6 7 8 9 10 11]);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([13]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([12]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri4').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri4').selection.set([15]);
model.component('comp1').mesh('mesh1').feature('ftri4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri5').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri5').selection.set([14]);
model.component('comp1').mesh('mesh1').feature('ftri5').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri6').create('size1', 'Size');

model.component('comp1').view('view1').axis.set('xmin', 127.1814193725586);
model.component('comp1').view('view1').axis.set('xmax', 1088.78369140625);
model.component('comp1').view('view1').axis.set('ymin', -344.0926513671875);
model.component('comp1').view('view1').axis.set('ymax', 107.66679382324219);

model.component('comp1').cpl('minop1').label('Borehole Wall Minimum Operator');
model.component('comp1').cpl('intop1').label('Borehole Wall Integration Operator');
model.component('comp1').cpl('intop1').set('axisym', true);
model.component('comp1').cpl('aveop1').label('Top Inlet Average Operator');
model.component('comp1').cpl('aveop1').set('axisym', true);
model.component('comp1').cpl('aveop2').label('Bottom Outlet Average Operator');
model.component('comp1').cpl('aveop2').set('axisym', true);
model.component('comp1').cpl('aveop3').label('Top Outlet Average Operator');
model.component('comp1').cpl('aveop3').set('axisym', true);

model.component('comp1').physics('ht').prop('ShapeProperty').set('order_temperature', 1);
model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock'; '0'; '0'; '0'; 'k_rock'; '0'; '0'; '0'; 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');
model.component('comp1').physics('ht').feature('solid1').label('Bedrock Solid');
model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');
model.component('comp1').physics('ht').feature('hf1').set('q0', '-q_geothermal');
model.component('comp1').physics('ht').feature('hf1').label('Surface Heat Flux BC');
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');
model.component('comp1').physics('ht').feature('hf2').label('Bottom Heat Flux BC');
model.component('comp1').physics('ht').feature('solid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('k', {'k_soil'; '0'; '0'; '0'; 'k_soil'; '0'; '0'; '0'; 'k_soil'});
model.component('comp1').physics('ht').feature('solid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('rho', 'rho_soil');
model.component('comp1').physics('ht').feature('solid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('Cp', 'Cp_soil');
model.component('comp1').physics('ht').feature('solid2').label('Soil Solid');
model.component('comp1').physics('ht').feature('solid3').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('k', {'k_clay'; '0'; '0'; '0'; 'k_clay'; '0'; '0'; '0'; 'k_clay'});
model.component('comp1').physics('ht').feature('solid3').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('rho', 'rho_clay');
model.component('comp1').physics('ht').feature('solid3').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('Cp', 'Cp_clay');
model.component('comp1').physics('ht').feature('solid3').label('Claystone Solid');
model.component('comp1').physics('ht').feature('fluid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('k', {'1000'; '0'; '0'; '0'; '1000'; '0'; '0'; '0'; 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('u', {'0'; '0'; 'ramp_function(t[1/a])*v_inner'});
model.component('comp1').physics('ht').feature('fluid1').label('Inner Fluid');
model.component('comp1').physics('ht').feature('solid4').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('k', {'k_pipe'; '0'; '0'; '0'; 'k_pipe'; '0'; '0'; '0'; 'k_pipe'});
model.component('comp1').physics('ht').feature('solid4').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('rho', 'rho_pipe');
model.component('comp1').physics('ht').feature('solid4').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('Cp', 'Cp_pipe');
model.component('comp1').physics('ht').feature('solid4').label('Pipe Solid');
model.component('comp1').physics('ht').feature('fluid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('k', {'1000'; '0'; '0'; '0'; '1000'; '0'; '0'; '0'; 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('u', {'0'; '0'; '-ramp_function(t[1/a])*v_outer'});
model.component('comp1').physics('ht').feature('fluid2').label('Outer Fluid');
model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-ramp_function(t[1/a])*delta_T');
model.component('comp1').physics('ht').feature('temp2').label('Top Inlet Temperature BC');
model.component('comp1').physics('ht').feature('temp3').set('T0', 'T_bottom');
model.component('comp1').physics('ht').feature('temp3').label('Bottom Inlet Temperature BC');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.001);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '2.5[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg2').feature('se1').set('sizeexpr', 'mesh_size(z)');
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri2').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri3').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri4').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmax', 5);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.create('pg1', 'PlotGroup3D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup1D');
model.result.create('pg4', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('con1', 'Contour');
model.result('pg3').create('glob1', 'Global');
model.result('pg4').create('lngr1', 'LineGraph');
model.result('pg4').feature('lngr1').set('xdata', 'expr');
model.result('pg4').feature('lngr1').selection.set([38 39 40]);
model.result('pg4').feature('lngr1').set('expr', 'z');

model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', '0 50');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', 0.001);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'0 50' '1e-6[a]'});
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('tlist', '0 50');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('dDef').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.result.dataset('rev1').label('Revolution 2D');
model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result('pg1').label('Temperature, 3D (ht)');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('colortable', 'ThermalLight');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label('Isothermal Contours (ht)');
model.result('pg2').feature('con1').label('Contour');
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('colortable', 'ThermalLight');
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('resolution', 'normal');
model.result('pg3').set('xlabel', 'Time (a)');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').feature('glob1').set('expr', {'T_min'});
model.result('pg3').feature('glob1').set('unit', {'degC'});
model.result('pg3').feature('glob1').set('descr', {''});
model.result('pg4').set('looplevelinput', {'last'});
model.result('pg4').set('xlabel', 'Temperature (degC)');
model.result('pg4').set('ylabel', 'z-coordinate (m)');
model.result('pg4').set('xlabelactive', false);
model.result('pg4').set('ylabelactive', false);
model.result('pg4').feature('lngr1').set('xdataunit', 'degC');
model.result('pg4').feature('lngr1').set('resolution', 'normal');

out = model;
