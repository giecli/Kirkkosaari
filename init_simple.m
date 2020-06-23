% =========================================================================
% Initializes model with only bedrock.
% =========================================================================

function model = init_simple(params)

import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

model = ModelUtil.create('Model');

model.param.set('R_model', '500[m]');
model.param.set('H_model', '2500[m]');
model.param.set('r_borehole', sprintf('%f[m]', 0.5*params.d_borehole));
model.param.set('L_borehole', sprintf('%f[m]', params.L_borehole));
model.param.set('T_surface', sprintf('%f[degC]', params.T_surface));
model.param.set('q_geothermal', sprintf('%f[W/m^2]', params.q_geothermal));
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
model.param.set('Q_extraction', sprintf('%f[W]', params.Q_extraction));
model.param.set('h_buffer', '1[m]');
model.param.set('A_wall', '2*pi*r_borehole*L_borehole+pi*r_borehole^2');

model.func.create('an1', 'Analytic');
model.func('an1').set('funcname', 'T_initial');
model.func('an1').set('expr', 'T_surface-q_geothermal/k_rock*z');
model.func('an1').set('args', {'z'});
model.func('an1').set('argunit', 'm');
model.func('an1').set('fununit', 'K');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', {'r_borehole+h_buffer' '0';  ...
'R_model' '0';  ...
'R_model' '-H_model';  ...
'0' '-H_model';  ...
'0' '-L_borehole-h_buffer';  ...
'r_borehole+h_buffer' '-L_borehole-h_buffer'});
model.component('comp1').geom('geom1').create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol2').set('table', {'r_borehole' '0';  ...
'r_borehole+h_buffer' '0';  ...
'r_borehole+h_buffer' '-L_borehole-h_buffer';  ...
'0' '-L_borehole-h_buffer';  ...
'0' '-L_borehole';  ...
'r_borehole' '-L_borehole'});
model.component('comp1').geom('geom1').create('topsel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('topsel').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('topsel').label('Ground Surface Selection');
model.component('comp1').geom('geom1').feature('topsel').set('xmin', 'r_borehole-1[mm]');
model.component('comp1').geom('geom1').feature('topsel').set('xmax', 'R_model+1[mm]');
model.component('comp1').geom('geom1').feature('topsel').set('ymin', '-1[mm]');
model.component('comp1').geom('geom1').feature('topsel').set('ymax', '+1[mm]');
model.component('comp1').geom('geom1').feature('topsel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').create('botsel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('botsel').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('botsel').label('Bottom Boundary Selection');
model.component('comp1').geom('geom1').feature('botsel').set('xmin', '-1[mm]');
model.component('comp1').geom('geom1').feature('botsel').set('xmax', 'R_model+1[mm]');
model.component('comp1').geom('geom1').feature('botsel').set('ymin', '-H_model-1[mm]');
model.component('comp1').geom('geom1').feature('botsel').set('ymax', '-H_model+1[mm]');
model.component('comp1').geom('geom1').feature('botsel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').create('wallsel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('wallsel').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('wallsel').label('Borehole Wall Selection');
model.component('comp1').geom('geom1').feature('wallsel').set('xmin', '-1[mm]');
model.component('comp1').geom('geom1').feature('wallsel').set('xmax', 'r_borehole+1[mm]');
model.component('comp1').geom('geom1').feature('wallsel').set('ymin', '-L_borehole-1[mm]');
model.component('comp1').geom('geom1').feature('wallsel').set('ymax', '+1[mm]');
model.component('comp1').geom('geom1').feature('wallsel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run;

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('T_min', 'minop1(T)');
model.component('comp1').variable('var1').set('T_ave', 'aveop1(T)');

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('minop1').selection.named('geom1_wallsel');
model.component('comp1').cpl('aveop1').selection.named('geom1_wallsel');

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').selection.named('geom1_wallsel');
model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp1').selection.named('geom1_topsel');
model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').selection.named('geom1_botsel');

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('edg1').selection.named('geom1_wallsel');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([1]);

model.component('comp1').cpl('aveop1').set('axisym', true);

model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock'; '0'; '0'; '0'; 'k_rock'; '0'; '0'; '0'; 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');
model.component('comp1').physics('ht').feature('hf1').set('q0', '-Q_extraction/A_wall');
model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_surface');
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');

model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '20[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').create('i1', 'Iterative');
model.sol('sol1').feature('t1').create('i2', 'Iterative');
model.sol('sol1').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', '0 50');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'0 50' '1e-6[a]'});
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('tlist', '0 50');
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 5);
model.sol('sol1').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('stabacc', 'aacc');
model.sol('sol1').feature('t1').feature('fc1').set('aaccdim', 5);
model.sol('sol1').feature('t1').feature('fc1').set('aaccmix', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('aaccdelay', 1);
model.sol('sol1').feature('t1').feature('d1').label('Direct, Heat Transfer Variables (ht)');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('i1').label('AMG, Heat Transfer Variables (ht)');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('i2').label('GMG, Heat Transfer Variables (ht)');
model.sol('sol1').feature('t1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
