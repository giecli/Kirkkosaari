function model = init_complex(params)

fprintf(1, 'Entered init_complex.\n');

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

% -------------------------------------------------------------------------
% Creates the initial temperature function.
% -------------------------------------------------------------------------

model.func.create('an1', 'Analytic');
model.func('an1').set('funcname', 'T_initial');
model.func('an1').set('expr', 'T_surface-q_geothermal/k_rock*z');
model.func('an1').set('args', {'z'});
model.func('an1').set('argunit', 'm');
model.func('an1').set('fununit', 'K');

% -------------------------------------------------------------------------
% Creates the geometry and runs it.
% -------------------------------------------------------------------------

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'0' '-H_soil'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'R_model' 'H_soil'});

model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0' '-H_soil-H_clay'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'R_model' 'H_clay'});

model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'0' '-H_model'});
model.component('comp1').geom('geom1').feature('r3').set('size', {'R_model' 'H_model-H_soil-H_clay'});

model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'0' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'r_borehole' 'L_borehole'});

model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1' 'r2' 'r3'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'r4'});

model.component('comp1').geom('geom1').run('fin');

% -------------------------------------------------------------------------
% Creates selections and updates the geometry with them.
% -------------------------------------------------------------------------

model.component('comp1').geom('geom1').create('wall_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('wall_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('wall_selection').label('Borehole Wall Selection');
model.component('comp1').geom('geom1').feature('wall_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('wall_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('wall_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('wall_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('wall_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('surface_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('surface_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('surface_selection').label('Ground Surface Selection');
model.component('comp1').geom('geom1').feature('surface_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('surface_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('surface_selection').set('ymin', 0);
model.component('comp1').geom('geom1').feature('surface_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('surface_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('bottom_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('bottom_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('bottom_selection').label('Bottom Boundary Selection');
model.component('comp1').geom('geom1').feature('bottom_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('bottom_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('bottom_selection').set('ymin', '-H_model');
model.component('comp1').geom('geom1').feature('bottom_selection').set('ymax', '-H_model');
model.component('comp1').geom('geom1').feature('bottom_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('soil_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('soil_selection').label('Soil Domain Selection');
model.component('comp1').geom('geom1').feature('soil_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('soil_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('soil_selection').set('ymin', '-H_soil');
model.component('comp1').geom('geom1').feature('soil_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('soil_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('clay_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('clay_selection').label('Clay Domain Selection');
model.component('comp1').geom('geom1').feature('clay_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('clay_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('clay_selection').set('ymin', '-H_soil-H_clay');
model.component('comp1').geom('geom1').feature('clay_selection').set('ymax', '-H_soil');
model.component('comp1').geom('geom1').feature('clay_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').run;

% -------------------------------------------------------------------------
% Creates variables and component couplings.
% -------------------------------------------------------------------------

model.component('comp1').variable.create('var1');

model.component('comp1').variable('var1').set('T_min', 'minop1(T)');
model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').label('Borehole Wall Minimum Operator');
model.component('comp1').cpl('minop1').selection.named('geom1_wall_selection');

model.component('comp1').variable('var1').set('T_ave', 'aveop1(T)');
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').label('Borehole Wall Average Operator');
model.component('comp1').cpl('aveop1').selection.named('geom1_wall_selection');
model.component('comp1').cpl('aveop1').set('axisym', true);

model.component('comp1').variable('var1').set('A_wall', 'intop1(1)');
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').label('Borehole Wall Integration Operator');
model.component('comp1').cpl('intop1').selection.named('geom1_wall_selection');
model.component('comp1').cpl('intop1').set('axisym', true);

% -------------------------------------------------------------------------
% Creates the physics.
% -------------------------------------------------------------------------

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');

model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock'; '0'; '0'; '0'; 'k_rock'; '0'; '0'; '0'; 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');

model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');

model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').label('Borehole Wall Heat Flux');
model.component('comp1').physics('ht').feature('hf1').selection.named('geom1_wall_selection');
model.component('comp1').physics('ht').feature('hf1').set('q0', '-Q_extraction/A_wall');

model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp1').selection.named('geom1_surface_selection');
model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_surface');

model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').selection.named('geom1_bottom_selection');
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');

model.component('comp1').physics('ht').create('solid2', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid2').selection.named('geom1_soil_selection');

model.component('comp1').physics('ht').feature('solid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('k', {'k_soil'; '0'; '0'; '0'; 'k_soil'; '0'; '0'; '0'; 'k_soil'});
model.component('comp1').physics('ht').feature('solid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('rho', 'rho_soil');
model.component('comp1').physics('ht').feature('solid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('Cp', 'Cp_soil');

model.component('comp1').physics('ht').create('solid3', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid3').selection.named('geom1_clay_selection');

model.component('comp1').physics('ht').feature('solid3').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('k', {'k_clay'; '0'; '0'; '0'; 'k_clay'; '0'; '0'; '0'; 'k_clay'});
model.component('comp1').physics('ht').feature('solid3').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('rho', 'rho_clay');
model.component('comp1').physics('ht').feature('solid3').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('Cp', 'Cp_clay');

% -------------------------------------------------------------------------
% Creates the mesh and runs it.
% -------------------------------------------------------------------------

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.named('geom1_wall_selection');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '20[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);

model.component('comp1').mesh('mesh1').run;

% -------------------------------------------------------------------------
% Creates the study and solution.
% -------------------------------------------------------------------------

tlist = sprintf('0 %d', params.t_simulation);

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', tlist);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {tlist '1e-6[a]'});
model.sol('sol1').feature('t1').set('control', 'user');
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('tlist', tlist);
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);

fprintf(1, 'Exiting init_complex.\n');
