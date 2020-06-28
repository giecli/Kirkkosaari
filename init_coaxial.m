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

model.component('comp1').curvedInterior(false);

model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

% -------------------------------------------------------------------------
% Sets up model parameters.
% -------------------------------------------------------------------------

model.param.set('R_model', '500[m]');
model.param.set('H_model', '2500[m]');

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

model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'0' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'r_borehole' 'L_borehole'});

model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1' 'r2' 'r3'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'r4'});

model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').set('pos', {'0' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r5').set('size', {'r_inner' 'L_borehole'});

model.component('comp1').geom('geom1').create('r6', 'Rectangle');
model.component('comp1').geom('geom1').feature('r6').set('pos', {'r_inner' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r6').set('size', {'r_outer-r_inner' 'L_borehole'});

model.component('comp1').geom('geom1').create('r7', 'Rectangle');
model.component('comp1').geom('geom1').feature('r7').set('pos', {'r_outer' '-L_borehole'});
model.component('comp1').geom('geom1').feature('r7').set('size', {'r_borehole-r_outer' 'L_borehole'});

model.component('comp1').geom('geom1').run('fin');

% -------------------------------------------------------------------------
% Creates selections and updates the geometry with them.
% -------------------------------------------------------------------------

model.component('comp1').geom('geom1').create('soil_domain_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('soil_domain_selection').label('Soil Domain Selection');
model.component('comp1').geom('geom1').feature('soil_domain_selection').set('xmin', 'r_borehole');
model.component('comp1').geom('geom1').feature('soil_domain_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('soil_domain_selection').set('ymin', '-H_soil');
model.component('comp1').geom('geom1').feature('soil_domain_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('soil_domain_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('clay_domain_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('clay_domain_selection').label('Claystone Domain Selection');
model.component('comp1').geom('geom1').feature('clay_domain_selection').set('xmin', 'r_borehole');
model.component('comp1').geom('geom1').feature('clay_domain_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('clay_domain_selection').set('ymin', '-H_soil-H_clay');
model.component('comp1').geom('geom1').feature('clay_domain_selection').set('ymax', '-H_soil');
model.component('comp1').geom('geom1').feature('clay_domain_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('ground_surface_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ground_surface_selection').label('Ground Surface Selection');
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('xmin', 'r_borehole');
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('ymin', 0);
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('ground_surface_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('bottom_boundary_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').label('Bottom Boundary Selection');
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('xmax', 'R_model');
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('ymin', '-H_model');
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('ymax', '-H_model');
model.component('comp1').geom('geom1').feature('bottom_boundary_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('inner_fluid_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').label('Inner Fluid Selection');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').set('xmin', '0');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').set('xmax', 'r_inner');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').set('ymax', '0');
model.component('comp1').geom('geom1').feature('inner_fluid_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('pipe_domain_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').label('Pipe Domain Selection');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').set('xmin', 'r_inner');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').set('xmax', 'r_outer');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').set('ymax', '0');
model.component('comp1').geom('geom1').feature('pipe_domain_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('outer_fluid_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').label('Outer Fluid Selection');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').set('xmin', 'r_outer');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').set('ymax', '0');
model.component('comp1').geom('geom1').feature('outer_fluid_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('borehole_wall_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('borehole_wall_selection').label('Borehole Wall Selection');
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('xmin', 'r_borehole');
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('borehole_wall_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('top_inlet_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('top_inlet_selection').label('Top Inlet Selection');
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('xmin', 'r_outer');
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('ymin', 0);
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('top_inlet_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('bottom_outlet_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').label('Bottom Outlet Selection');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('xmin', 'r_outer');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('ymax', '-L_borehole');
model.component('comp1').geom('geom1').feature('bottom_outlet_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('bottom_inlet_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').label('Bottom Inlet Selection');
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('xmax', 'r_inner');
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('ymax', '-L_borehole');
model.component('comp1').geom('geom1').feature('bottom_inlet_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('top_outlet_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('top_outlet_selection').label('Top Outlet Selection');
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('xmax', 'r_inner');
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('ymin', 0);
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('top_outlet_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('horizontal_edges_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').label('Horizontal Edges Selection');
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('xmax', 'r_borehole');
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('ymin', 0);
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('horizontal_edges_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('vertical_edge1_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').label('Vertical Edge #1 Selection');
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('xmin', 0);
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('xmax', 0);
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('vertical_edge1_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('vertical_edge2_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').label('Vertical Edge #2 Selection');
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('xmin', 'r_inner');
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('xmax', 'r_inner');
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('vertical_edge2_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').create('vertical_edge3_selection', 'BoxSelection');
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('entitydim', 1);
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').label('Vertical Edge #3 Selection');
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('xmin', 'r_outer');
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('xmax', 'r_outer');
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('ymin', '-L_borehole');
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('ymax', 0);
model.component('comp1').geom('geom1').feature('vertical_edge3_selection').set('condition', 'allvertices');

model.component('comp1').geom('geom1').run;

% -------------------------------------------------------------------------
% Creates component couplings and variables.
% -------------------------------------------------------------------------

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').label('Borehole Wall Minimum Operator');
model.component('comp1').cpl('minop1').selection.named('geom1_borehole_wall_selection');

model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').label('Top Inlet Average Operator');
model.component('comp1').cpl('aveop1').selection.named('geom1_top_inlet_selection');
model.component('comp1').cpl('aveop1').set('axisym', true);

model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl('aveop2').label('Bottom Outlet Average Operator');
model.component('comp1').cpl('aveop2').selection.named('geom1_bottom_outlet_selection');
model.component('comp1').cpl('aveop2').set('axisym', true);

model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('aveop3').label('Top Outlet Average Operator');
model.component('comp1').cpl('aveop3').selection.named('geom1_top_outlet_selection');
model.component('comp1').cpl('aveop3').set('axisym', true);

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').label('Borehole Wall Integration Operator');
model.component('comp1').cpl('intop1').selection.named('geom1_borehole_wall_selection');
model.component('comp1').cpl('intop1').set('axisym', true);

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('T_inlet', 'aveop1(T)');
model.component('comp1').variable('var1').set('T_bottom', 'aveop2(T)');
model.component('comp1').variable('var1').set('T_outlet', 'aveop3(T)');
model.component('comp1').variable('var1').set('Q_wall', 'intop1(-ht.ndflux)');

% -------------------------------------------------------------------------
% Creates the physics.
% -------------------------------------------------------------------------

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');

model.component('comp1').physics('ht').prop('ShapeProperty').set('order_temperature', 1);

model.component('comp1').physics('ht').feature('solid1').label('Bedrock Solid');
model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('k', {'k_rock'; '0'; '0'; '0'; 'k_rock'; '0'; '0'; '0'; 'k_rock'});
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho_rock');
model.component('comp1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('Cp', 'Cp_rock');

model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T_initial(z)');

model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp1').label('Ground Surface Temperature');
model.component('comp1').physics('ht').feature('temp1').selection.named('geom1_ground_surface_selection');
model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_surface');

model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf2').label('Geothermal Heat Flux');
model.component('comp1').physics('ht').feature('hf2').selection.named('geom1_bottom_boundary_selection');
model.component('comp1').physics('ht').feature('hf2').set('q0', 'q_geothermal');

model.component('comp1').physics('ht').create('solid2', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid2').label('Soil Solid');
model.component('comp1').physics('ht').feature('solid2').selection.named('geom1_soil_domain_selection');

model.component('comp1').physics('ht').feature('solid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('k', {'k_soil'; '0'; '0'; '0'; 'k_soil'; '0'; '0'; '0'; 'k_soil'});
model.component('comp1').physics('ht').feature('solid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('rho', 'rho_soil');
model.component('comp1').physics('ht').feature('solid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid2').set('Cp', 'Cp_soil');

model.component('comp1').physics('ht').create('solid3', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid3').label('Claystone Solid');
model.component('comp1').physics('ht').feature('solid3').selection.named('geom1_clay_domain_selection');

model.component('comp1').physics('ht').feature('solid3').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('k', {'k_clay'; '0'; '0'; '0'; 'k_clay'; '0'; '0'; '0'; 'k_clay'});
model.component('comp1').physics('ht').feature('solid3').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('rho', 'rho_clay');
model.component('comp1').physics('ht').feature('solid3').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid3').set('Cp', 'Cp_clay');

model.component('comp1').physics('ht').create('fluid1', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid1').label('Inner Fluid');
model.component('comp1').physics('ht').feature('fluid1').selection.named('geom1_inner_fluid_selection');
model.component('comp1').physics('ht').feature('fluid1').set('u', {'0'; '0'; 'v_inner'});

model.component('comp1').physics('ht').feature('fluid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('k', {'k_fluid' '0' '0' '0' 'k_fluid' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid1').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid1').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid1').set('gamma', '1');

model.component('comp1').physics('ht').create('fluid2', 'FluidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('fluid2').label('Outer Fluid');
model.component('comp1').physics('ht').feature('fluid2').selection.named('geom1_outer_fluid_selection');
model.component('comp1').physics('ht').feature('fluid2').set('u', {'0'; '0'; '-v_outer'});

model.component('comp1').physics('ht').feature('fluid2').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('k', {'k_fluid' '0' '0' '0' 'k_fluid' '0' '0' '0' 'k_fluid'});
model.component('comp1').physics('ht').feature('fluid2').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('rho', 'rho_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('Cp', 'Cp_fluid');
model.component('comp1').physics('ht').feature('fluid2').set('gamma_mat', 'userdef');
model.component('comp1').physics('ht').feature('fluid2').set('gamma', '1');

model.component('comp1').physics('ht').create('solid4', 'SolidHeatTransferModel', 2);
model.component('comp1').physics('ht').feature('solid4').label('Pipe Solid');
model.component('comp1').physics('ht').feature('solid4').selection.named('geom1_pipe_domain_selection');

model.component('comp1').physics('ht').feature('solid4').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('k', {'k_pipe'; '0'; '0'; '0'; 'k_pipe'; '0'; '0'; '0'; 'k_pipe'});
model.component('comp1').physics('ht').feature('solid4').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('rho', 'rho_pipe');
model.component('comp1').physics('ht').feature('solid4').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('Cp', 'Cp_pipe');

model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp2').label('Top Inlet Temperature');
model.component('comp1').physics('ht').feature('temp2').selection.named('geom1_top_inlet_selection');
model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_outlet-delta_T');

model.component('comp1').physics('ht').create('temp3', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp3').label('Bottom Inlet Temperature');
model.component('comp1').physics('ht').feature('temp3').selection.named('geom1_bottom_inlet_selection');
model.component('comp1').physics('ht').feature('temp3').set('T0', 'T_bottom');

% -------------------------------------------------------------------------
% Creates the mesh and runs it.
% -------------------------------------------------------------------------

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.named('geom1_horizontal_edges_selection');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', '5[mm]');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg2').selection.named('geom1_borehole_wall_selection');
model.component('comp1').mesh('mesh1').feature('edg2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hmax', '5[mm]');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').create('edg3', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg3').selection.named('geom1_vertical_edge1_selection');
model.component('comp1').mesh('mesh1').feature('edg3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmax', '5[mm]');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').create('edg4', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg4').selection.named('geom1_vertical_edge2_selection');
model.component('comp1').mesh('mesh1').feature('edg4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('hmax', '5[mm]');
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').create('edg5', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg5').selection.named('geom1_vertical_edge3_selection');
model.component('comp1').mesh('mesh1').feature('edg5').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg5').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg5').feature('size1').set('hmax', '5[mm]');
model.component('comp1').mesh('mesh1').feature('edg5').feature('size1').set('hmaxactive', true);
return
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.named('geom1_inner_fluid_selection');

return

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');

model.result.table('tbl1').comments('Line Integration 1');

model.capeopen.label('Thermodynamics Package');

model.component('comp1').variable('var1').label('Component Variables');

model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_outlet-delta_T');
model.component('comp1').physics('ht').feature('temp1').label('Inlet Temperature BC');
model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_bottom');
model.component('comp1').physics('ht').feature('temp2').label('Bottom Temperature BC');
model.component('comp1').physics('ht').feature('hf1').set('q0', 'q_geothermal');
model.component('comp1').physics('ht').feature('hf1').label('Bottom Heat Flux BC');
model.component('comp1').physics('ht').feature('hf2').set('q0', '-q_geothermal');
model.component('comp1').physics('ht').feature('hf2').label('Top Heat Flux BC');

model.component('comp1').mesh('mesh1').label('Mesh');
model.component('comp1').mesh('mesh1').feature('ftri1').label('Bedrock Mesh');
model.component('comp1').mesh('mesh1').feature('ftri1').set('method', 'del');
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
model.sol('sol1').feature('t1').feature.remove('dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.result.dataset.create('cln1', 'CutLine2D');
model.result.dataset.create('cln2', 'CutLine2D');
model.result.numerical.create('int1', 'IntLine');
model.result.numerical('int1').selection.set([12]);
model.result.numerical('int1').set('probetag', 'none');
model.result.create('pg5', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup1D');
model.result.create('pg4', 'PlotGroup1D');
model.result.create('pg6', 'PlotGroup2D');
model.result.create('pg7', 'PlotGroup1D');
model.result.create('pg8', 'PlotGroup1D');
model.result('pg5').create('lngr1', 'LineGraph');
model.result('pg5').create('lngr2', 'LineGraph');
model.result('pg5').feature('lngr1').set('xdata', 'expr');
model.result('pg5').feature('lngr1').set('expr', 'z');
model.result('pg5').feature('lngr2').set('xdata', 'expr');
model.result('pg5').feature('lngr2').set('expr', 'z');
model.result('pg3').create('glob1', 'Global');
model.result('pg4').create('glob1', 'Global');
model.result('pg6').create('con1', 'Contour');
model.result('pg6').create('surf1', 'Surface');
model.result('pg6').feature('con1').set('expr', 'T-T_initial(z)');
model.result('pg6').feature('surf1').set('expr', 'T-T_initial(z)');
model.result('pg7').create('glob1', 'Global');
model.result('pg8').create('lngr1', 'LineGraph');
model.result('pg8').feature('lngr1').set('xdata', 'expr');
model.result('pg8').feature('lngr1').selection.set([12]);
model.result('pg8').feature('lngr1').set('expr', 'z');
model.result.export.create('plot1', 'Plot');
model.result.export.create('img1', 'Image');

model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', '0 50');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', '1e-3');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'0 10^range(-3,0.1,2)' '1e-6[a]'});
model.sol('sol1').feature('t1').set('control', 'user');
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('tlist', '0 10^range(-3,0.1,2)');
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-6');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('tout', 'tsteps');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 5);
model.sol('sol1').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('stabacc', 'aacc');
model.sol('sol1').feature('t1').feature('fc1').set('aaccdim', 5);
model.sol('sol1').feature('t1').feature('fc1').set('aaccmix', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('aaccdelay', 1);
model.sol('sol1').feature('t1').feature('d1').label('PARDISO (ht)');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.result.dataset('cln1').set('genpoints', {'0.5*r_inner' '0'; '0.5*r_inner' '-L_borehole'});
model.result.dataset('cln2').set('genpoints', {'0.5*(r_outer+r_borehole)' '0'; '0.5*(r_outer+r_borehole)' '-L_borehole'});
model.result.dataset('cln2').set('spacevars', {'cln1x'});
model.result.dataset('cln2').set('normal', {'cln1nx' 'cln1ny'});
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').set('expr', {'ht.ndflux'});
model.result.numerical('int1').set('unit', {'W'});
model.result.numerical('int1').set('descr', {'Normal conductive heat flux'});
model.result.numerical('int1').set('intsurface', true);
model.result.numerical('int1').setResult;
model.result('pg5').label('Fluid Temperature Profiles');
model.result('pg5').set('data', 'none');
model.result('pg5').set('titletype', 'manual');
model.result('pg5').set('title', '76mm borehole, 50mm pipe: thickness 10mm k=0.05[W/(m*K)]');
model.result('pg5').set('xlabel', 'Temperature (degC)');
model.result('pg5').set('ylabel', 'z-coordinate (m)');
model.result('pg5').set('xlabelactive', false);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').feature('lngr1').set('data', 'cln1');
model.result('pg5').feature('lngr1').set('looplevelinput', {'manual'});
model.result('pg5').feature('lngr1').set('looplevel', [1]);
model.result('pg5').feature('lngr1').set('xdataunit', 'degC');
model.result('pg5').feature('lngr1').set('linestyle', 'cycle');
model.result('pg5').feature('lngr1').set('linecolor', 'red');
model.result('pg5').feature('lngr1').set('legend', true);
model.result('pg5').feature('lngr1').set('smooth', 'internal');
model.result('pg5').feature('lngr1').set('resolution', 'normal');
model.result('pg5').feature('lngr2').set('data', 'cln2');
model.result('pg5').feature('lngr2').set('looplevelinput', {'manual'});
model.result('pg5').feature('lngr2').set('looplevel', [1]);
model.result('pg5').feature('lngr2').set('xdataunit', 'degC');
model.result('pg5').feature('lngr2').set('linestyle', 'cyclereset');
model.result('pg5').feature('lngr2').set('linecolor', 'blue');
model.result('pg5').feature('lngr2').set('smooth', 'internal');
model.result('pg5').feature('lngr2').set('resolution', 'normal');
model.result('pg3').label('Heat Rate');
model.result('pg3').set('xlabel', 'Time (a)');
model.result('pg3').set('axislimits', true);
model.result('pg3').set('xmin', 7.979083569138569E-8);
model.result('pg3').set('xmax', 125.29543017236963);
model.result('pg3').set('ymin', 65235.681457469276);
model.result('pg3').set('ymax', 152078.58585926064);
model.result('pg3').set('xlog', true);
model.result('pg3').set('xlabelactive', false);
model.result('pg3').feature('glob1').set('expr', {'Q_wall'});
model.result('pg3').feature('glob1').set('unit', {'W'});
model.result('pg3').feature('glob1').set('descr', {''});
model.result('pg4').label('Temperature Difference');
model.result('pg4').set('xlabel', 'Time (a)');
model.result('pg4').set('axislimits', true);
model.result('pg4').set('xmin', -1.088155727175178);
model.result('pg4').set('xmax', 101.08815572717518);
model.result('pg4').set('ymin', 19.927556159440023);
model.result('pg4').set('ymax', 19.927556159440265);
model.result('pg4').set('xlabelactive', false);
model.result('pg4').feature('glob1').set('expr', {'T_outlet-T_inlet'});
model.result('pg4').feature('glob1').set('unit', {'K'});
model.result('pg4').feature('glob1').set('descr', {''});
model.result('pg6').label('Temperature Disturbance');
model.result('pg6').set('looplevel', [1]);
model.result('pg6').feature('con1').set('levelmethod', 'levels');
model.result('pg6').feature('con1').set('levels', '-0.1 -0.01 -0.001');
model.result('pg6').feature('con1').set('contourlabels', true);
model.result('pg6').feature('con1').set('coloring', 'uniform');
model.result('pg6').feature('con1').set('resolution', 'normal');
model.result('pg6').feature('surf1').set('rangecoloractive', true);
model.result('pg6').feature('surf1').set('rangecolormin', -0.5);
model.result('pg6').feature('surf1').set('rangecolormax', 0.5);
model.result('pg6').feature('surf1').set('colortable', 'WaveLight');
model.result('pg6').feature('surf1').set('resolution', 'normal');
model.result('pg7').set('xlabel', 'Time (a)');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').feature('glob1').set('expr', {'T_outlet'});
model.result('pg7').feature('glob1').set('unit', {'degC'});
model.result('pg7').feature('glob1').set('descr', {''});
model.result('pg8').set('looplevelinput', {'last'});
model.result('pg8').set('xlabel', 'Temperature (degC)');
model.result('pg8').set('ylabel', 'z-coordinate (m)');
model.result('pg8').set('xlabelactive', false);
model.result('pg8').set('ylabelactive', false);
model.result('pg8').feature('lngr1').set('xdataunit', 'degC');
model.result('pg8').feature('lngr1').set('resolution', 'normal');
model.result.export('plot1').set('plotgroup', 'pg7');
model.result.export('plot1').set('plot', 'glob1');
model.result.export('plot1').set('filename', 'C:\Users\kpiippon\Documents\projektit\Kimmo\Outlet_T_coaxial2.txt');
model.result.export('img1').set('sourceobject', 'pg5');
model.result.export('img1').set('size', 'manualweb');
model.result.export('img1').set('unit', 'px');
model.result.export('img1').set('lockratio', 'off');
model.result.export('img1').set('width', '2400');
model.result.export('img1').set('height', '1800');
model.result.export('img1').set('resolution', '300');
model.result.export('img1').set('zoomextents', 'off');
model.result.export('img1').set('antialias', 'on');
model.result.export('img1').set('options1d', 'on');
model.result.export('img1').set('options2d', 'on');
model.result.export('img1').set('options3d', 'off');
model.result.export('img1').set('title1d', 'on');
model.result.export('img1').set('title2d', 'on');
model.result.export('img1').set('title3d', 'on');
model.result.export('img1').set('legend1d', 'on');
model.result.export('img1').set('legend2d', 'on');
model.result.export('img1').set('legend3d', 'on');
model.result.export('img1').set('axes1d', 'on');
model.result.export('img1').set('axes2d', 'on');
model.result.export('img1').set('logo1d', 'on');
model.result.export('img1').set('logo2d', 'on');
model.result.export('img1').set('logo3d', 'on');
model.result.export('img1').set('showgrid', 'on');
model.result.export('img1').set('axisorientation', 'on');
model.result.export('img1').set('grid', 'on');
model.result.export('img1').set('fontsize', '9');
model.result.export('img1').set('background', 'color');
model.result.export('img1').set('gltfincludelines', 'on');
model.result.export('img1').set('qualitylevel', '92');
model.result.export('img1').set('qualityactive', 'off');
model.result.export('img1').set('imagetype', 'png');
model.result.export('img1').set('target', 'file');
model.result.export('img1').set('addsuffix', 'off');
model.result.export('img1').set('lockview', 'off');
model.result.export('img1').set('customcolor', [1 1 1]);
model.result.export('img1').set('width', 1000);
model.result.export('img1').set('height', 600);
model.result.export('img1').set('resolution', 96);
model.result.export('img1').set('sizedesc', '265 x 159 mm');
model.result.export('img1').set('pngfilename', '\\gtknas01\gtk_project\KLDATA$\Energiakaivos\COMSOL-malleja\Koaksiaali-U-putki-vertailuja\coaxial_v3_T_profiles.png');
model.result.export('img1').set('logo1d', false);
model.result.export('img1').set('logo2d', false);
model.result.export('img1').set('zoomextents', true);
model.result.export('img1').set('title1d', false);
model.result.export('img1').set('logo3d', false);
model.result.export('img1').set('size', 'manualweb');
model.result.export('img1').set('unit', 'px');
model.result.export('img1').set('height', '600');
model.result.export('img1').set('width', '1000');
model.result.export('img1').set('lockratio', 'off');
model.result.export('img1').set('resolution', '96');
model.result.export('img1').set('antialias', 'on');
model.result.export('img1').set('zoomextents', 'on');
model.result.export('img1').set('fontsize', '9');
model.result.export('img1').set('customcolor', [1 1 1]);
model.result.export('img1').set('background', 'color');
model.result.export('img1').set('gltfincludelines', 'on');
model.result.export('img1').set('title1d', 'off');
model.result.export('img1').set('legend1d', 'on');
model.result.export('img1').set('logo1d', 'off');
model.result.export('img1').set('options1d', 'on');
model.result.export('img1').set('title2d', 'on');
model.result.export('img1').set('legend2d', 'on');
model.result.export('img1').set('logo2d', 'off');
model.result.export('img1').set('options2d', 'on');
model.result.export('img1').set('title3d', 'on');
model.result.export('img1').set('legend3d', 'on');
model.result.export('img1').set('logo3d', 'off');
model.result.export('img1').set('options3d', 'off');
model.result.export('img1').set('axisorientation', 'on');
model.result.export('img1').set('grid', 'on');
model.result.export('img1').set('axes1d', 'on');
model.result.export('img1').set('axes2d', 'on');
model.result.export('img1').set('showgrid', 'on');
model.result.export('img1').set('target', 'file');
model.result.export('img1').set('qualitylevel', '92');
model.result.export('img1').set('qualityactive', 'off');
model.result.export('img1').set('imagetype', 'png');
model.result.export('img1').set('lockview', 'off');

out = model;
