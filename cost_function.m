function cost = cost_function(x)
global model
model.param.set('Q_extraction', sprintf('%f[W]', x));
model.sol('sol1').runAll;
T_min = mphglobal(model, 'T_min', 'unit', 'degC');
cost = abs(T_min(end));
