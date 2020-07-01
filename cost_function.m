function cost = cost_function(model, x)
model.param.set('Q_extraction', sprintf('%f[W]', x(1)));
model.param.set('Q_fluid', sprintf('%f[L/s]', x(2)));
model.sol('sol1').runAll;
T = mphglobal(model, 'T_min', 'unit', 'degC', 'solnum', 'end');
cost = abs(T(end));
