function cost = cost_function(model, x)
model.param.set('Q_extraction', sprintf('%f[W]', x));
model.sol('sol1').runAll;
T = mphglobal(model, 'T_ave', 'unit', 'degC');
cost = abs(T(end));
