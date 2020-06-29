from pylab import *

H_soil = 20
H_clay = 50
H_model = 2500

# Soil layer...

z = linspace(0, -H_soil, 21)

expr = (0.05-z/H_soil*2.45)

print(expr)

# Soil layer...

z = linspace(-H_soil, -H_soil-H_clay, 21)

expr = 2.5+22.5*(-z-H_soil)/H_clay

print(expr)
