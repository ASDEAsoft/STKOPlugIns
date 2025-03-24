
# frictionModel Coulomb
frictionModel Coulomb 111 1.1

# frictionModel VelDepMultiLinear
frictionModel VelDepMultiLinear 211 -vel 0.1 0.2 -frn 2.1 2.2

# frictionModel VelDependent
frictionModel VelDependent 311 0.22 0.33 0.44

# frictionModel VelNormalFrcDep
frictionModel VelNormalFrcDep 411 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08

# frictionModel VelPressureDep
frictionModel VelPressureDep 511 1.1 1.2 1.3 1.4 1.5 1.6

uniaxialMaterial Elastic 10 100.1 
uniaxialMaterial Elastic 20 200.2 
uniaxialMaterial Elastic 30 300.3 
uniaxialMaterial Elastic 40 400.4 
uniaxialMaterial Elastic 50 400.5 

node 1 0 0 0
node 2 0 0 100
node 3 100 0 0
node 4 100 0 100

element TripleFrictionPendulum 123 1 2 111 311 511 10 20 30 40 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1
element TripleFrictionPendulum 555 3 4 211 411 511 10 20 50 40 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1