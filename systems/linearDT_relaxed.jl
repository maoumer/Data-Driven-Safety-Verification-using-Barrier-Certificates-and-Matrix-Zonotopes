# Install these libraries beforehand
using MAT, LazySets, Polyhedra, DynamicPolynomials, LinearAlgebra

# "x0min","x0max","xumin","xumax","xmin","xmax","umin","umax","Acenter","Bcenter","uncertaintyMin","uncertaintyMax"
matlab = matread("vars_5.mat"); # change to respective dimensions
x0min = matlab["x0min"]; #min of initial states
x0max = matlab["x0max"]; #max of initial states
xumin = matlab["xumin"]; #min of unsafe states
xumax = matlab["xumax"]; #max of unsafe states
xmin = matlab["xmin"]; #min of state space
xmax = matlab["xmax"]; #max of state space
umin = matlab["umin"]; #min of input
umax = matlab["umax"]; #max of input
Acenter = matlab["Acenter"];
Bcenter = matlab["Bcenter"];
Dmin = matlab["uncertaintyMin"];
Dmax = matlab["uncertaintyMax"];

m = size(Bcenter)[2]; # inputs dimension
dim = size(xmin)[1]
@polyvar x[1:dim] u d[1:dim]

vars = [x;u;d]

# system dynamics - discrete-time linear system model with noise
function dyn(x,u,d,A,B)
     f = A*x + B*u + d
     return f
end

# safety: semi-algebraic set description, >=0 by default
go = [x.-x0min; x0max.-x] # initial states
gu = [x.-xumin; xumax.-x]  # unsafe states
urange = [u.-umin; umax.-u] # input
drange = [d.-Dmin; Dmax.-d] # full uncertainty outside center model (from generators + additive noise as zonotope)
g = [x.-xmin; xmax.-x; drange] # state space
