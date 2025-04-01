using MAT, LazySets, Polyhedra, DynamicPolynomials, LinearAlgebra

# "Amin","Bmin","Amax","Bmax","x0min","x0max","xumin","xumax","xmin","xmax","umin","umax","wcenter","wmin","wmax"
matlab = matread("./vars_2.mat"); # 2D system
Amin = matlab["Amin"]; #min of set of A
Bmin = matlab["Bmin"]; #min of set of B
Amax = matlab["Amax"]; #max of set of A
Bmax = matlab["Bmax"]; #max of set of B
x0min = matlab["x0min"]; #min of initial states
x0max = matlab["x0max"]; #max of initial states
xumin = matlab["xumin"]; #min of unsafe states
xumax = matlab["xumax"]; #max of unsafe states
xmin = matlab["xmin"]; #min of state space
xmax = matlab["xmax"]; #max of state space
umin = matlab["umin"]; #min of input
umax = matlab["umax"]; #max of input
wcenter = matlab["wcenter"];
wmin = matlab["wmin"]; #min of additive noise
wmax = matlab["wmax"]; #max of additive noise

dim = size(xmin)[1]
@polyvar x[1:dim] u w[1:dim] A[1:dim,1:dim] B[1:dim]

vars = [x;u;w;B;vec(A)]

# system dynamics - discrete-time linear system model with noise
function dyn(x,u,w,A,B)
     f = A*x + B*u + w
     return f
end

# safety: semi-algebraic set description, >=0 by default
g = [x.-xmin; xmax.-x] # state space
go = [x.-x0min; x0max.-x] # initial states
gu = [x.-xumin; xumax.-x]  # unsafe states
urange = [u.-umin; umax.-u] # input
wrange = [w.-wmin; wmax.-w] # additive noise
Arange = [A.-Amin; Amax.-A]
Brange = [B.-Bmin; Bmax.-B]
uwAB = [urange; wrange; Brange; vec(Arange)] 
