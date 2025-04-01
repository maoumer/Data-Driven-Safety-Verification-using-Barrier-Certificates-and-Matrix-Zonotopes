# hard = run the SOS optimization using the harder condition

# sos_zonotope_hard - searches for barrier certificate for data driven set of models for discrete time systems
# x(k+1) = Ax(k) + Bu(k) + w(k) using hard condition
#
# Other jl-files required: linearDT_hard
# Subfunctions: none
# MAT-files required: vars_2.mat
# 
# Author:       Mohammed Adib Oumer
# Written:      February 4, 2025
# Last update:  March 21, 2025
# Last revision:---

# include important libraries - install beforehand. Additional libraries in other JL file
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using LinearAlgebra
using TSSOS # important for SOS, https://github.com/wangjie212/TSSOS

ϵ = 10^(-5) # B(x) ≥ ϵ instead of B(x) > 0
sos_tol = 1 # the maximum degree of unknown SOS polynomials = deg + sos_tol 
error = 5   # precision digit places

function bc_hard(deg)
    # synthesize BC by using the standard formulation
    # deg: degree of BC template
    model = Model(optimizer_with_attributes(Mosek.Optimizer)) #COSMO.Optimizer
    set_optimizer_attribute(model, MOI.Silent(), true) # uncomment for verbose output
    
    f = dyn(x,u,w,A,B);
    BC, Bc, Bb = add_poly!(model, x, deg) # generate polynomial template with given variables and degree
    Bf = BC(x => f) #B(f(x)) https://juliapackages.com/p/dynamicpolynomials

    # 1st condition, >= 0 by default, -B >= 0
    model,_ = add_psatz!(model, -BC, x, go, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    # 2nd condition
    model,_ = add_psatz!(model, BC-ϵ , x, gu, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    # 3rd and last condition
    model,_ = add_psatz!(model, BC-Bf, vars, [g;uwAB], [], div(maxdegree(BC-Bf)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    
    optimize!(model) # solve for coefficients
    status = termination_status(model)
    Bc = value.(Bc)  # get the values of each coefficient
    for i in eachindex(Bc)
        Bc[i] = round(Bc[i]; digits = error) # round to order of error
    end

    # objv = objective_value(model)
    # dual = dual_objective_value(model)
    # if(abs(objv-dual)<=1e-4 && string(status) == "SLOW_PROGRESS") # MOSEK specific
    #     status = "optimal"
    # end

    # status might be optimal but if all Bc approx 10^{-error}, it's essentially 0.
    return (status,Bc'*Bb,Bc,Bb) # optimization status and barrier certificate function
end



# Simulation
names = ["linearDT_hard"]
system_data = ["vars: ","f: ","go: ","gu: ","g: ","urange: ","wrange: ","Brange: ","Arange: "]
max_deg = 3
for name in names
    include("./systems/"*name*".jl"); # load system dynamics info
    f = dyn(x,u,w,A,B)
    file = open("./systems/"*name*"_system.txt", "w"); # open file to write/save system dynamics info
    for (j,k) in zip(system_data, [vars,f,go,gu,g,urange,wrange,Brange,Arange]) # save these data in readable form
        write(file, j*"{")
        for i = 1:length(k)-1
            write(file, string(k[i])*", ")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)

    # print sufficient condition results for standard barrier certificate
    file = open("./systems/"*name*"_bc_dim"*string(dim)*".txt", "w");
    for deg = 1:max_deg
        stats = @timed data = bc_hard(deg) # time execution
        status, BC, Bcoef, Bmonoms = data
        write(file, "poly deg: "*string(deg)*"\n")
        write(file, "status: "*string(status)*"\n")
        write(file, "B: "*string(BC)*"\n")
        write(file, "Coefficients: "*Base.replace(string(Bcoef),"e"=>"*10^")*"\n")
        write(file, "Monomials: "*string(Bmonoms)*"\n") 
        write(file, "time: "*string(stats.time)*"\n\n") 
    end
    close(file)
    println("Finished "*name)
end