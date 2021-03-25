
using JuMP, Gurobi, Distributions, Dates;

function loadAndSolveProblem(p,K,N,solver)

    #######################################################
    #
    # El problema que estamos considerando es la asignacion de aviones a rutas para
    # maximizar el beneficio obtenido de la programacion de rutas
    # 
    # Basado en el articulo "A Langrangian  Relaxation Approach to Assigning Aircraft to Routes in Hub
    # and Spoke Networks
    #
    #######################################################
    
    # input:
    # solver : CPLEX, GUROBI
    # pij: beneficion debido a la asignacion de la ruta i al avion j
    # K: Conjuntos de periodos de tiempos
    # Njk: conjuntos de rutas i que podrian utilizar el avion j durante el periodo k si es asignado al avion j
    
    # output: 
    # Status: Optima, Infactible, No Acotada
    # Valor Objetivo:
    # Solucion X:
    
    # definimos el nombre del objeto
    prob = Model(solver=solver);

    # calculamos el numero de rutas y aviones
    num_route, num_aircraft = size(p)

    # Xij: 1 si la ruta i es asignada al avion j
    @variable(prob, x[1:num_route,1:num_aircraft] >=0, Bin)
  
    # funcion objetivo: MAX
    @objective(prob, :Max, sum(p[i,j]*x[i,j] for i=1:num_route,j=1:num_aircraft)); 
    
    # res 1
    # Asegura que cada ruta i es asignada a lo menos un avion
    for i in 1:num_route
        @constraint(prob, sum(x[i,j] for j in 1:num_aircraft) <= 1);
    end
    
    # res 2
    # Asegura que cada avion j es asignada a lo mas una ruta durante un periodo de tiempo k 
    for j in 1:num_aircraft
        for k in collect(keys(K))
            @constraint(prob, sum(x[i,j] for i in 1:num_route if in(i,N[j,k])) <= 1);
        end
    end

    # Resolvemos el problema
    status = solve(prob)
                    
    #print(prob)

    return status, getobjectivevalue(prob), getvalue(x)

end

@everywhere function loadAndSolveSubProblem(p,M,β,K,i,solver)
    
    # input:
    # p array dim 1 x num_aircraft: beneficio de un avion j sobre una ruta dada
    # M array dim num_route x num_aircraft: Conjuntos de periodos de tiempo k durante el cual el avion j podria estar ocupado si este fuera asignado a volar la ruta i
    # β dual array dim num_aircraft x num_periods: Variables duales de la restriccion (res 2)
    # K array dim 1 x num_periods: Conjuntos de periodos de tiempos
    # solver: CPLEX, GUROBI
    
    # output:
    # Status
    # dual
    # objetive
    # X solution
    
    # definimos el nombre del objeto
    prob = Model(solver=solver);

    # calculamos el numero de aviones
    num_aircraft = length(p)

    # Es 1 si se asigna el avion j, sino 0
    @variable(prob, x[1:num_aircraft] >=0)
        
    # Funcion objetivo
    @objective(prob, :Max, sum((p[j]-sum(β[j,k] for k in K if in(k,M[j,i])))*x[j] for j=1:num_aircraft)); 
        
    # res 1
    # Asegura que cada ruta i es asignada a lo menos un avion
    @constraint(prob, res,sum(x[j] for j in 1:num_aircraft) <= 1);

    # resolvemos
    status = solve(prob)
                    
    #print(prob)

    return status, getdual(res), getobjectivevalue(prob), getvalue(x)

end

function genData(Avg_route,num_route,num_aircraft)

    #####################################################
    #
    # Esta rutina se usa para generar la data
    #
    #####################################################
    
    # input:
    # Avg_Route: mean route duration (6,8,10) hours
    # num_route: number of routes (50, 75, 100)
    # num_aircraft: number of aircraft.
    
    # output:
    # profit
    # departure
    # arrival
    # K   : set of time periods
    # Njk : set of routes i that would utilize aircraft j during period k if assigned to aircraft j
    # M   : set of time periods k during which aircraft j would be busy if it were assigned to fly route i
    
    d = DiscreteUniform(0,48) # departure time for each route between 0 and 48 hours. sec 3 
    departure = [floor(Dates.DateTime(Dates.now()+Dates.Hour(rand(d,1)[1])+Dates.Minute(15)), Dates.Minute(15)) for i=1:num_route]

    departure = copy(sort(departure)) # sort
    
    # distribution of route duration
    d = DiscreteUniform(0.5*Avg_route,1.5*Avg_route)

    # hours
    aux = [rand(d,1)[1] for i=1:num_route] 
    arrival = [floor(Dates.DateTime(departure[i]+Dates.Hour(aux[i])), Dates.Minute(15)) for i=1:num_route]
    
    # profit per route hour
    d = DiscreteUniform(10,30) 
    profit = [rand(d,1)[1] for i=1:num_route, j=1:num_aircraft]
    
    K = Dict{Int64,Array{Any}}()
    
    K[1] = [departure[1]]
    
    Ea  = 1;
    add = 1;
    for r = 2:num_route
        if departure[r] >= arrival[Ea]
            push!(K[add],departure[r])
            K[add+1] = [K[add][2]];
            add = add + 1;
            Ea = r;
        end
    end
    
    push!(K[add],arrival[Ea])
        
    N = Dict{Tuple{Int64,Int64},Array{Int64}}()
    
    for j in 1:num_aircraft
        for k in collect(keys(K))
            N[j,k] = [1];
            for r in 1:num_route
                if (departure[r] >= K[k][1] && departure[r] < K[k][2]) || (arrival[r] > K[k][1] && departure[r] < K[k][2])
                    if ~in(r,N[j,k])
                        push!(N[j,k],r)
                    end
                    if k >= 2 && in(1,N[j,k])
                        shift!(N[j,k])
                    end
                end
            end
        end
    end
                
    M = Dict{Tuple{Int64,Int64},Array{Int64}}()
    
    for (k,v) in N
        for i in v
            try
                M[k[1],i]
            catch error
                if isa(error, KeyError)
                    M[k[1],i] = [k[2]];
                end
            end
            if ~isempty(M[k[1],i]) && ~in(k[2],M[k[1],i])
                push!(M[k[1],i],k[2])
            end
        end
    end

    
 return departure, arrival, profit, K, N, M
        
end

function comprobar_fact(X,N)
    
    ########################################
    #
    # Esta rutina comprueba si la solucion obtenida de la relajacion es factible o no.
    # Si la solucion Lagrangiana es infactible es porque algun avion fue asignado a mutiples rutas
    # en el mismo periodo k
    #
    ########################################
    
    # input
    # X: solucion
    # N
    
    # Output
    # X
    # conflict
    # sum_conflict
    
    conflict = Dict{Tuple{Int64,Int64},Array{Int64}}((j,k) => [0] for j in 1:num_aircraft, k in 1:num_periods)
    
    sum_conflict = []
    for j in 1:num_aircraft
        for k in 1:num_periods
            num_conflict = sum(X[i,j] for i in N[j,k])
            push!(sum_conflict,if  num_conflict > 0  num_conflict - 1 else 0 end)
            #println("aircraft $(j) period $(k) conflic $(if  num_conflict > 0  num_conflict - 1 else 0 end)")
            for i in N[j,k]
                if X[i,j] == 1
                    push!(conflict[j,k],i)
                    if conflict[j,k][1] == 0
                        shift!(conflict[j,k])
                    end
                end
            end
        end
    end
     
    return conflict,sum(sum_conflict)
end

function heuristic_fac(conflict,X,N,num_aircraft)

    #########################################
    #
    # Esta heuristica busca una solcion factible (SB)
    # Usa dos forma: reasignar y desasignar
    #
    #########################################
        
    for (k,v) in conflict
        v_aux = copy(v)
        if length(v) >= 2
            for i in 1:length(v)-1
                #println(i)
                adis = [sum(X[ii,dis] for ii in N[dis,k[2]]) for dis in 1:num_aircraft] # busco avion disponible
                padis = indmin(adis)
                #println("Posible avion disponible $(padis)")
                if adis[padis] == 0 # al menos algun avion j esta disponible, posible cambio a la ruta i
                
                    avion = find( x->(x == 0), adis) #
                    #println(avion) #
                
                    while (~isempty(avion)) #
                 
                    #padis = avion[1] #
                    padis = sample(avion,1,replace=false)[1]
                
                    ss = indmin(p[v_aux,k[1]] - p[v_aux,padis]) # busco la minima diferencia entre los v
                    #println("Se puede cambiar la ruta $(v_aux[ss]) al avion $(padis) en el periodo $(k[2])?")
                    best_act_conflict = []
                    for all in 1:num_periods
                        push!(best_act_conflict,sum(X[id,padis] for id in N[padis,all]))
                        if best_act_conflict[end] >= 1
                            best_act_conflict[end] = best_act_conflict[end] - 1
                        end
                    end
                    sbest_act_conflict = sum(best_act_conflict)
                    #println("conflictos $(sbest_act_conflict)")
                    X[v_aux[ss],k[1]] = 0    # quito un conflicto del avion actual
                    X[v_aux[ss],padis] = 1   # valor propuesto
                    new_act_conflict = []
                    for all in 1:num_periods
                        push!(new_act_conflict,sum(X[id,padis] for id in N[padis,all]))
                        if new_act_conflict[end] >= 1
                            new_act_conflict[end] = new_act_conflict[end] - 1
                        end
                    end
                    snew_act_conflict = sum(new_act_conflict)
                    #println("nuevo conflictos $(snew_act_conflict)")
                    if  snew_act_conflict <= sbest_act_conflict
                        #println("mejora")
                        break #
                    else
                        #println("no mejora")
                        X[v_aux[ss],padis] = 0 # quito el valor propuesto
                    end
                
                    shift!(avion) #
                
                    end
                    
                    deleteat!(conflict[k], conflict[k] .== [v_aux[ss]])
                    deleteat!(v_aux, v_aux .== [v_aux[ss]])
                else
                    #println("El avion $(padis) no esta disponible")
                    ss = indmin(p[v_aux,k[1]]) # busco el de menor beneficio
                    X[v_aux[ss],k[1]] = 0 # desasigno el valor de la ruta i al avion j
                    deleteat!(conflict[k], conflict[k] .== [v_aux[ss]]) # borro de conflicto
                    deleteat!(v_aux, v_aux .== [v_aux[ss]]) # borro elemento del vector auxiliar
                end
            end
        end
    end
    
    conflict,total_conflict = comprobar_fact(X,N);
    
    # una ruta no puede ser asignada mas de dos veces
    
    route_satisfechas = []
    aircraft_used = []

    for j in 1:num_aircraft
        for i in 1:num_route
            if X[i,j] == 1
                if ~in(j,aircraft_used)
                    push!(aircraft_used,j)
                end
                if ~in(i,route_satisfechas)
                    push!(route_satisfechas,i);
                else
                    X[i,j] = 0
                end
            end
        end
    end
    
    if total_conflict == 0
        #println("Se construyo una solucion factible")
        SB = sum(p[i,j]*X[i,j] for i in 1:num_route,j in 1:num_aircraft)
    else
        #println("No se pudo construir una solucion factible")
    end
    
    return SB,X,conflict,total_conflict
    
end
    

# se define el solver, numero de aviones, numero de rutas y la duracion promedio de las rutas

env = Gurobi.Env()

sub_solver   = GurobiSolver(env,OutputFlag=0)
num_aircraft = 40
num_route    = 800
avg_route    = 8

d,a,p,K,N,M  = genData(avg_route,num_route,num_aircraft);

# se define el numero de periodos
num_periods  = maximum(collect(keys(K)));




# initialize:

πo         = 2
iter       = 0
iter_L     = 0
max_iter   = 500
max_iter_L = 40
contar_L   = 0
πo_min     = 5e-7
SB         = 0
best       = 0
XB         = Dict{Tuple{Int64,Int64},Int64}((i,j) => 0 for i in 1:num_route, j in 1:num_aircraft)
t0         = time();

β = Dict{Tuple{Int64,Int64},Float64}()

for j in 1:num_aircraft
    for k in 1:num_periods
        β[j,k] = maximum(p)
    end
end

UB = sum(β[j,k] for j in 1:num_aircraft,k in 1:num_periods)
LB = 0

println("iter \tSB \tUB \tGAP")
println("$(iter) \t$(SB) \t$(UB) \t$((UB-SB)/UB*100)")

G0 = Dict{Tuple{Int64,Int64},Float64}((j,k) => 0 for j in 1:num_aircraft,k in 1:num_periods)

for con in 1:max_iter
    
    output = pmap(loadAndSolveSubProblem, 
        [p[i,:] for i = 1:num_route], 
        [M for i = 1:num_route],
        [β for i = 1:num_route],
        [collect(keys(K)) for i=1:num_route],
        [i for i=1:num_route],
        [sub_solver for i = 1:num_route]
        );
    
    # get solution
    X = Dict{Tuple{Int64,Int64},Int64}((i,j) => output[i][4][j] for i in 1:num_route, j in 1:num_aircraft)
    
    UB = min(UB,sum((p[i,j] - sum(β[j,k] for k in 1:num_periods if k in M[j,i]))*X[i,j] for i=1:num_route,j in 1:num_aircraft) + sum(β[j,k] for k in 1:num_periods, j in 1:num_aircraft))

    G = Dict{Tuple{Int64,Int64},Float64}()
                    
    for j in 1:num_aircraft
        for k in 1:num_periods
            G[j,k] = 1 - sum(X[i,j] for i in N[j,k]) + 0.25*G0[j,k]
        end
    end
    
    #T = πo * (UB - SB)/(sum(G[j,k] for j in 1:num_aircraft, k in 1:num_periods))^2
    T = πo * (UB - SB)/(sum(abs(G[j,k]) for j in 1:num_aircraft, k in 1:num_periods))
                                
    for j in 1:num_aircraft
        for k in 1:num_periods
            β[j,k] = max(0,β[j,k] - T*G[j,k])
        end
    end
    
    conflict,total_conflict = comprobar_fact(X,N);    
    
    if total_conflict == 0
        #println("La solucion de la relajacion es factible")
        LB = max(LB,sum(p[i,j]*X[i,j] for i in 1:num_route,j in 1:num_aircraft));
        if LB >= max(SB,LB)
            SB = LB
            XB = X
            best = SB
        end
    else
        #println("La solucion de la relajacion es infactible")
        SB_NEW,X_NEW,conflict,total_conflict = heuristic_fac(conflict,X,N,num_aircraft)
        
        # verifico si hay aviones y rutas que estan disponibles para ser asignadas
        route_satisfechas = []
        aircraft_used = []

        for j in 1:num_aircraft
            for i in 1:num_route
                if X[i,j] == 1
                    if ~in(j,aircraft_used)
                        push!(aircraft_used,j)
                    end
                    push!(route_satisfechas,i);
                end
            end
        end
                                
        act_route = setdiff(1:num_route,route_satisfechas)
        act_aircraft = setdiff(1:num_aircraft,aircraft_used)

        if ~isempty(act_route) && ~isempty(act_aircraft)
            for j in act_aircraft
                route_satisfechas = []
                time_end = Dates.DateTime(d[act_route[1]]-Dates.Hour(5))
                for i in act_route
                    if d[i] >= time_end
                        X_NEW[i,j] = 1
                        push!(route_satisfechas,i)
                        time_end = a[i]
                    end
                end
                act_route = setdiff(copy(act_route),route_satisfechas)
            end
        end

        conflict,total_conflict = comprobar_fact(X_NEW,N);
    
        SB_NEW,X_NEW,conflict,total_conflict = heuristic_fac(conflict,X_NEW,N,num_aircraft)
                                    
        if total_conflict >= 1
            println("Infactible solucion")
            break
        end
        if SB_NEW >= max(SB,LB)
            SB = SB_NEW
            XB = X_NEW
            best = SB
        end
    end

    if best == SB
        contar_L = contar_L + 1;
        if contar_L == max_iter_L
            contar_L = 0
            πo = πo/2
        end
        #println(contar_L)
    else
        contar_L = 0;
    end                                

    println("$(con) \t$(SB) \t$(UB) \t$((UB-SB)/UB*100)")

    X = XB
    G0 = Dict{Tuple{Int64,Int64},Float64}((j,k) => G[j,k] for j in 1:num_aircraft,k in 1:num_periods)

    if πo_min > πo || (UB-SB)/UB*100 < 1e-5
        break
    end
    
end
                            
timefinal = time() - t0
println("Elapsed Time: ", timefinal, " seconds")




@time output = pmap(loadAndSolveProblem, 
        [p for i = 1:1], 
        [K for i = 1:1],
        [N for i = 1:1],
        [sub_solver for i = 1:1]
        )



route_satisfechas = []
aircraft_used = []

for j in 1:num_aircraft
    println("Aircraft $(j)")
    for i in 1:num_route
        if XB[i,j] == 1
            if ~in(j,aircraft_used)
                push!(aircraft_used,j)
            end
            if ~in(i,route_satisfechas)
                push!(route_satisfechas,i);
            else
                println("Error la misma ruta asignada dos veces")
                break
            end
            println("route $(i) - \t$(d[i]) - $(a[i])")
        end
    end
end

println("satisfied route: $(length(route_satisfechas)/num_route*100)%")
println("Aircraft used: $(length(aircraft_used)/num_aircraft*100)%")


     
