#title Portfolio Optimization for Electric Utilities (POUTIL,SEQ=342)

#We discuss a portfolio optimization problem occurring in the energy
#market. Energy distributing public services have to decide how much
#of the requested energy demand has to be produced in their own power
#plant, and which complementary amount has to be bought from the spot
#market and from load following contracts.

#This problem is formulated as a mixed-integer linear programming
#problem and implemented in GAMS. The formulation is applied to real data
#of a German electricity distributor.

#Most equations contain the reference number of the formula in the
#publication.

#Rebennack, S, Kallrath, J, and Pardalos, P M, Energy Portfolio
#Optimization for Electric Utilities: Case Study for Germany. In
#Bj�rndal, E, Bj�rndal, M, Pardalos, P.M. and R�nnqvist, M Eds,.
#Springer, pp. 221-246, 2010.

#Keywords: mixed integer linear programming, energy economics, portfolio optimization,
#          unit commitment, economic dispatch, power plant control, day-ahead market

using JuMP, Gurobi;
using GLPKMathProgInterface;

#solver = GurobiSolver()
# OutputFlag=0

solver = GLPKSolverMIP(presolve=true)

# Set
t = [string("t$(i)") for i=1:96];

# Parameters

# electric power forecast
PowerForecast = Dict{String,Integer}(
"t1" => 287,"t2" => 275,"t3" => 262,"t4" => 250,"t5" => 255,"t6" => 260,"t7" => 265,"t8" => 270,"t9" => 267,"t10"=> 265,
"t11"=> 262,"t12"=> 260,"t13"=> 262,"t14"=> 265,"t15"=> 267,"t16"=> 270,"t17"=> 277,"t18"=> 285,"t19"=> 292,"t20"=> 300,
"t21"=> 310,"t22"=> 320,"t23"=> 330,"t24"=> 340,"t25"=> 357,"t26"=> 375,"t27"=> 392,"t28"=> 410,"t29"=> 405,"t30"=> 400,
"t31"=> 395,"t32"=> 390,"t33"=> 400,"t34"=> 410,"t35"=> 420,"t36"=> 430,"t37"=> 428,"t38"=> 427,"t39"=> 426,"t40"=> 425,
"t41"=> 432,"t42"=> 440,"t43"=> 447,"t44"=> 455,"t45"=> 458,"t46"=> 462,"t47"=> 466,"t48"=> 470,"t49"=> 466,"t50"=> 462,
"t51"=> 458,"t52"=> 455,"t53"=> 446,"t54"=> 437,"t55"=> 428,"t56"=> 420,"t57"=> 416,"t58"=> 412,"t59"=> 408,"t60"=> 405,
"t61"=> 396,"t62"=> 387,"t63"=> 378,"t64"=> 370,"t65"=> 375,"t66"=> 380,"t67"=> 385,"t68"=> 390,"t69"=> 383,"t70"=> 377,
"t71"=> 371,"t72"=> 365,"t73"=> 368,"t74"=> 372,"t75"=> 376,"t76"=> 380,"t77"=> 386,"t78"=> 392,"t79"=> 398,"t80"=> 405,
"t81"=> 408,"t82"=> 412,"t83"=> 416,"t84"=> 420,"t85"=> 413,"t86"=> 407,"t87"=> 401,"t88"=> 395,"t89"=> 386,"t90"=> 377, 
"t91"=> 368,"t92"=> 360,"t93"=> 345,"t94"=> 330,"t95"=> 315,"t96"=> 300 
);


# Power Plant (PP)
cPPvar = 25.0;    # variable cost of power plant [euro / MWh]
pPPMax = 300.0;   # maximal capacity of power plant      [MW]

# Set
 m = [string("m$(i)") for i=1:8];    # stage of the power plant
iS = [string("iS$(i)") for i=0:8];   # interval for constant PP operation
iI = [string("iI$(i)") for i=0:16];  # length of idle time period

# title Spot Market (SM)
cBL = 32.0; # cost for one base load contract [euro / MWh]
cPL = 41.0; # cost for one peak load contract [euro / MWh]

#Parameter IPL(t) indicator function for peak load contracts
IPL = Dict{String,Integer}();

for i=1:length(t)
    IPL[t[i]] = if (i>= 33 && i<= 80) 1 else 0 end;
end 

# title Load following Contract (LFC)
pLFCref = 400.0 # power reference level for the LFC

# Set
b = [string("b$(i)") for i=1:3];  # support points of the zone prices

# amount of energy at support point b
eLFCbY = Dict{String,Integer}(
"b1"=> 54750,
"b2"=> 182500,
"b3"=> 9000000
);

# specific energy price in segment b
cLFCvar = Dict{String,Integer}(
"b1"=> 80.0,
"b2"=> 65.0,
"b3"=> 52.0
);

# daily border of energy volumes for LFC
eLFCb = Dict{String,Float64}();

# calculate the daily borders of the energy volumes for the zones
for k in eLFCbY
    eLFCb[k[1]] = eLFCbY[k[1]]/365;
end

# accumulated cost for LFC up to segment b
cLFCs = Dict{String,Float64}();

# calculate the accumulated cost
cLFCs["b1"]         = 0;
cLFCs["b2"]         = cLFCvar["b1"]*eLFCb["b1"];

for i=3:length(b)
    cLFCs[b[i]] = cLFCs[b[i-1]] + cLFCvar[b[i-1]]*(eLFCb[b[i-1]] - eLFCb[b[i-2]]);
end


# variables

poutil = Model(solver=solver)

# total cost
@variable(poutil, c >= 0);

# cost of PP usage
@variable(poutil, cPP >= 0);

# power withdrawn from power plant
@variable(poutil, pPP[t] >= 0);

# indicate if the PP is in stage m at time t
@variable(poutil, delta[m,t] >= 0, Bin);

# indicator for segment b (for zone prices)
@variable(poutil, mu[b] >= 0, Bin);

# indicate if there is a PP stage change
@variable(poutil, chiS[t] >= 0);

# indicate if the PP left the idle stage
@variable(poutil, chiI[t] >= 0);

# cost of energy from SM
@variable(poutil, cSM >= 0);

# power from the spot market
@variable(poutil, pSM[t] >= 0);

# quantity of base load contracts
@variable(poutil, 0 <= alpha <= maximum([PowerForecast[i] for i=t]), Int);

# quantity of peak load contracts
@variable(poutil, 0 <= beta <= maximum([PowerForecast[i] for i=t]), Int);

# cost of LFC which is the enery rate
@variable(poutil, cLFC >= 0);

# total energy amount of LFC
@variable(poutil, eLFCtot >= 0);

# energy from LFC in segment b
@variable(poutil, eLFCs[b] >= 0);

# power from the LFC
@variable(poutil, 0 <= pLFC[t] <= pLFCref);

# CONSTRAINT

# the objective function: total cost; eq. (6)
@constraint(poutil, obj, c == cPP + cSM + cLFC);

# meet the power demand for each time period exactly; eq. (23)
@constraint(poutil, demand[i=t], pPP[i] + pSM[i] + pLFC[i] == PowerForecast[i]);

# (fix cost +) variable cost * energy amount produced; eq. (7) & (8)
@constraint(poutil, PPcost, cPP == cPPvar*sum(0.25*pPP[i] for i=t));

# power produced by the power plant; eq. (26)
@constraint(poutil, PPpower[i=t], pPP[i] == pPPMax * sum(0.1*(s + 2)*delta[m[s],i] for s=2:length(m)));

# the power plant is in exactly one stage at any time; eq. (25)
@constraint(poutil, PPstage[i=t], sum(delta[s,i] for s=m) == 1);

# next constraints model the minimum time period a power plant is in the
# same state and the constraint of the minimum idle time
# we need variable 'chiS' to find out when a status change takes place
# eq. (27)
@constraint(poutil, PPchiS1[i=2:length(t),j=1:length(m)],chiS[t[i]] >= delta[m[j],t[i]] - delta[m[j],t[i-1]]);

# second constraint for 'chiS' variable; eq. (28)
@constraint(poutil, PPchiS2[i=2:length(t),j=1:length(m)],chiS[t[i]] >= delta[m[j],t[i-1]] - delta[m[j],t[i]]);

# control the minimum change time period; eq. (29)
@constraint(poutil, PPstageChange[i=1:(length(t)-length(iS)+2)], sum(chiS[t[i + s]]  for s=0:(length(iS)-2)) <= 1);

# indicate if the plant left the idle state; eq. (30)
@constraint(poutil, PPstarted[i=2:length(t)],chiI[t[i]] >= delta["m1",t[i-1]] - delta["m1",t[i]]);

# control the minimum idle time period:
# it has to be at least Nk2 time periods long; eq. (31)
@constraint(poutil, PPidleTime[i=1:(length(t)-length(iI)+2)], sum(chiI[t[i + s]]  for s=0:(length(iI)-2)) <= 1);

# cost for the spot market; eq. (12)
# consistent of the base load (alpha) and peak load (beta) contracts
@constraint(poutil, SMcost, cSM == 24*cBL*alpha + 12*cPL*beta);

# Spot Market power contribution; eq. (9)
@constraint(poutil, SMpower[i=t], pSM[i] == alpha + IPL[i]*beta);

# cost of the LFC is given by the energy rate; eq. (14) & (21)
@constraint(poutil, LFCcost, cLFC == sum(cLFCs[i]*mu[i] + cLFCvar[i]*eLFCs[i] for i=b));

# total energy from the LFC; eq. (16)
# connect the eLFC(t) variables with eLFCtot
@constraint(poutil, LFCenergy, eLFCtot == 0.25*sum(pLFC));

# indicator variable 'mu':
# we are in exactly one price segment b; eq. (18)
@constraint(poutil, LFCmu, sum(mu) == 1);

# connect the 'mu' variables with the total energy amount; eq. (19)
@constraint(poutil, LFCenergyS, eLFCtot == sum(eLFCb[b[i-1]]*mu[b[i]]  for i=2:length(b)) + sum(eLFCs));

# accumulated energy amount for segment "b1"; eq. (20)
@constraint(poutil, LFCemuo, eLFCs["b1"] <= eLFCb["b1"]*mu["b1"]);

# accumulated energy amount for all other segments (then "b1"); eq. (20)
@constraint(poutil, LFCemug[i=2:length(b)],eLFCs[b[i]] <= (eLFCb[b[i]] - eLFCb[b[i-1]])*mu[b[i]]);

@objective(poutil, Min, c);

