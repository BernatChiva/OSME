# f.javier.heredia(at)upc.edu

 
# Parameters
param nT;
set T = 1..nT;             # Set of time periods
param nB;
set B = 1..nB;
set G;                     # Set of generation units
set D;                     # Set of demand (loads)
set L within B cross B;    # Transmission lines
set refB within B;         # Reference bus.
param sb;                # Base power (MVAr)
param x     { L };       # Line reactance (p.u.)
param smax  { L };       # Line capacity (MW)

set GB{B} default {};      # Generation units at bus "k"                                       
set DB{B} default {};      # Loads units at bus "k"

param nbG;
param nbD;
set bG = 1..nbG;           # Set of generation bid blocks
set bD = 1..nbD;           # Set of demand bid blocks

#set S within G default {}; # Set of slack generators. 
param pgmin { G };         # Minimum power output (MW)
param pgmax { G };         # Minimum power output (MW)

param rd  { G };           # Downward ramp (MW/h)
param ru  { G };           # Upward ramp (MW/h)
param tU    { G };         # minimum on time.
param tD    { G };         # minimum off time.
param pg0 { i in G };      # Initial power output (MW)
param u0    { G };         # Initial state

param pbG { G, bG, T};
param lbG { G, bG, T};

param pbD { D, bD, T};
param lbD { D, bD, T};

# Decision variables
var P    { (k,l) in L, T } >= -smax[k,l], <= smax[k,l];  # Load flow line "(k,l)" (MW)
var Th {B, T};
var PG  { G, bG, T } >=0;                               # Matched generation bid unit "i" block "b" (MW)
var PGT   { G, {0} union T } >=0;                         # Total matched generation bid of unit "i"  (MW)
var U    { G, {0} union T} binary;                       # Unit state 
var PD  { D, bD, T } >=0;                               # Matched demand bid unit "i" block "b" (MW)
var PDT   { D, T } >=0;                                   # Total matched demand bid of unit "i"  (MW)



# Objective function
var SW{t in T}  =
	sum{i in D, b in bD} lbD[i,b,t]*PD[i,b,t]- sum{i in G, b in bG} lbG[i,b,t]*PG[i,b,t];

# Objective function
maximize Total_SW :   sum{t in T} SW[t];
 
# Matched generation and demand:
s.t. Matched_G { i in G, b in bG, t in T } : PG[i,b,t] <= pbG[i,b,t];
s.t. Matched_D { i in D, b in bD, t in T } : PD[i,b,t] <= pbD[i,b,t];
# Total matched generation
s.t. Total_Matched_G {i in G, t in T}: PGT[i,t] = sum{b in bG}PG[i,b,t];
# Total matched demand
s.t. Total_Matched_D {i in D, t in T}: PDT[i,t] = sum{b in bD}PD[i,b,t];
# Nodal market equilibrium
s.t. NPFE { k in B, t in T} :
	sum{ (k,l) in L } P[k,l,t] - sum{ (l,k) in L} P[l,k,t] =
	sum{ i in GB[k]} PGT[i,t] - sum{ i in DB[k]} PDT[i,t];
# Line power flow equations
s.t. LPFE{(k, l) in L, t in T}: P[k,l,t] = sb*(Th[k, t] - Th[l, t])/x[k,l];
# Ramp limits
s.t. Ramp { i in G, t in T } : -rd[i] <= PGT[i,t] - PGT[i, t-1] <= ru[i];
# Min power output
s.t. PGMin { i in G, t in T } : pgmin[i]*U[i,t] <= PGT[i,t];
# Max power output
s.t. PGMax { i in G, t in T } : PGT[i,t] <= pgmax[i]*U[i,t];
# Reference bus
s.t. RefB { k  in refB, t in T } : Th[k,t] = 0.0;
# Initial power output
s.t. Initial_PO {i in G}: PGT[i,0] = pg0[i];

######### Unit commitment ##############
var V    { G, T } binary;                                # Start-up variable.
var W    { G, T } binary;								# shut-down variable.
# Unit commitment: as in Morales et al. OR Spectrum (2015) DOI 10.1007/s00291-015-0400-4
#s.t. UCdef {i in G, t in T} : U[i,t] - U[i,t-1] = V[i,t] - W[i,t];
#s.t. UCtU {i in G, t in tU[i]..nT} : sum{l in (t-tU[i]+1)..t} V[i,l] <= U[i,t];
#s.t. UCtD {i in G, t in tD[i]..nT} : sum{l in (t-tD[i]+1)..t} W[i,l] <= 1-U[i,t];
# Initial state:
#s.t. U_0  { i in G }: U[i,0]  = u0[i];

########## MIC ########################
param price{T} default 40;
param nit >= 0;							 # number of iterations of the heuristic
param perc_mic default 0.9;  # percentage of the initial offer to use or not the mic constraints
param aux_price{T} default 0;
param aux_price2{T} default 0;  # to avoid bouncing back
var U_total{G} binary; # 1 si la maquina sha engegat algun instant
s.t. Turned_on {i in G, t in T}: U_total[i] >= U[i,t];			
s.t. MIC{i in G}: 1 - U_total[i] >=1 -(sum{b in bG, t in T}price[t]*PG[i,b,t])/
									(perc_mic*sum{b in bG, t in T}lbG[i,b,t]*pbG[i,b,t]);	
									# aqui es pot sumar epsilon al denominador per evitar 
									# que peti quan perc_mic = 0							