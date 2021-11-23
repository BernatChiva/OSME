#######################################################
# Generalized Unit Commitment. Model file.
# F.-Javier Heredia, http://gnom.upc.edu/heredia
# Code under Creative Commons License CC BY-NC-ND 3.0
# http://creativecommons.org/licenses/by-nc-nd/3.0/
#######################################################

# Parameters
param nB;
param nT;
set G;                     # Set if generation units
set D;                     # Set of demand (loads)
set B = 1..nB;             # Set of busses
set T = 1..nT;             # Set of time periods

set refB within B;         # Reference bus.
set GB{B} default {};      # Generation units at bus "k"                                       
set DB{B} default {};      # Loads units at bus "k"
set L within B cross B;    # Transmission lines
 
param cq    { G };         # Quadratic cost coeficients (BTu/hMW2)
param cl    { G };         # Linear cost coeficients (BTu/MWh)
param cb    { G };         # Constant cost coeficients (BTu/h)
param csud  { G };         # Constant start-up/down cost coeficients (BTu)
param pgmin { G };         # Minimum power output (MW)
param pgmax { G };         # Minimum power output (MW)
param rd    { G };         # Downward ramp (MW/h)
param ru    { G };         # Upward ramp (MW/h)
param pg0   { i in G };    # Initial power output (MW)
param u0    { G };         # Initial state
param tU    { G };         # minimum on time.
param tD    { G };         # minimum off time.

param pd    { D, T };      # Total system load (MW)
param pdB   { k in B, t in T} := sum{i in DB [k]} pd[i,t];

param maxLS >=0, <=1 default 0.2;  # Maximum load curtailment in Load shedding.
param pLS >= 0;                    # Price of the load curtailment [€/MWh].

param x     { L };       # Line reactance (p.u.)
param r     { L };       # Line reactance (p.u.)
param smax  { L };       # Line capacity (MW)
param sb;                # Base power (MVAr)

# Decision variables
var PG   { G, {0} union T };                             # Power generation of unit "i" (MW)
var P    { (k,l) in L, T } >= -smax[k,l], <= smax[k,l];  # Load flow line "(k,l)" (MW)
var Th   { B, T };                                       # Voltage angle (p.u.)
var U    { G, {0} union T} binary;                       # Unit state 
var V    { G, T } binary;                                # Start-up variable.
var W    { G, T } binary;                                # shut-down variable.
var LS { i in D, t in T } >=0, <= maxLS;                 # Load shedding: load curtailment 
                                                         # as a fraction of the total load  

# Objective function
minimize Total_Cost : 
	sum{ i in G, t in T } (cq[i]*PG[i,t]^2+cl[i]*PG[i,t] + cb[i]*U[i,t] + csud[i]*(V[i,t]+W[i,t]))
	+pLS*sum{ i in D, t in T} LS[i,t]*pd[i,t];

# Nodal Power Flow Equations with load shedding:
s.t. NPFE { k in B, t in T} :
	sum{ (k,l) in L } P[k,l,t] - sum{ (l,k) in L} P[l,k,t] =
	sum{ i in GB[k]} PG[i,t] - sum{ i in DB[k]} (1-LS[i,t])*pd[i,t];

# Line Power Flow Equations
s.t. LPFE { (k,l)  in L, t in T } : P[k,l,t] = sb*( Th[k,t]- Th[l,t] ) / x[k,l];

# Reference bus
s.t. RefB { k  in refB, t in T } : Th[k,t] = 0.0;

# Min power output
s.t. PGMin { i in G, t in T } : pgmin[i]*U[i,t] <= PG[i,t];

# Max power output
s.t. PGMax { i in G, t in T } : PG[i,t] <= pgmax[i]*U[i,t];

# Ramp limits
s.t. Ramp { i in G, t in T} : -rd[i] <= PG[i,t] - PG[i,t-1] <= ru[i];

# Unit commitment: as in Morales et al. OR Spectrum (2015) DOI 10.1007/s00291-015-0400-4
s.t. UCdef {i in G, t in T} : U[i,t] - U[i,t-1] = V[i,t] - W[i,t];
s.t. UCtU {i in G, t in tU[i]..nT} : sum{l in (t-tU[i]+1)..t} V[i,l] <= U[i,t];
s.t. UCtD {i in G, t in tD[i]..nT} : sum{l in (t-tD[i]+1)..t} W[i,l] <= 1-U[i,t];

# Initial state:
s.t. U_0  { i in G }: U[i,0]  = u0[i];
s.t. PG_0 { i in G }: PG[i,0] = pg0[i];