[PROB]

1. Minimal PBPK Metamodeling of Target Engagement in Skin 

[PARAM] @annotated

ka : 0.18 : absorption rate constant
CL : 0.154 : clearance rate
F : 0.729 : biovailability


// Lymph flow rates
// units in liters per day
L : 2.9 : Total lymph flow rate
Ls: 0.247 : skin lymph flow rate
L1: 0.71 : muscle lymph flow rate
L2: 1.943 : Leaky tissue lymph flow rate

// Volume of compartments
// units in liters

Vp : 2.6 : plasma volume
V2 : 4.37 : leaky tissue volume
Vs : 1.81 : skin volume
V1 : 6.3 : muscle volume
VL : 2.6 : lymph volume

// reflection coefficients

sk : 0.630 : skin
s1 : 0.95 : muscle
s2 : 0.363 : leaky
sL : 0.2 : lymph

// PD (enzyme parameters)

kdeg_skin : 2.44 : degrad skin
kdeg_plasma : 45.5 : degrad plasma

kint_skin : 0.34 : inte skin
kint_plasma : 1.24 : inte plasma

E0_skin : 0.28 : baseline skin
E0_plasma : 0.015 : baseline plasma

ksyn : 0.683 : syn of IL-17a
Kd : 129 : dissociation constant

n : 5: hill coeff
Emax : 1 : Max response
EC50 : 1 : Half response
[CMT] @annotated

AUC : Area under curve?

SERUM: SERUM compartment
SKIN : SKIN compartment
MUSC : MUSCLE compartment
LEAKY: Leaky compartment
ABS : Absorption compartment
LYMPH: LYMPH compartment
T1 : Total IL7A serum compartment
T2 : Total IL17A skin compartment
T3 : Drug-IL7A serum compartment
T4 : DRug-IL7a skin compartment


[ODE]

// Four algebraic equations


double Cp_f = 0.5*((SERUM - T1 - Kd) + (pow(pow(SERUM - T1 - Kd,2) + 4*SERUM*Kd,0.5)));
double Csk_f = 0.5*((SKIN - T2 - Kd) + (pow(pow(SKIN - T2 - Kd,2) + 4*SKIN*Kd,0.5)));

double AR_sk = T2*Csk_f/(Kd + Csk_f);
double AR_p = T1*Cp_f/(Kd + Cp_f);

// Response for skin and plasma

double R_skin = E0_skin + pow(SKIN,n)*(Emax - E0_skin)/(pow(SKIN,n) + pow(EC50,n));
double R_plasma = E0_plasma + pow(SERUM,n)*(Emax - E0_plasma)/(pow(SERUM,n) + pow(EC50,n));

dxdt_ABS = -ka*ABS;

dxdt_SERUM = (ka*F*ABS - Cp_f*Ls*(1-sk) - Cp_f*L1*(1-s1) - Cp_f*L2*(1-s2) - Cp_f*CL - kint_plasma*AR_p*Vp)/Vp;

dxdt_SKIN = (Cp_f*Ls*(1-sk) - Csk_f*Ls*(1-sL) - kint_skin*AR_sk*Vs)/Vs;

dxdt_MUSC = (Cp_f*L1*(1-s1) - MUSC*L1*(1-sL))/V1;

dxdt_LEAKY = (Cp_f*L2*(1-s2) - LEAKY*L2*(1-sL))/V2;

dxdt_LYMPH = (Csk_f*Ls*(1-sL) + MUSC*L1*(1-sL) + LEAKY*L2*(1-sL) - LYMPH*L)/VL;

dxdt_T1 = ksyn - kdeg_plasma*(T1-AR_p) - kint_plasma*AR_p;
dxdt_T2 = ksyn - kdeg_skin*(T2 - AR_sk) - kint_skin*AR_sk;

dxdt_AUC = SERUM/Vp;

if(SS_ADVANCE) dxdt_AUC = 0;

[TABLE] capture Cserum = SERUM;
capture resp = R_skin;
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
