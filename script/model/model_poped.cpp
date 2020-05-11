[ param ]
CL = 1, VMAX = 10, KM = 10, V1 = 8, Q = 10, V2 = 100

[ cmt ] CENT PERIPH

[ main ]
double ke  = CL/V1;
double k12 = Q/V1;
double k21 = Q/V2;

[ ode ]
double CP = CENT/V1;

dxdt_CENT = k21*PERIPH - k12*CENT - VMAX*CP/(KM + CP) - ke*CENT;
dxdt_PERIPH = k12*CENT - k21*PERIPH;

[ capture ]
CP
