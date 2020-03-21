[ pkmodel ] cmt="DEPOT CENT PERI", depot=TRUE, trans = 11

[ param ]
CL=20, VC=70, Q=1, VP=30, KA = 0.25, WT= 70

[ main ]
double CLi = CL*pow(WT/70,0.75);
double V2i = VC*(WT/70);
double Qi  = Q*pow(WT/70,0.75);
double V3i = VP*(WT/70);
double KAi = KA;

[ table ]
capture CP = CENT/V2i ;