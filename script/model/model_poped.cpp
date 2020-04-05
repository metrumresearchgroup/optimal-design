[ pkmodel ] cmt="DEPOT CENT PERI", depot=TRUE, trans = 11

[ param ]
CL=10, V2=100, Q=1, V3=30, KA = 0.25, WT= 70,
  wt_cl = 0.75, wt_v = 1

[ main ]
double CLi = CL * pow(WT/70, wt_cl);
double V2i = V2 * pow(WT/70, wt_v);
double Qi  = Q * pow(WT/70, wt_cl);
double V3i = V3 * pow(WT/70, wt_v);
double KAi = KA;

[ table ]
capture CP = CENT/V2i ;