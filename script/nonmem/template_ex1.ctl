$PROB RUN# run_num

$INPUT ID TIME EVID MDV CMT AMT SS II DV WT

$DATA data_fname IGNORE=@

$SUB ADVAN2 TRANS2

$PK
  CL_WT     = LOG(WT / 70) * THETA(4)
  V_WT      = LOG(WT / 70) * THETA(5)

  CL       = THETA(1) * EXP(CL_WT + ETA(1))
  V        = THETA(2) * EXP(V_WT  + ETA(2))
  KA       = THETA(3) * EXP(ETA(3))
  F1       = 1

  S2       = V

$ERROR
  IPRED = F
  Y     = IPRED*(1 + EPS(1)) + EPS(2)

$THETA
  (0,10)    ; CL
  (0,100)   ; V
  (0,0.25)  ; KA

  0.75 FIX  ; CL_WT
  1 FIX     ; V_WT

$OMEGA
  0.08      ; CL
  0.1       ; V
  0.2       ; KA

$SIGMA
  0.05      ; [P] Proportional
  1         ; [A] Additive

$SIMULATION (run_num)
$ESTIMATION METHOD=1 INTER PRINT=1 MSFO=./run_num.msf

