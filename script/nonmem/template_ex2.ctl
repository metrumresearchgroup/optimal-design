$PROB RUN# run_num

$INPUT ID TIME EVID MDV CMT AMT DV

$DATA data_fname IGNORE=@

$SUB ADVAN13 TOL=10
$MODEL NCOMPARTMENTS = 2

$PK
  CL      = THETA(1) * EXP(ETA(1))
  VMAX    = THETA(2) * EXP(ETA(2))
  KM      = THETA(3) * EXP(ETA(3))
  V1      = THETA(4) * EXP(ETA(4))
  Q       = THETA(5) * EXP(ETA(5))
  V2      = THETA(6) * EXP(ETA(6))

  KE      = CL/V1
  K12     = Q/V1
  K21     = Q/V2

$DES
  CP      = A(1)/V1
  DADT(1) = K21*A(2) - K12*A(1) - VMAX*CP/(KM + CP) - KE*A(1)
  DADT(2) = K12*A(1) - K21*A(2)

$ERROR
  IPRED   = A(1)/V1
  Y       = IPRED*(1 + EPS(1))
  IF (ICALL.EQ.4) THEN
    IF (Y < 0.001) MDV = 1
  ENDIF

$THETA

  (0,0.5,5)  ; CL
  (0,20,100) ; VMAX
  (0,1.2,5)  ; KM
  (0,2.5,20) ; V1
  (0,10,50)  ; Q
  (0,4,20)   ; V2

$OMEGA
  0.2        ; CL
  0.2        ; VMAX
  0 FIX      ; KM
  0.1        ; V1
  0 FIX      ; Q
  0 FIX      ; V2

$SIGMA
  0.15      ; [P] Proportional

$SIMULATION (run_num)
$ESTIMATION METHOD=1 INTER PRINT=1 SIGL=10 NSIG=3 MSFO=./run_num.msf
