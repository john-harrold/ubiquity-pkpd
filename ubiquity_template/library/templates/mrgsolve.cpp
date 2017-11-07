#
# Notes: 
#  how are variables in the dataset (covariates) implemented
#  how are covariates applied to initial conditions
#  what optional formats are there for omega and sigma? 
#  how do we define the ETA secondary parameters
#    - does the colon at the end of the block represent a comment?
#    - can we have 0 off diags with the eta matrix?
#
$PARAM
<SYSTEM_PARAMS>


$INIT
<STATE_INIT>

$MAIN
<IIV_SP>

<SSP>

<STATE_IC>

<DSP>

$OMEGA
#$OMEGA @annotated @block
#ECL: 0.09 : ETA on clearance
#EVC: 0.001 0.19 : ETA on volume
#EKA: 0.001 0.001 0.45 : ETA on absorption rate constant


$SIGMA

$ODE
<ODES>
