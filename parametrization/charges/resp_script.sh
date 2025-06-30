#!/bin/bash 

antechamber="$HOME/miniconda3/envs/AmberTools22/bin/antechamber"
respgen="$HOME/miniconda3/envs/AmberTools22/bin/respgen"
espgen="$HOME/miniconda3/envs/AmberTools22/bin/espgen"
resp="$HOME/miniconda3/envs/AmberTools22/bin/resp"

$antechamber -i "TTDimer.mol2" -fi "mol2" -o "TTDimer.ac" -fo "ac" -pf "yes" -nc 0
$respgen -i "TTDimer.ac" -o "TTDimer.respin1" -f "resp1" -a "TTDimer_Constant_Charges.dat"
$respgen -i "TTDimer.ac" -o "TTDimer.respin2" -f "resp2" -a "TTDimer_Constant_Charges.dat"
$espgen -i "TTDimer.gesp" -o "TTDimer.esp"
$resp -O -i "TTDimer.respin1" -o "TTDimer.respout1" -e "TTDimer.esp" -q "QIN" -t "TTDimer_qout_s1" -p "TTDimer_punch_s1" -s "TTDimer_esout_s1"
$resp -O -i "TTDimer.respin2" -o "TTDimer.respout2" -e "TTDimer.esp" -q "TTDimer_qout_s1" -t "TTDimer_qout_s2" -p "TTDimer_punch_s2" -s "TTDimer_esout_s2"
$antechamber -i "TTDimer.ac" -fi "ac" -o "TTDimer_RESP.mol2" -fo "mol2" -c "rc" -cf "TTDimer_qout_s2" -pf "yes" -at "amber"

