#!/bin/bash
parameters=( 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3 )
# sadly, the for loop solution did not work :(
# TODO implement for loop
  (echo "[SHELL] Running test with parameter ${parameters[0]}"; python3 2DProbabilityViralTransmission.py 0.003; echo "[SHELL] test with parameter ${parameters[0]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[1]}"; python3 2DProbabilityViralTransmission.py 0.01; echo "[SHELL] test with parameter ${parameters[1]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[2]}"; python3 2DProbabilityViralTransmission.py 0.03; echo "[SHELL] test with parameter ${parameters[2]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[3]}"; python3 2DProbabilityViralTransmission.py 0.1; echo "[SHELL] test with parameter ${parameters[3]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[4]}"; python3 2DProbabilityViralTransmission.py 0.3; echo "[SHELL] test with parameter ${parameters[4]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[5]}"; python3 2DProbabilityViralTransmission.py 1; echo "[SHELL] test with parameter ${parameters[5]} is done.") &
  (echo "[SHELL] Running test with parameter ${parameters[6]}"; python3 2DProbabilityViralTransmission.py 3; echo "[SHELL] test with parameter ${parameters[6]} is done.") &
#  ps

echo "[SHELL]"
ps -AF | grep "python3 2DProbabilityViralTransmission.py"
