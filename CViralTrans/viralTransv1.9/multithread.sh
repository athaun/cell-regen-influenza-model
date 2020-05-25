!/bin/bash
parameters=( 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3 )
# sadly, the for loop solution did not work :(
# TODO implement for loop
  (echo "[SHELL] Running test with parameter ${parameters[0]}"; nvcc ViralTransmission.cu -o program0.out && ./program0.out 0.003; echo "[SHELL] test with parameter ${parameters[0]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[1]}"; nvcc ViralTransmission.cu -o program1.out && ./program1.out 0.01; echo "[SHELL] test with parameter ${parameters[1]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[2]}"; nvcc ViralTransmission.cu -o program2.out && ./program2.out 0.03; echo "[SHELL] test with parameter ${parameters[2]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[3]}"; nvcc ViralTransmission.cu -o program3.out && ./program3.out 0.1; echo "[SHELL] test with parameter ${parameters[3]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[4]}"; nvcc ViralTransmission.cu -o program4.out && ./program4.out 0.3; echo "[SHELL] test with parameter ${parameters[4]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[5]}"; nvcc ViralTransmission.cu -o program5.out && ./program5.out 1.0; echo "[SHELL] test with parameter ${parameters[5]} is done.")
  (echo "[SHELL] Running test with parameter ${parameters[6]}"; nvcc ViralTransmission.cu -o program6.out && ./program6.out 3.0; echo "[SHELL] test with parameter ${parameters[6]} is done.")
#  ps

echo "[SHELL]"
ps -AF | grep "program.out"
