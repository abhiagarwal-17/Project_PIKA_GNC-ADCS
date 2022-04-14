# Project_PIKA_GNC-ADCS
Princeton MAE undergraduates designing a scientific mission to the moon. This repository provides packages to do GNC/ ADCS analysis for the same. More details about project PIKA can be obtained by contacting the team. 

Use `PIKA_ADCS_scripts.m` to calculate external disturbance torques on a cylindrical spacecraft and to size reaction wheels when for a cylindrical 3-axis stabilized spacecraft. 

Use `PIKA_Actuator.m` to solve coupled differential equations and find the Euler parameters and angular spin rate of spacecraft when corrections are made for thruster misalignment from center of mass. Same file can then be used to justify calculations on attitude control thrusters. 

Use `PIKA_Stationkeeping.m` to estimate delta v values in low lunar orbit for stationkeeping. 

Use `PIKA_ADCS_ControlPlay.m` to simulate open-loop (uncontrolled) and closed-loop (LQR) dynamical models for the on-orbit spacecraft.
