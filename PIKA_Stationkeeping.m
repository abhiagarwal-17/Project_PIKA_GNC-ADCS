%% Constants
clc; 
clear; 
m_moon = 7.347e22; %kg
R_moon = 1.737e6; %m
G = 6.67e-11; %SI Units
mu_moon = G * m_moon; 
h_orbit = 50 * 1000; %m (altitude of circular orbit) 
R_orbit = R_moon + h_orbit; 
v_orbit = sqrt(mu_moon/R_orbit); 

frequency = 27; % days

% SRP values
% all planets plus asteroid tempelates
% SRP: spherical with 0.0141 m^2/kg estimate of area/mass
%% Calculation of orbital velocity change
% parameters of orbit after 25 days 
perigee = 45.98 * 1000 + R_moon; 
apogee = 51.42 * 1000 + R_moon; 

perigee = 30 * 1000 + R_moon; 
apogee = 70 * 1000 + R_moon; 
a = (apogee + perigee)/2; 

v_perigee = sqrt(mu_moon * (2/perigee - 1/a)); 
v_apogee = sqrt(mu_moon * (2/apogee - 1/a)); 

% case 1 (fire at perigee and then at apogee)
apogee_ideal = R_orbit; 
a_ideal_c1 = (apogee_ideal + perigee)/2; 

v_perigee_ideal_c1 = sqrt(mu_moon * (2/perigee - 1/a_ideal_c1));
v_apogee_ideal_c1 = sqrt(mu_moon * (2/apogee_ideal - 1/a_ideal_c1));

dv1_c1 = v_perigee_ideal_c1 - v_perigee; 
dv2_c1 = v_apogee_ideal_c1 - v_orbit;  

total_c1 = abs(dv1_c1) + abs(dv2_c1); 

total_c1 = (365/frequency) * total_c1; 


% case 2 (fire at apogee and then at perigee) 
perigee_ideal = R_orbit; 
a_ideal_c2 = (perigee_ideal +apogee)/2; 

v_apogee_ideal_c2 = sqrt(mu_moon * (2/apogee - 1/a_ideal_c2));
v_perigee_ideal_c2 = sqrt(mu_moon * (2/perigee_ideal - 1/a_ideal_c2)); 

dv1_c2 = v_apogee_ideal_c2 - v_apogee; 
dv2_c2 = v_perigee_ideal_c2 - v_orbit; 

total_c2 = abs(dv1_c2) + abs(dv2_c2); 

total_c2 = (365/frequency) * total_c2; 

