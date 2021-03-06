%% Constants 
clc; 
clear; 
mass = 811; % kg
R = 1.575/2; % m (radius of the ESPA ring) 
H = 1.524; % m (height of the ESPA ring) 

Iz = 0.5 * mass * R^2; 
Ix = 0.25 * mass * R^2 + (1/12) * mass * H^2; 
Iy = 0.25 * mass * R^2 + (1/12) * mass * H^2; 

% Iz = 593.9; 
% Ix = 283.4; 
% Iy = 310.5; 

Thr = 490; % N (thrust exerted by main thruster) 

offset = deg2rad(0.3); % deg (angle by which the center of mass is offset)

delta_theta = deg2rad(1); 

total_time = 12.74 * 60; % seconds (total time of the burn)
%% Calculating initial estimate of spin velocity
Thrz = Thr * cos(offset); 
Thry = Thr * sin(offset); 

T_x = Thry * (H/2);  

omega_z = (T_x  * total_time)/ (Iz * delta_theta); 


%% Differential Solver
omega_z = 8; 
omega_0 = [0; 0; omega_z];
time = linspace(0, total_time, 200);
[t, w] = ode45(@vdp2, time, omega_0, [], Ix, Iy, Iz, T_x);
figure(1)
plot(t,w(:,3))
title('\omega_z v/s t');
xlabel('Time t (seconds)');
ylabel('\omega (rad/s)');

figure(2)
plot(t, w(:,1),t,w(:,2)); 
title('\omega v/s t');
xlabel('Time t (seconds)');
ylabel('\omega (rad/s)');
legend('w_1','w_2');

%%
broad_time = [0];
phi = [0];
theta = [0];
psi = [0];

K = length(w); 
for i = 1:K
 t_const = total_time/K;
 curr_time = linspace((i-1)*t_const, i*t_const, 10);
 init_cond = [psi(end), theta(end), phi(end)];

 [t, s] = ode45(@vdp3, curr_time, init_cond,[],[w(i, 1), w(i, 2), w(i, 3)]);

 broad_time = cat(1, broad_time, t);

 psi = cat(1, psi, s(:,1));
 theta = cat(1, theta, s(:,2));
 phi = cat(1, phi, s(:,3));
 
 disp(i)
end

%% taking modulo
phi = mod(phi, 2*pi); 
theta = mod(theta, 2*pi); 
psi = mod(psi, 2*pi); 


%% switching values close to 360 degrees
for i = 1:length(psi)
     if psi(i) > pi
         psi(i) = psi(i) - 2*pi; 
     end
     if theta(i) > pi
         theta(i) = theta(i) - 2*pi; 
     end
     if phi(i) > pi
         phi(i) = phi(i) - 2*pi; 
     end
    
end

%% plotting

figure(3)
plot(broad_time,rad2deg(psi));
title('Euler Angle: \psi');     
xlabel('Time t (seconds)');
ylabel('\psi (deg)');

figure(4)
plot(broad_time,rad2deg(theta));
title('Euler Angle: \theta');
xlabel('Time t (seconds)');
ylabel('\theta (deg)');

figure(5)
plot(broad_time,rad2deg(phi));
title('Euler Angle: \phi');
xlabel('Time t (seconds)');
ylabel('\phi (deg)');


%% Actuator Sizing
omega_z = 5.5; 
t = 0.02; %seconds

T = 4 * 52 * sin(deg2rad(47.22)) * (1.575/2); 

num_needed = sqrt((Iz * omega_z)/(t * T));

time_needed = num_needed * 10/1000; 

%% Functions
% y1 = psi
% y2 = theta
% y3 = phi

function dwdt = vdp2(t, w, Ix, Iy, Iz, T_x)
I1 = Ix;
I2 = Iy;
I3 = Iz;
dwdt = [(-1 * (I3 - I2) * w(2) * w(3) + T_x)/I1; (-1*(I1 - I3)*w(3)*w(1))/I2; 0];
end

function dydt = vdp3(t, y, omega)
dydt = [omega(2) * (sin(y(3))/cos(y(2))) + omega(3) * (cos(y(3))/cos(y(2))); omega(2) * cos(y(3)) - omega(3) * sin(y(3)); omega(1) + omega(2) * (sin(y(3)) * tan(y(2))) + omega(3) * (cos(y(3)) * tan(y(2)))];
end




