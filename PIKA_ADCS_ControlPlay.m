%% State Space Realizations

%% Dynamical Modelling and State Space Realizations for the ADCS/GNC Subsystem 

close all 
clear vars

%Moments of Inertia calculated using parallel axis theorem

Ix = 209.5303;
Iy = 229.4092;
Iz = 438.9395;

%initial error [rad/s] --> 1 degree error deviation 
wx0 =    0.0017;
wy0 =    0.0017;
wz0 =    0.0017;

%5000 for 0.0017 rad/s
thex0 = 0.0873;
they0 = 0.0873;
thez0 = 0.0873;

%linearize dynamics for angular velocity
Ao = [0 (-wz0*(Iz-Iy))/Ix (-wy0*(Iz-Iy))/Ix; (-wz0*(Ix-Iy))/Iy 0 (-wx0*(Ix-Iz))/Iy; (-wy0*(Iy-Iz))/Iz (-wx0*(Iy - Ix))/Iz 0 ];
Bo = eye(3); 
Do = 0;
Co = [1 1 1];

A11 = zeros(3,3);
A21 = zeros(3,3);
A12 = eye(3,3);
B1 = zeros(3,3);

%state transformation required to apply control input on angle rather than
%angular velocity
An = [A11 A12; Ao A21];
Bn = [B1; Bo];
Cn = [ones(1,3) ones(1,3)];



%use the lqr method to stabilize the system 
M = 100; 
Q1 = Cn.'*M*Cn;
R = 10;
K2 = lqr(An,Bn,Q1,R);


%Discrete simulation of the open loop system

%timesteps
T = 6000;
%state vector and initial conditions
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.002;

for i = 1:6000
    u = 0;
    w = randn(6,1).*10^-2;
    x(:,i+1) = x(:,i) + t*(An*x(:,i) + w);
end

T = linspace(1,6000,6001);
figure()
plot(T,x(4,:),T,x(5,:),T,x(6,:));
title('Open Loop Dynamics');
legend('\theta_1','\theta_2','\theta_3','location','northeastoutside');
xlabel('Timesteps');
ylabel('Angle Error (rad)');

%Discrete simulation of closed loop system
T = 6000;
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.002;

%Stability envelope of +/- 0.1 degrees
for i = 1:6000
    if x(1:3,i) >= 0.00174533 | x(1:3,i) <= -0.00174533
        Kcl = K2;
    elseif x(1:3,i) <= (0.00174533)/2 & x(1:3,i) <= -(0.00174533)/2
        Kcl = 0;
    end
    
    Anew = An - Bn*K2;
    u = -Kcl*x(:,i);
    w = randn(6,1).*10^-2;
    x(:,i+1) = x(:,i) + t*(Anew*x(:,i) + Bn*u + w);
end

T = linspace(1,6000,6001);
figure()
plot(T,x(4,:),T,x(5,:),T,x(6,:));
yline(0.00174533,'r:');
yline(-0.00174533,'r:');
legend('\theta_1','\theta_2','\theta_3','location','northeastoutside');
title('Closed Loop Dynamics');
xlabel('Timesteps');
ylabel('Angle Error (rad)');

I = [Ix 0 0; 0 Iy 0; 0 0 Iz];

%T = I^1*K2*x(:,:);

%Max and mean torques for the manuever
T = I^-1*K2(:,1:3)*(Anew(1:3,1:3)*x(1:3,:) + Bn(4:6,:)*K2(:,1:3)*x(1:3,:));
Tmax = max(T,[],'all');
Tmean = mean(T,'all');

%for the future: solve non-linear equations with linear controller to
%validate controller in the non-linear regime

h = Tmean*4.734e+7;

%Number of saturations in 1.5 years
timesat = h/4;


