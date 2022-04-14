%% State Space Realizations


close all 
clear vars

%let's switch to a state space realization --> this will allows us to
%create a system with greater stability margins and state estimation 


   
Ix = 209.5303;
Iy = 229.4092;
Iz = 438.9395;

syms wx wy wz L1 L2 L3

wxd  = (-wy*wz*(Iz-Iy) + L1)/Ix;
wyd = (-wx*wz*(Ix-Iz) + L2)/Iy;
wzd = (-wx*wy*(Iy-Ix) + L3)/Iz;

Asym1 = jacobian([wxd 0 0],[wx,wy,wz]);
Asym2 = jacobian([0 wyd 0],[wx,wy,wz]);
Asym3 = jacobian([0 0 wzd],[wx,wy,wz]);


test = [wx 0 0; 0 wy 0; 0 0 wz];
 
%initial condition --> 1 degree error deviation 
wx0 = 0.0174533;
wy0 = 0.0174533;
wz0 = 0.0174533;

% wx0 = 0.00872665; %0.5 deg
% wy0 = 0.00872665;
% wz0 = 0.00872665;
% 
% wx0 = 0.00523599; %0.3 deg
% wy0 = 0.00523599;
% wz0 = 0.00523599;

%0.0174533; 1 deg

A = [0 (-wz0*(Iz-Iy))/Ix (-wy0*(Iz-Iy))/Ix; (-wz0*(Ix-Iy))/Iy 0 (-wx0*(Ix-Iz))/Iy; (-wy0*(Iy-Iz))/Iz (-wx0*(Iy - Ix))/Iz 0 ];
B = eye(3); 
%[1; 1; 1];
D = 0;
C = [1 1 1];

%open-loop system
%clearly unstable
sys2 = ss(A,B,C,D);
figure()
nyquist(sys2);
title('Open loop Nyquist Diagram');
ylabel('Angle (rad)');
opt = stepDataOptions('StepAmplitude',10^-4);
step(sys2,opt)
title('Open Loop Step Diagram');
ylabel('Angle (rad)');

%use the lqr method to stabilize the system 
M = 100; 
%can I accept > 0 steady state error and propose a steady state range
%50000;
Q1 = C.'*M*C;
R = 500;
K2 = lqr(A,B,Q1,R); %--> %turn into a disturbance rejection --> inject dsiturbance
%reference tracking so that out satellite will always try to adjust to some
%reference value, probably determined from a sensor. 

%Try a numerical integration


T = 10000;
x = zeros(3,T);
x(:,1) = [wx0;wy0;wz0];
t = 0.002;

for i = 1:10000
    u = 0;
    w = randn(3,1);
    x(:,i+1) = x(:,i) + t*(A*x(:,i) + w);
end

%open loop
T = linspace(1,10000,10001);
figure()
plot(T,x(1,:),T,x(2,:),T,x(3,:));
title('Open Loop Dynamics');
legend('\omega_1','\omega_2','\omega_3','location','northeastoutside');
xlabel('Timesteps');
ylabel('Angle Error (rad)');


%added control ennvelope %cone in a cone --> emphasize on meeting mission
%requirements --> could start with how analysis affects the main mission
%requirements

T = 100000;
x = zeros(3,T);
x(:,1) = [wx0;wy0;wz0];

for i = 1:100000
    if x(:,i) >= 0.00174533 | x(:,i) <= -0.00174533
        Kcl = K2;
    elseif x(:,1) <= (0.00174533)/2 & x(:,1) <= -(0.00174533)/2
        Kcl = 0;
    end
    
    Anew = A - B*K2;
    u = -Kcl*x(:,i);
    w = randn(3,1)*10^-5;
    x(:,i+1) = x(:,i) + t*(Anew*x(:,i) + B*u + w);
end

%for the future: solve non linear equations with linear solution

T = linspace(1,100000,100001);
figure()
plot(T,x(1,:),T,x(2,:),T,x(3,:));
yline(0.0174533,'r:');
yline(-0.0174533,'r:');
legend('\omega_1','\omega_2','\omega_3','location','northeastoutside');
title('Closed Loop Dynamics');
xlabel('Timesteps');
ylabel('Angle Error (rad)');


I = [Ix 0 0; 0 Iy 0; 0 0 Iz];

%T = I^1*K2*x(:,:); %why do I have to scale this by I^-1? --> units??
T = I^-1*K2*(Anew*x + B*K2*x(:,:));

%T = I*K2*x(:,:);

% figure()
% plot(linspace(1,100000,100001),cumtrapz(x(1,:)))

%% Dynamical Modelling and State Space Realizations for the ADCS/GNC Subsystem 


close all 
clear vars

%Moments of Inertia calculated using parallel axis theorem
Ix = 209.5303;
Iy = 229.4092;
Iz = 438.9395;

%initial error [rad/s] --> 1 degree error deviation 
wx0 = 0.05;
wy0 = 0.05;
wz0 = 0.05;
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
M = 10; 
Q1 = Cn.'*M*Cn;
R = 100;
K2 = lqr(An,Bn,Q1,R);


%Discrete simulation of the open loop system

%timesteps
T = 10000;
%state vector and initial conditions
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.002;

for i = 1:10000
    u = 0;
    w = randn(6,1);
    x(:,i+1) = x(:,i) + t*(An*x(:,i) + w);
end

T = linspace(1,10000,10001);
figure()
plot(T,x(4,:),T,x(5,:),T,x(6,:));
title('Open Loop Dynamics');
legend('\theta_1','\theta_2','\theta_3','location','northeastoutside');
xlabel('Timesteps');
ylabel('Angle Error (rad)');

%Discrete simulation of closed loop system
T = 10000;
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.002;

%Stability envelope of +/- 0.1 degrees
for i = 1:10000
    if x(1:3,i) >= 0.00174533 | x(1:3,i) <= -0.00174533
        Kcl = K2;
    elseif x(1:3,i) <= (0.00174533)/2 & x(1:3,i) <= -(0.00174533)/2
        Kcl = 0;
    end
    
    Anew = An - Bn*K2;
    u = -Kcl*x(:,i);
    w = randn(6,1)*10^-5;
    x(:,i+1) = x(:,i) + t*(Anew*x(:,i) + Bn*u + w);
end

T = linspace(1,10000,10001);
figure()
plot(T,x(4,:),T,x(5,:),T,x(6,:));
yline(0.0174533,'r:');
yline(-0.0174533,'r:');
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
