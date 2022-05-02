%% Dynamical Modelling and State Space Realizations for the ADCS/GNC Subsystem 

close all 
clear vars

%Moments of Inertia calculated using parallel axis theorem

Ix = 283.3587;
Iy = 310.4956;
Iz = 593.8545;

%initial error [rad/s] --> 1 degree error deviation 
wx0 =    0.0017;
wy0 =    0.0017;
wz0 =    0.0017;


%initial angular error (5 deg)
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
%B1 = zeros(3,3);
B1 = eye(3,3);

%state transformation required to apply control input on angle rather than
%angular velocity

An = [Ao A11; A12 A21];

Bn = [B1; Bo];

Cn = eye(6,6);

%use the lqr method to stabilize the system 
M = 500; 
Q1 = Cn.'*M*Cn;
R = 50;

K2 = lqr(An,Bn,Q1,R);

%Discrete simulation of the open loop system

%timesteps
T = 6000;
%state vector and initial conditions
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.01;

for i = 1:6000
    u = 0;
    w1 = randn(6,1).*10^-15;
    w2 = [0;0;0;0;0;Iz^-1*1.25*10^-4];
    x(:,i+1) = x(:,i) + t*(An*x(:,i) + w2);
end

T = linspace(1,6000,6001);
figure()
plot(T,x(4,:),T,x(5,:),T,x(6,:));
title('Open Loop Dynamics');
legend('\theta_1','\theta_2','\theta_3','location','northeastoutside');
xlabel('Timesteps (s)');
ylabel('Angle Error (rad)');

for j = 1:5
%Discrete simulation of closed loop system
T = 6000;
x = zeros(6,T);
x(:,1) = [wx0;wy0;wz0;thex0;they0;they0];
t = 0.01;



%Stability envelope of +/- 0.1 degrees
for i = 1:6000 %use linspace with width = delT linspace(1,6000,600001);
    if x(4:6,i) >= 0.00174533 | x(4:6,i) <= -0.00174533
        Kcl = K2;
    elseif x(4:6,i) <= (0.00174533)/2 & x(4:6,i) <= -(0.00174533)/2
        Kcl = zeros(3,6);
    end
    
    Anew = An - Bn*Kcl;
    u(1:3,i) = -Kcl*x(:,i);
    w1 = randn(6,1).*10^-15;
    w2 = [0;0;0;0;0;Iz^-1*10^-4];
    x(:,i+1) = x(:,i) + t*(Anew*x(:,i) + w2);
end

T = linspace(1,6000,6001);
% figure()
% plot(T,x(4,:),T,x(5,:),T,x(6,:));
% yline(0.00174533,'r:');
% yline(-0.00174533,'r:');
% legend('\theta_1','\theta_2','\theta_3','location','northeastoutside');
% title('Closed Loop Dynamics');
% xlabel('Timesteps');
% ylabel('Angle Error (rad)');

I = [Ix 0 0; 0 Iy 0; 0 0 Iz];

%T = I^1*K2*x(:,:);

%Max and mean torques for the manuever
%T = I^-1*K2(:,1:3)*(Anew(1:3,1:3)*x(1:3,:) + Bn(4:6,:)*K2(:,1:3)*x(1:3,:));
%T = I^-1*K2(:,1:3)*(An(1:3,1:3)*x(1:3,:) + Bn(4:6,:)*K2(:,1:3)*x(1:3,:));
%T = I^-1*Bn(1:3,1:3)*K2(:,1:3)*x(1:3,:);

%To = u(1:3,:);
To = u;

%T = (An(1:3,1:3)*x(1:3,:) + Bn(4:6,:)*K2(:,1:3)*x(1:3,:));

Tmax = max(To,[],'all');
Tmean = mean(To,'all');
%for the future: solve non-linear equations with linear controller to
%validate controller in the non-linear regime

h = Tmean*3.154e+7;
%h = Tmean*4.734e+7;

%Number of saturations in 1 years
timesat = h/4;
timesatt(1,j) = timesat;
Tomax(1,j) = Tmax;
Tomean(1,j) = Tmean;

end

T = linspace(1,6000,6001);
figure()
h = plot(T,x(4,:),T,x(5,:),T,x(6,:));
h(4) = yline(0.00174533,'r:');
h(5) = yline(-0.00174533,'r:');
h(6) = yline(0.0087,'g:');
h(7) = yline(-0.0087,'g:');

legend([h(1) h(2) h(3) h(4) h(6)],'\theta_1','\theta_2','\theta_3',...
    'controller stability envelope (\pm 0.1^o)','mission pointing requirement (\pm 0.5^o)','location','northeastoutside');
title('Closed Loop Dynamics');
xlabel('Timesteps (s)');
ylabel('Angle Error (rad)');

meansat = mean(timesatt);
meanmax = mean(Tomax);
meanmean = mean(Tomean);
