clear;
close all;
clc;

% for spacecraft with spin stabilization

I3 = .5*500*1.575^2/4;
I1 = .25*500*1.575^2/4 + 1/12*500*1.0668^2;
I2 = I1;
we3 = 8.1; % angular velocity rad/s
SRP = 1;

A = [(I2-I3)*we3/I1 0 0; 0 (I3-I2)*we3/I2 0; 0 0 0];
B = diag(ones(1, 3));
Q = [10 0 0; 0 10 0; 0 0 1];
R = 1;
C = Q;
D = 0;
[K, ~, ~] = lqr(A, B, Q, R);

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];
sys_cl = ss(Ac, B, C, D);
% c2d and d2c for continuous to discrete

t = 0:0.01:30;
r = zeros(length(t), 3);
r(1,:) = 1;
lsim(sys_cl, r, t);
grid on;

impulse(sys_cl)
%step(sys_cl)
