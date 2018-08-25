%   Optimal Control Problem

clc;
clear all;
close all;

time = 100;
dt = 0.01;
t = 0:dt:time;
npoints = time/dt + 1;

% Compartments %

% Case of optimal control %
s = zeros(npoints, 1);
in = zeros(npoints, 1);
r = zeros(npoints, 1);
cost = zeros(npoints, 1);

% Case of no vaccination %
s2 = zeros(npoints, 1);
in2 = zeros(npoints, 1);
r2 = zeros(npoints, 1);
cost2 = zeros(npoints, 1);

% Case of constant vaccination %
s3 = zeros(npoints, 1);
in3 = zeros(npoints, 1);
r3 = zeros(npoints, 1);
cost3 = zeros(npoints, 1);

% Initial Conditions and Constants %
s(1,1) = 0.95;
in(1,1) = 0.05;

s2(1,1) = 0.95;
in2(1,1) = 0.05;

s3(1,1) = 0.95;
in3(1,1) = 0.05;

k = 0.5;
p = 0.4;
N = 1;
v = 0.3796479517*exp(-0.1242921619*t);
vv3 = 0.06;
u = 0.1;
A = 0.03089708305;

for i=2:npoints
   s(i, 1) = s(i-1, 1) + ((-k*p/N).*(in(i-1, 1).*s(i-1, 1)) - v(1, i-1).*s(i-1, 1)).*dt;
   in(i, 1) = in(i-1, 1) + ((k*p/N).*(in(i-1, 1).*s(i-1, 1)) - u.*in(i-1, 1)).*dt;
   r(i, 1) = r(i-1, 1) + (u.*in(i-1, 1) + v(1, i-1).*s(i-1, 1)).*dt;
   
   s2(i, 1) = s2(i-1, 1) + ((-k*p/N).*(in2(i-1, 1).*s2(i-1, 1))).*dt;
   in2(i, 1) = in2(i-1, 1) + ((k*p/N).*(in2(i-1, 1).*s2(i-1, 1)) - u.*in2(i-1, 1)).*dt;
   r2(i, 1) = r2(i-1, 1) + (u.*in2(i-1, 1)).*dt;
   
   s3(i, 1) = s3(i-1, 1) + ((-k*p/N).*(in3(i-1, 1).*s3(i-1, 1)) - vv3.*s3(i-1, 1)).*dt;
   in3(i, 1) = in3(i-1, 1) + ((k*p/N).*(in3(i-1, 1).*s3(i-1, 1)) - u.*in3(i-1, 1)).*dt;
   r3(i, 1) = r3(i-1, 1) + (u.*in3(i-1, 1) + vv3.*s3(i-1, 1)).*dt;
   
   cost(i, 1) = in(i-1, 1) + cost(i-1, 1) + (A/2)*(v(1, i-1).^2);
   cost2(i, 1) = in2(i-1, 1) + cost2(i-1, 1);
   cost3(i, 1) = in3(i-1, 1) + cost3(i-1, 1) + (A/2)*(vv3.^2);
end

figure(1);
plot(t, s, 'b');
hold on;
plot(t, in, 'r');
hold on;
plot(t, r, 'g');
xlabel('Time (Units)');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
title('SIR Model with optimal control');

figure(2);
plot(t, in2, 'r', 'LineStyle', '--', 'LineWidth', 3);
hold on;
plot(t, in, 'g', 'LineStyle', '-.', 'LineWidth', 3);
grid on;
xlabel('Time (Units)');
ylabel('Population');
legend('I without control', 'I with control');
title('Comparison of I');


figure(3);
plot(t, s2, 'r', 'LineStyle', '--', 'LineWidth', 3);
hold on;
plot(t, s, 'g', 'LineStyle', '-.', 'LineWidth', 3);
grid on;
xlabel('Time (Units)');
ylabel('Population');
legend('S without control', 'S with control');
title('Comparison of S');

figure(4);
plot(t, r2, 'r', 'LineStyle', '--', 'LineWidth', 3);
hold on;
plot(t, r, 'g', 'LineStyle', '-.', 'LineWidth', 3);
grid on;
xlabel('Time (Units)');
ylabel('Population');
legend('R without control', 'R with control');
title('Comparison of R');


figure(5);
plot(t, cost2, 'r', 'LineStyle', '--', 'LineWidth', 3);
hold on;
plot(t, cost, 'g', 'LineStyle', '-.', 'LineWidth', 3);
grid on;
xlabel('Time (Units)');
title('Objective Functional Comparison');
legend('Case without control', 'Optimal Case');

figure(6);
plot(t, 0.3796479517*exp(-0.1242921619*t), 'r', 'LineStyle', '--', 'LineWidth', 3);
grid on;
xlabel('Time (Units)');
ylabel('u(t)');
title('Control');
