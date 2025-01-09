% Author: gjh
% date: 2024/11/26

%% clc
clc; clear; close all;

%% init
global M m l g sigma1 sigma2 Tp sita
M = 1;
omega = 1;
m = 0.5;
l = 1;
g = 9.81;
sigma1 = 3;
sigma2 = 3;
Tp = 0.5;
x_init1 = [5; -5];
x_init2 = [3; -3];
x_init3 = [1; -1];
sita_hat = 0;
sita = omega / M;
f_init1 = cat(1, x_init1, sita_hat, 0);
f_init2 = cat(1, x_init2, sita_hat, 0);
f_init3 = cat(1, x_init3, sita_hat, 0);

%% Plotting x1
figure;

tic;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init1);
plot(t, x(:, 1), 'b');
hold on;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init2);
plot(t, x(:, 1), 'r');
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init3);
plot(t, x(:, 1), 'g');
hold off;
title('Response of x1');
xlabel('Time (s)');
ylabel('x1');
legend('(x_1(0),x_2(0)) = (5, -5)', '(x_1(0),x_2(0)) = (3, -3)', '(x_1(0),x_2(0)) = (1, -1)');

%% Plotting x2
figure;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init1);
plot(t, x(:, 2), 'b');
hold on;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init2);
plot(t, x(:, 2), 'r');
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init3);
plot(t, x(:, 2), 'g');
hold off;
title('Response of x2');
xlabel('Time (s)');
ylabel('x2');
legend('(x_1(0),x_2(0)) = (5, -5)', '(x_1(0),x_2(0)) = (3, -3)', '(x_1(0),x_2(0)) = (1, -1)');

%% Plotting Control Input u
figure;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init1);
[~, U] = cellfun(@(t, x) Ap(t, x.'), num2cell(t), num2cell(x, 2), 'uni', 0);
plot(t, cell2mat(U), 'b');
hold on;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init2);
[~, U] = cellfun(@(t, x) Ap(t, x.'), num2cell(t), num2cell(x, 2), 'uni', 0);
plot(t, cell2mat(U), 'r');
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init3);
[~, U] = cellfun(@(t, x) Ap(t, x.'), num2cell(t), num2cell(x, 2), 'uni', 0);
plot(t, cell2mat(U), 'g');
hold off;
title('Control Input u');
xlabel('Time (s)');
ylabel('u');
legend('(x_1(0),x_2(0)) = (5, -5)', '(x_1(0),x_2(0)) = (3, -3)', '(x_1(0),x_2(0)) = (1, -1)');
toc;

%% Plotting sita_hat
figure;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init1);
plot(t, x(:, 3), 'b');
hold on;
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init2);
plot(t, x(:, 3), 'r');
[t, x] = ode45(@Ap, 0:0.01:0.6, f_init3);
plot(t, x(:, 3), 'g');
hold off;
title('Estimate of Uncertain Parameter sita_hat');
xlabel('Time (s)');
ylabel('sita_hat');
legend('(x_1(0),x_2(0)) = (5, -5)', '(x_1(0),x_2(0)) = (3, -3)', '(x_1(0),x_2(0)) = (1, -1)');

%% Function Ap
function [dx, u] = Ap(t, x)
global M m l g sigma1 sigma2 Tp sita
x1 = x(1);
x2 = x(2);
sita_hat = x(3);
dx = zeros(length(x), 1);

fai1 = -x1 + sita_hat * x2 + 1/2 * m * g * l * sin(x1 / M) - sigma1 * x2 / (Tp - t) - sigma1 * x1 / (Tp - t)^2;

if t <= Tp
    z1 = x1;
    alpha1 = -sigma1 * z1 / (Tp - t);
    z2 = x2 - alpha1;
    u = -sigma2 * z2 / (Tp - t) + fai1;
    sita_hat_dot = -2 * z2 * x2;
else
    u = 0;
    sita_hat_dot = 0;
end

x1_dot = x2;
x2_dot = u - sita * x2 - 1/2 * m * g * l * sin(x1 / M);

dx(1) = x1_dot;
dx(2) = x2_dot;
dx(3) = sita_hat_dot;
end