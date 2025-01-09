% Author:gjh
% date：2024/11/26

%% clc
clc;clear;close all;

%% init

global sigma1 sigma2 sigma3 Tp sita

sigma2 = 5;
sigma1 = 5;
sigma3 = 5;
Tp = 5;
x_init1 = [1; -2; 1];
x_init2 = [-1; 2; -1];
x_init3 = [1; -1; -1];
sita_hat = 0;
sita = 0.5;
f_init1 = cat(1,x_init1,sita_hat);
f_init2 = cat(1,x_init2,sita_hat);
f_init3 = cat(1,x_init3,sita_hat);

options = odeset('RelTol',1e-6, 'AbsTol',1e-6);

tspan = 0:0.01:6;

figure

%% run1
tic
[t,x] = ode45(@Ap,tspan,f_init1, options);
plot(t,x(:,1))
hold on
%% run2
[t,x] = ode45(@Ap,tspan,f_init2, options);
plot(t,x(:,1))
%% run3
[t,x] = ode45(@Ap,tspan,f_init3, options);
plot(t,x(:,1))
hold off

figure

%% run1
[t,x] = ode45(@Ap,tspan,f_init1, options);
plot(t,x(:,2))
hold on
%% run2
[t,x] = ode45(@Ap,tspan,f_init2, options);
plot(t,x(:,2))
%% run3
[t,x] = ode45(@Ap,tspan,f_init3, options);
plot(t,x(:,2))
hold off

figure


%% run1
[t,x] = ode45(@Ap,tspan,f_init1, options);
plot(t,x(:,3))
hold on
%% run2
[t,x] = ode45(@Ap,tspan,f_init2, options);
plot(t,x(:,3))
%% run3
[t,x] = ode45(@Ap,tspan,f_init3, options);
plot(t,x(:,3))
hold off


figure


%% u
%% run1
[t,x] = ode45(@Ap,tspan,f_init1, options);
[~,U] = cellfun(@(t,x)Ap(t,x.'), num2cell(t), num2cell(x,2),'uni',0);
plot(t,cell2mat(U))
hold on
%% run2
[t,x] = ode45(@Ap,tspan,f_init2, options);
[~,U] = cellfun(@(t,x)Ap(t,x.'), num2cell(t), num2cell(x,2),'uni',0);
plot(t,cell2mat(U))
%% run3
[t,x] = ode45(@Ap,tspan,f_init3, options);
toc
[~,U] = cellfun(@(t,x)Ap(t,x.'), num2cell(t), num2cell(x,2),'uni',0);
plot(t,cell2mat(U))
hold off

%% f
function [dx,u] = Ap(t, x)

global sigma1 sigma2 sigma3 Tp sita

x1 = x(1);
x2 = x(2);
x3 = x(3);
sita_hat = x(4);
dx = zeros(length(x),1);

% Smoothing function
f1 = x1; 
f2 = 0; 
f3 = x3^2; 

% % Partial derivative
% 
% z1 = x1;
% alpha1 = -sigma1*z1/(Tp-t) - sita_hat*f1;
% z2 = x2 - alpha1;
% 
% tau1 = z1*f1;
% 
% dalpha1_x1 = -sigma1/(Tp-t) - sita_hat;
% tau2 = tau1 + z2*(f2 - dalpha1_x1*f1);
% dalpha1_sita_hat = - f1;
% dalpha1_t = -z1*sigma1/(Tp-t)^2;
% alpha2_z2 = -sigma1/(Tp-t) - (sigma1/(Tp-t)+sita_hat)*x1*x1;
% dalpha2_x1 = -1-(z1*x1+(sigma1/(Tp-t)+sita_hat)*z2*x1)- ...
%              x1*(2*x1+(sigma1/(Tp-t)+sita_hat)*z2)- ...
%              (sigma1/(Tp-t)+sita_hat)*sita_hat - sigma1/((Tp-t)^2) + ...
%              alpha2_z2*(sigma1/(Tp-t)+sita_hat);
% dalpha2_x2 = -(sigma1/(Tp-t)+sita_hat) + alpha2_z2;
% dalpha2_t = -sigma1*z2/((Tp-t)^2) - sigma1*(x2+x1*sita_hat)/((Tp-t)^2)- ...
%             2*sigma1*z1/((Tp-t)^3) - (sigma1*z2*x1^2)/((Tp-t)^2)- ...
%             2*sigma1*z1/(Tp-t)^3 + alpha2_z2*sigma1*z1/((Tp-t)^2);
% dalpha2_sita_hat = -z2*x1^2-(x2+sita_hat*x1) -(sigma1/(Tp-t)+sita_hat)*x1 + alpha2_z2*x1;
% 
% % Formula derivation
% 
% alpha2 = -sigma2*z2/(Tp-t) - z1 -...
%     sita_hat'*f2 + dalpha1_sita_hat*tau2 +...
%     dalpha1_x1*(x2+sita_hat'*f1) + dalpha1_t;
% 
% z3 = x3 - alpha2;
% 
% tau = z1*x1 - z2*dalpha1_x1*x1 + z3*(x3^2 - dalpha2_x1*x1);
% 
% fai2 = -z2 - sita_hat*x3^2 + ...
%     dalpha2_x1*(x2+sita_hat*x1) + dalpha2_x2*x3 +...
%     dalpha2_t + dalpha2_sita_hat*tau - ...
%     z2*dalpha1_sita_hat*(x3-dalpha2_x1*x1);

if 0 <= t && t < Tp
    % Partial derivative
    
    z1 = x1;
    alpha1 = -sigma1*z1/(Tp-t) - sita_hat*f1;
    z2 = x2 - alpha1;
    
    tau1 = z1*f1;
    
    dalpha1_x1 = -sigma1/(Tp-t) - sita_hat;
    tau2 = tau1 + z2*(f2 - dalpha1_x1*f1);
    dalpha1_sita_hat = - f1;
    dalpha1_t = -z1*sigma1/(Tp-t)^2;
    alpha2_z2 = -sigma1/(Tp-t) - (sigma1/(Tp-t)+sita_hat)*x1*x1;
    dalpha2_x1 = -1-(z1*x1+(sigma1/(Tp-t)+sita_hat)*z2*x1)- ...
                 x1*(2*x1+(sigma1/(Tp-t)+sita_hat)*z2)- ...
                 (sigma1/(Tp-t)+sita_hat)*sita_hat - sigma1/((Tp-t)^2) + ...
                 alpha2_z2*(sigma1/(Tp-t)+sita_hat);
    dalpha2_x2 = -(sigma1/(Tp-t)+sita_hat) + alpha2_z2;
    dalpha2_t = -sigma1*z2/((Tp-t)^2) - sigma1*(x2+x1*sita_hat)/((Tp-t)^2)- ...
                2*sigma1*z1/((Tp-t)^3) - (sigma1*z2*x1^2)/((Tp-t)^2)- ...
                2*sigma1*z1/(Tp-t)^3 + alpha2_z2*sigma1*z1/((Tp-t)^2);
    dalpha2_sita_hat = -z2*x1^2-(x2+sita_hat*x1) -(sigma1/(Tp-t)+sita_hat)*x1 + alpha2_z2*x1;
    
    % Formula derivation
    
    alpha2 = -sigma2*z2/(Tp-t) - z1 -...
        sita_hat'*f2 + dalpha1_sita_hat*tau2 +...
        dalpha1_x1*(x2+sita_hat'*f1) + dalpha1_t;
    
    z3 = x3 - alpha2;
    
    tau = z1*x1 - z2*dalpha1_x1*x1 + z3*(x3^2 - dalpha2_x1*x1);
    
    fai2 = -z2 - sita_hat*x3^2 + ...
        dalpha2_x1*(x2+sita_hat*x1) + dalpha2_x2*x3 +...
        dalpha2_t + dalpha2_sita_hat*tau - ...
        z2*dalpha1_sita_hat*(x3-dalpha2_x1*x1);
    u = -sigma3*z3/(Tp-t) + fai2;
    sita_hat_dot = tau;
else
    u = 0;
    sita_hat_dot = 0;
end

x1_dot = x2 + sita*x1;
x2_dot = x3;
x3_dot = u + sita*x3^2;

dx(1) = x1_dot;
dx(2) = x2_dot;
dx(3) = x3_dot;
dx(4) = sita_hat_dot;
end