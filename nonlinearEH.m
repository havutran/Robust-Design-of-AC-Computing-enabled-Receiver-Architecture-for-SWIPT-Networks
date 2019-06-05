function [theta_new] = nonlinearEH (theta,M)
theta=theta/1000;
M=M/1000;

a = 1500;
b = 0.0022;

%a = 6.400;
%b = 0.0022;

theta_new = (M/(1+exp(-a*(theta-b))) - M/(1+exp(a*b)))/(1-1/(1+exp(a*b)));

theta_new = theta_new*1000;
end

