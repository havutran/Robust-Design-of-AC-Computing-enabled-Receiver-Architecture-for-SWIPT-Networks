function [theta_new] = nonlinear_linear_threshold (theta)
theta=theta/1000;
%M = 0.024;
%a = 150;
%b = 0.014;
M = 3.9/1000;
a = 1500;
b = 0.0022;
theta_new = abs(b - 1/a.*log ( (exp(a*b)*(M - theta))./(exp(a*b)*theta + M) ));
theta_new = theta_new*1000;
end