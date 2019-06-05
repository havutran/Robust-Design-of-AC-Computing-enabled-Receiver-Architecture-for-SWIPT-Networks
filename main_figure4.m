%MATLAB code to reproduce Fig.4 in our works on AC-enabled SWIPT research
%(imperfect CSI)
%Written by Ha Vu Tran and Georges Kaddoum

clear;
dist = 4; %transmission distance
n_u = 1; %one user
n_t = 4; %six transmit antennas
sigma1 = sqrt(10^-(3.5)); %Gaussian noise
sigma0 = sqrt(10^-11); %Gaussian noise
Kdb = 6; %Rician factor
MEH = 3.9; %the maximum of harvested energy - nonlinear model

P_0db = linspace(12,17,10); %dbm
P_0   = 10.^(P_0db./10); %convert dbm into mW

epsilon = 0.2; %Threshold of EH DC
[hat_epsilon] = nonlinear_linear_threshold (epsilon); %mapping it to the threshold for linear EH 

psi = [0 0.02 0.06]; %Error factor

S1 = zeros(length(psi),length(P_0));

loop = 10000;

    %------------- AC computing -----------------------

theta = 0.00027; % this value is set according to the results of Fig. 6 in [10]
for l = 1:loop
    [Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb); %generating channels

    w = opt_beamformer(Hpl);  %Optimal beamformer

    for k = 1:length(psi)
        for i =1: length(P_0)
    
        Gamma(k) = P_0(i)*(1-sqrt(psi(k)))^2*abs(w'*Hpl)^2;

        phi(k,i) = hat_epsilon/(theta+hat_epsilon);  %optimal phi

        rho(k,i) = hat_epsilon/(Gamma(k)*(phi(k,i))); %optimal rho

        R(k,i) = log2 (1 + Gamma(k)/(sigma0^2 + sigma1^2/(1-rho(k,i))));   %Data rate
        end
    end

    S1 = S1 + R;
l
end

R = S1/loop;



S2 = zeros(length(psi),length(P_0));

    %------------- DC computing -----------------------
    
theta2 = 0.04764; % this value is set according to the results of Fig. 6 in [10] 
[theta2] = nonlinear_linear_threshold (theta2);

for l = 1:loop
    
    [Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb); %generating channels

    w = opt_beamformer(Hpl);  %Optimal beamformer

    for k = 1:length(psi)
        for i =1: length(P_0)
    
        Gamma2(k) = P_0(i)*(1-sqrt(psi(k)))^2*abs(w'*Hpl)^2;

        phi(k,i) = hat_epsilon/(theta2+hat_epsilon) ;   %optimal phi

        rho(k,i) = hat_epsilon/(Gamma2(k)*(phi(k,i)));  %optimal rho

        R2(k,i) = log2 (1 + Gamma2(k)/(sigma0^2 + sigma1^2/(1-rho(k,i))));   %Data rate
        end
    end

    S2 = S2 + R2;
l
end

R2 = S2/loop;

figure(4)
plot(P_0db, R(1,:), '-b', P_0db, R(2,:), '-or', P_0db, R(3,:), '-*k', P_0db, R2(1,:), '--b', P_0db, R2(2,:), '--or', P_0db, R2(3,:), '--*k')
grid on
ylabel('Rate (bit/s/Hz)')
xlabel('P_0 (dBm)')
legend('AC computing - \psi = 0','AC computing - \psi = 0.02','AC computing - \psi = 0.06', 'DC computing - \psi = 0','DC computing - \psi = 0.02','DC computing - \psi = 0.06')



























