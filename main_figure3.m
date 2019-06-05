%MATLAB code to reproduce Fig.3 in our works on AC-enabled SWIPT research
%(perfect CSI)
%Written by Ha Vu Tran and Georges Kaddoum

clear;
dist = 4; %transmission distance
n_u = 1; %one user
n_t = 4; %four transmit antennas
sigma1 = sqrt(10^-(3.5)); %Gaussian noise
sigma0 = sqrt(10^-11); %Gaussian noise
Kdb = 6; %Rician factor
MEH = 3.9; %maximum of harvested energy - nonlinear model
loop = 10000; %the number of channel realizations
rho = linspace(0.00001, 0.99999, 99); %\rho - splitting factor

S1 = zeros(1,length(rho));
S2 = zeros(1,length(rho));
S3 = zeros(1,length(rho));



for k = 1:loop   
    [Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb); %generating channels
    P_0db = 10; %dbm
    P_0   = 10.^(P_0db/10); %convert dbm to mW
    
    %------------- DC computing -----------------------
    theta = 0.04764; %set threshold - nonlinear EH DC - this value is set according to the results of Fig. 6 in [10] 

    w = opt_beamformer(Hpl);  %Optimal beamformer

    Gamma = P_0*abs(w'*Hpl)^2; %Note that \psi = 0 because we are considering perfect CSI for this case
    
    %R-E region
    for i =1: length(rho)
        R(i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(i)))); %rate
        
        EHDC_noloss(i) = rho(i)*Gamma; % beamed energy at the user
        
        [hat_theta] = nonlinear_linear_threshold (theta); %mapping it to the threshold for linear EH 

        hat_EHDC(i) = rho(i)*Gamma - hat_theta; %subtracting energy used for DC computing

        [EHDC(i)] = nonlinearEH (hat_EHDC(i),MEH ); %mapping to nonlinear EH
    end

    S1 = S1 + EHDC;
    S3 = S3 + EHDC_noloss;

    %------------- AC computing -----------------------
    theta1 = 0.00027; %note that we use the AC directly for AC computing - this value is set according to the results of Fig. 6 in [10]
    
    %R-E region
    for i =1: length(rho)
        R1(i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(i)))); %rate

        hat_EHDC1(i) = rho(i)*Gamma*(1-theta1/(rho(i)*Gamma)); %linear EH

        [EHDC1(i)] = nonlinearEH (hat_EHDC1(i),3.9); %nonlinear EH
    end
    S2 = S2 + EHDC1;
k
end

EHDC = S1/loop;
EHDC1 = S2/loop;
EHDC_noloss = S3/loop;

figure(3)
plot( R, EHDC_noloss,'-.k',R, EHDC, '-r', R1,EHDC1,'--b')
grid on
xlabel('Rate (bit/s/Hz)')
ylabel('EH (mW/s)')
legend('EH without AD-to-DC conversion loss', 'non-linear EH with DC computing', 'non-linear EH with AC computing')




