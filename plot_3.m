clear;

dist = 4;
n_u = 1;
n_t = 4;
sigma1 = sqrt(10^-(3.5));
sigma0 = sqrt(10^-11);
Kdb = 6;
loop = 1000;

%rho = [ linspace(0.1,0.9,9) linspace(0.90001,0.9999998,10) 0.9999999];
rho = linspace(0.001, 0.999999, 99);

S1 = zeros(1,length(rho));
S2 = zeros(1,length(rho));
S3 = zeros(1,length(rho));



for k = 1:loop
[Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb);

Pdb = 10; %dbm
P   = 10.^(Pdb/10);
theta = 0.04764;
[theta] = nonlinearEH (theta);



%optimal beamformer
w = conj(Hpl)/norm(Hpl);

Gamma = P*abs(w.'*Hpl)^2;

for i =1: length(rho)
%Data rate
R(i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(i))));

EHDC(i) = rho(i)*Gamma*(1-theta/(rho(i)*Gamma));

EHDC_linear(i) = rho(i)*Gamma;

[EHDC(i)] = inverse_nonlinearEH (EHDC(i),3.9 );
end
S1 = S1 + EHDC;
S3 = S3 + EHDC_linear;


theta1 = 0.00027;


for i =1: length(rho)
%Data rate
R1(i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(i))));

EHDC1(i) = rho(i)*Gamma*(1-theta1/(rho(i)*Gamma));

[EHDC1(i)] = inverse_nonlinearEH (EHDC1(i),3.9);
end
S2 = S2 + EHDC1;

k
end

EHDC = S1/loop;
EHDC1 = S2/loop;
EHDC_linear = S3/loop;

figure(4)
plot( R, EHDC_linear,'-.k',R, EHDC, '-r', R1,EHDC1,'--b')
grid on
xlabel('Rate (bits/channel use)')
ylabel('EH (mW/s)')
legend('input energy', 'non-linear EH', 'non-linear EH with AC computing')




