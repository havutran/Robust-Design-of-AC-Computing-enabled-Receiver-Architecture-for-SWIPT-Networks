clear;

dist = 5;
n_u = 1;
n_t = 6;
sigma1 = sqrt(10^-5);
sigma0 = sqrt(10^-7);
Kdb = 6;

[Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb);

Pdb = 15; %dbm
P   = 10.^(Pdb/10);
theta = 2;
epsilon = linspace(1,21,20);
[epsilon_new] = nonlinearEH (epsilon);


%optimal beamformer
w = conj(Hpl)/norm(Hpl);

Gamma = P*abs(w.'*Hpl)^2;

rho = linspace(0,1,11);

for i =1: length(rho)
%optimal phi
%phi(i) = 1/(theta*epsilon_new(i)) + 1/(epsilon_new(i)^2);

%optimal rho
%rho(i) = epsilon_new(i)/(Gamma*(1-phi(i)));

%Data rate
R(i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(i))));

EH(i) = rho(i)*Gamma/2;
end


figure(2)
plot(R,EH, '-r')
grid on
xlabel('Rate')
ylabel('EH (mW/s)')
legend('')




