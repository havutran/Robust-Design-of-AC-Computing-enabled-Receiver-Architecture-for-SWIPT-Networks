clear;

dist = 4;
n_u = 1;
n_t = 6;
sigma1 = sqrt(10^-(3.5));
sigma0 = sqrt(10^-11);
Kdb = 6;


Pdb = linspace(10,15,10); %dbm
P   = 10.^(Pdb./10);
theta = 0.00027;
epsilon = [0.2 0.4 0.6];
[epsilon_new] = nonlinearEH (epsilon);
S1 = zeros(length(epsilon),length(P));

loop = 10000;

for l = 1:loop
[Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb);

%optimal beamformer
w = conj(Hpl)/norm(Hpl);


for k = 1:length(epsilon)
for i =1: length(P)
    
Gamma = P(i)*abs(w.'*Hpl)^2;

%optimal phi
phi(k,i) = epsilon_new(k)/(theta+epsilon_new(k)) ;

%optimal rho
rho(k,i) = epsilon_new(k)/(Gamma*(phi(k,i)));

%Data rate
R(k,i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(k,i))));
end
end

S1 = S1 + R;

end

R = S1/loop;



theta2 = 0.04764;
[theta2] = nonlinearEH (theta2);
S2 = zeros(length(epsilon),length(P));

loop = 10000;

for l = 1:loop
[Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb);

%optimal beamformer
w = conj(Hpl)/norm(Hpl);


for k = 1:length(epsilon)
for i =1: length(P)
    
Gamma = P(i)*abs(w.'*Hpl)^2;

%optimal phi
phi(k,i) = epsilon_new(k)/(theta2+epsilon_new(k)) ;

%optimal rho
rho(k,i) = epsilon_new(k)/(Gamma*(phi(k,i)));

%Data rate
R2(k,i) = log2 (1 + Gamma/(sigma0^2 + sigma1^2/(1-rho(k,i))));
end
end

S2 = S2 + R2;

end

R2 = S2/loop;

figure(1)
plot(Pdb, R(1,:), '-b', Pdb, R(2,:), '-or', Pdb, R(3,:), '-*k', Pdb, R2(1,:), '--b', Pdb, R2(2,:), '--or', Pdb, R2(3,:), '--*k')
grid on
ylabel('Rate (bits/channel use)')
xlabel('P (dBm)')
legend('EH_{DC} = 0.2 (mW) - AC computing','EH_{DC} = 0.4 (mW) - AC computing','EH_{DC} = 0.6 (mW) - AC computing', 'EH_{DC} = 0.2 (mW) - DC computing','EH_{DC} = 0.4 (mW) - DC computing','EH_{DC} = 0.6 (mW) - DC computing')




