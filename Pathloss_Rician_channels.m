function [Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb)
%n_u - the number of user
%n_t - the number of transmit antennas
%pathloss
pl = 2.6;

K = 10.^(Kdb/10); %Rician factor

%Rician channels
mu = sqrt( K/((K+1))); 

s = sqrt( 1/(2*(K+1)));

Hw = mu + s*(randn(n_t,n_u) + 1j*randn(n_t,n_u));

%with pathloss
Hpl = sqrt(1./(dist^pl))*Hw;
end
