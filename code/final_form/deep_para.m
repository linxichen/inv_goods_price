bbeta = 0.996; % if you don't know beta... good luck
ttau = 0.1; % search cost
aalpha = 0.3; % y = a*k^aalpha*l^v
v = 0.6; % labor share
aalpha0 = 2; % search elasticity
ddelta = 0.02;
pphi = 0.000000; % price of disinvestment relative to investment
rrhox = 0.95; % persistence of idio TFP
ppsi = 0.00; % quadratic cost of investment adjustment
rrhoz = 0.95; % persistence of agg TFP
ssigmaz = 0.01; % std of z innov
ssigmax_low = 0.01; % low std of x innov
ssigmax_high= 0.03; % high std of x innov
Pssigmax = [0.9 0.1; 0.1 0.9]; % Transition prob of ssigmax

% Notational difference
nu = v;
ppsi_n = 1;
xi_n = ppsi_n;
C_ss = 1;

