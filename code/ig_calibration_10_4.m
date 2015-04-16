%% Housekeeping
clear all
close all
clc

%% Parameters
aalpha = 0.3; % capital share
ggamma = 0.6; % labor share, DRTS here
ddelta = .02; % depreciate
bbeta = .996; % discount

%% Initialize 
N = 11; % number of firms
rrhox = 0.98;
ssigmax = 0.03;
% minA = 0.9; maxA = 1.1;
[A,P] = tauchen(N,0,rrhox,ssigmax,3);
A = exp(A');
% A = linspace(minA,maxA,N);
%A = normrnd(1,0.3,1,N);

%% Find division of labor
sol = fsolve(@(n) dividelabor(n,aalpha,bbeta,ggamma,ddelta,A),[ones(N,1)./N ;1]);
nss = real(sol(1:N)); wss = real(sol(end));

%% Backout division of capital and aggregate capital
kss = ((1/bbeta-1+ddelta)/aalpha)^(1/(aalpha-1)).*A'.^(1/(1-aalpha)).*nss.^(ggamma/(1-aalpha));
Kss = sum(kss);

%% Backout the rest of variables
yss = A'.*kss.^(aalpha).*nss.^(ggamma);
Yss = sum(yss);
Css = Yss - ddelta*Kss;
profit = yss - wss.*nss - ddelta*kss;
Yss