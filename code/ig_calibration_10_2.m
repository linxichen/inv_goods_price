% calibration of a RBC model, Yang Yu 10/2/2014
% assume there's no agg and idiocy shock

function ss = ig_calibration_10_2(x) % x is N+1 vector -- consumption and N grid of capital stock 
% para
aalpha = .7; % capital share
ddelta = .02; % depreciate
bbeta = .98; % discount
L0 = 1; % normalize the total labor supply to 1
N = 2; % number of grid
A = (.1:.1:.1*N)'; % idiocy productivity

% load variables
C = x(1,1);
K = x(2:N+1,1);
W = C/L0; % get wage from household budget constraint
L = ((1-aalpha).*A/W).^(1/aalpha).*K;% get labor input from firm's FOC

% N+1 equations that can solve for x
SS = zeros(N+1,1); % N+1 equations
eq = 0;
ii = 1;
while ii<N+1
    eq = eq+1;
    a = A(ii);% idiocyncratic prod. of firm ii
    k = K(ii);% capital stock of firm ii
    l = L(ii); % labor input of firm ii
    SS(eq,1) = a*aalpha*k^(aalpha-1)*l^(1-aalpha);% MPK
    SS(eq,2) = 1/bbeta-(1-ddelta);
    ii = ii+1;
end
% labor market clear
eq = eq+1;
SS(eq,1) = sum(L);
SS(eq,2) = L0;

% economy resource constraint
eq = eq+1;
SS(eq,1) = sum(A.*K.^aalpha.*L.^(1-aalpha));
SS(eq,2) = C+ddelta*sum(K);

ss = sum(abs(SS(:,1)-SS(:,2)));
    