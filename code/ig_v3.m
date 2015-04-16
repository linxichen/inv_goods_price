%% investment goods model 
% this version has a quadratic adjustment cost
% make the grid even
% introduce loss of disinvestment (lower price of resale)
clear all
clc
close

%% "Deep" Parameters
q = 2; % price of investment goods
bbeta = 0.996;
ttau = 0.3; % search cost
aalpha = 0.25; % y = a*k^aalpha*l^v
v = 0.5;
aalpha0 = 1.25; % search elasticity
ddelta = 0.02;
pphi = .7; % price of disinvestment relative to investment
rrhox = 0.98;
ssigmax = 0.03;
ppsi = 0; % quadratic cost of investment adjustment
%% Accuracy Controls
nk = 10; % number of grid points on capital stock
nx = 10; % number of grid points on idiosyncractic prod.
m = 3; % support is m s.d. away from mean
tol = 1e-3;

%% Initialize
W = ones(nk,nx);% value of matching with investment goods producer after paying the search cost
U = ones(nk,nx);% value of not going to search 
V_old = ones(nk,nx);% maximized value 
V_new = ones(nk,nx);
Kd_old = ones(nk,nx)/(nk*nx); % distribution of capital, the elements should add up to one
Kd_new = ones(nk,nx)/(nk*nx);
theta_old = 1; % measure of firms who decide to search
theta_new = 1; % measure of firms who decide to search
I = ones(nk,nx); % policy function
P = zeros(nk*nx,nk*nx); % This is the "joint" transition matrix, useful later

%% Exogenuous Things
[X,PX] = tauchen(nx,0,rrhox,ssigmax,m);
X = exp(X); % Convert back to level
X0 = X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = []; % calibrate this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the capital stock of next period if the firm does not do investment. 
K_tmr = K*(1-ddelta);
Ids = K;
for i = 1:length(K)
[~, ids] = min(abs(K_tmr(i)-K));% ids is the index of 
Ids(i) = ids; 
end
KK = ones(nk,nx); % index of tomorrow's capital stock
Y = ones(nk,nx);
L = ones(nk,nx);
d = ones(nk,nx);% this is the gross profit, i.e. investment cost is not substracted
% for each capital, compute the labor, production and gross profit
for ik = 1:nk
    for ix = 1:nx
        k = K(ik);
        x = X(ix);
        L(ik,ix) = (w/(x*v))^(1/(v-1))*k^(aalpha/(1-v));
        Y(ik,ix) = x*k^aalpha*L(ik,ix)^v;
        d(ik,ix) = Y(ik,ix)-w*L(ik,ix);
    end
end

PX_unc = PX^10000;
ppi = PX_unc(1,:);
%% Iteration Step
iter = 0;
err = 1;
while err > tol
    EV = V_old*PX'; % E(V(k_{t+1},x_{t+1}))
    U = d+bbeta*EV(Ids,:); % capital stock slide down                                
    for ik = 1:nk
        Ind = K(:,1)-(1-ddelta)*K(ik,1); % index of investment, allow for disinvest in this version
        AC = ppsi*(Ind/K(ik,1)).^2; % quadratic adjustment cost
        %% check. 
        W_temp = repmat(d(ik,:),nk,1)-q*repmat(Ind.*(Ind>=0),1,nx)...
            -pphi*q*repmat(Ind.*(Ind<0),1,nx)-repmat(AC,1,nx)+bbeta*EV; % d-q*i-qac+E(V(k_{t+1}))
        [cc,ii] = max(W_temp); % for each k and x of today, find capital of tomorrow that maximize the value
        W(ik,:) = cc'; % maximized value conditional on being matched.
        I(ik,:) = Ind(ii');
        KK(ik,:) = ii'; % destination of the capital. K_temp(ik,ix) indicates the capital stock of tomorrow of this catergory
    end
    mu = 1/(1+theta_old^(-aalpha0))^(1/aalpha0); % matching probability
    EV1 = -ttau+((1-mu)*U+mu*W); % value of going to search
    EV2 = U; % value of not going search
    V_new = EV1.*(EV1>=EV2)+EV2.*(EV2>EV1);
    %% match rate?
    Kd_temp = Kd_old.*(EV1>=EV2); % distribution that put zero on the innactive region, 
    theta_new = sum(sum(Kd_temp)); % measure of activity region
    KK_final = KK.*(EV1>=EV2)+repmat(KK(Ids),1,nx).*(EV1<EV2); % intended investment
    % use the policy function to get the stationary distribution of capital
    % define a large transition matrix
    %% check: depreciation
    for i = 1:nk*nx
        % i is the "composite" index for (i_k,i_x), now we convert it back to
        % separate tuple index (i_k,i_x) for current shock
        [i_k,i_x] = ind2sub([nk,nx],i);
        for j = 1:nk*nx
        % j is the composite index for tomorrow's shock
        [j_k,j_x] = ind2sub([nk,nx],j);
        % the transition prob. from [ik,ix] to [KK(ik,ix),ix2] is
        % PX(ix,ix2). Otherwise the prob. is zero
        %% check
        P(i,j) = mu*(j_k == KK_final(i_k,i_x))*PX(i_x,j_x)+(1-mu)*(j_k == K(Ids))*PX(i_x,j_x);
        % P(i,j) = (mu*(j_k == KK_final(i_k,i_x))*PX(i_x,j_x)+(1-mu)*(j_k == K(Ids))*PX(i_x,j_x))*(indicator of invesment)...
          %  +(j_k == K(Ids))*PX(i_x,j_x)*(indicator of not investment);
        end
    end
    %% check for this...
    P_unc = P^10000;
    Kd_new = reshape(P_unc(1,:),nk,nx);
    iter = iter + 1;
    err = max(max(abs(V_new(:)-V_old(:))))+max(max(abs(Kd_new-Kd_old)))+abs(theta_old-theta_new)
    V_old = V_new;
    Kd_old = Kd_new;
    theta_old = theta_new;
end

    
    
    
