%% investment goods model 
% this version has a quadratic adjustment cost
% make the grid even
% introduce loss of disinvestment (lower price of resale)
clear all
clc
close

%% "Deep" Parameters
bbeta = 0.996;
ttau = 0.00; % search cost
aalpha = 0.3; % y = a*k^aalpha*l^v
v = 0.6;
aalpha0 = 1.25; % search elasticity
ddelta = 0.02;
pphi = 0.0000001; % price of disinvestment relative to investment
rrhox = 0.98;
ssigmax = 0.08;
ppsi = 0; % quadratic cost of investment adjustment

%% Accuracy Controls
nk = 101; % number of grid points on capital stock
nx = 11; % number of grid points on idiosyncractic prod.
m = 3; % support is m s.d. away from mean
tol = 1e-3;

%% Initialize
W = ones(nk,nx); % value of matching with investment goods producer after paying the search cost
U = ones(nk,nx); % value of not going to search 
V_old = ones(nk,nx); % maximized value 
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
K = linspace(0.7*10,1.3*10,nk)'; % calibrated from ig_calibration_10_4
w = 3.63; % wage, taken as given, calibrated from ig_calibration_10_4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the capital stock of next period if the firm does not do investment. 
K_tmr = K*(1-ddelta);
Ids = K;
for i = 1:length(K)
[~, ids] = min(abs(K_tmr(i)-K));% ids is the index of 
Ids(i) = ids; 
end
KK = ones(nk,nx); % index of tomorrow's capital stock, potentially not sliding down
Y = ones(nk,nx);
L = ones(nk,nx); % labor decision, given wage and productivity
d = ones(nk,nx);% this is the gross profit, i.e. investment cost is not substracted
% for each capital, compute the labor, production and gross profit
for ik = 1:nk
    for ix = 1:nx
        k = K(ik);
        x = X(ix);
        L(ik,ix) = (w/(x*v))^(1/(v-1))*k^(aalpha/(1-v)); % From FOC wrt l
        Y(ik,ix) = x*k^aalpha*L(ik,ix)^v;
        d(ik,ix) = Y(ik,ix)-w*L(ik,ix);
    end
end

PX_unc = PX^10000;
ppi = PX_unc(1,:);

%% Iteration Step
%qq = linspace(0.1,5,100);
qq = 2;
kdstar = cell(length(qq),1);% each cell is a stationary dist for one q
ii0 = zeros(length(qq),1); % total investment 
for i_q = 1:length(qq)
q = qq(i_q);
iter = 0;
err = 1;
while err > tol
    EV = V_old*PX'; % E(V(k_{t+1},x_{t+1}))
    U = d+bbeta*EV(Ids,:); % capital stock slide down
    mu = 1/(1+theta_old^(-aalpha0))^(1/aalpha0); % matching probability
    for ik = 1:nk
        Ind = K-(1-ddelta)*K(ik,1); % index of investment, allow for disinvest in this version
        AC = ppsi*(Ind/K(ik,1)).^2; % quadratic adjustment cost
        %% check. 
        W_temp = repmat(d(ik,:),nk,1)-ttau-q*mu*repmat(Ind.*(Ind>=0),1,nx)...
            -pphi*q*mu*repmat(Ind.*(Ind<0),1,nx)-repmat(AC,1,nx)+bbeta*mu*EV+bbeta*(1-mu)*repmat(EV(Ids(ik),:),nk,1); % d-q*i-qac+E(V(k_{t+1}))
        [cc,ii] = max(W_temp); % for each k and x of today, find capital of tomorrow that maximize the value
        W(ik,:) = cc'; % maximized value conditional on being matched.
        I(ik,:) = Ind(ii'); % This is the intended investment level!!! Not the actually realized level, it could be in Ids
        KK(ik,:) = ii'; % destination of the capital. K_temp(ik,ix) indicates the capital stock of tomorrow of this catergory
    end
    EV1 = W; % value of going to search
    EV2 = U; % value of not going search
    V_new = EV1.*(EV1>=EV2)+EV2.*(EV2>EV1);
    err = norm(V_old-V_new,Inf)
    V_old = V_new;

    KK_final = KK.*(EV1>=EV2)+repmat(KK(Ids),1,nx).*(EV1<EV2); % intended investment
    % use the policy function to get the stationary distribution of capital
    
    % Find the new distribution by definition, slow but accurate version
    Kd_new = zeros(nk,nx);
    for i_kplus = 1:nk
        for i_xplus = 1:nx
            for i_x = 1:nx
                list_success_invest = (W(:,i_x)>U(:,i_x)).*(KK(:,i_x)==i_kplus);
                list_failed_invest  = (W(:,i_x)>U(:,i_x)).*(Ids==i_kplus);
                list_intend_deprec  = (W(:,i_x)<U(:,i_x)).*(Ids==i_kplus);
                mass_success_invest = mu*sum(Kd_old(:,i_x).*list_success_invest);
                mass_failed_invest  = (1-mu)*sum(Kd_old(:,i_x).*list_failed_invest);
                mass_intend_deprec  = sum(Kd_old(:,i_x).*list_intend_deprec);
                Kd_new(i_kplus,i_xplus) = Kd_new(i_kplus,i_xplus) ...
                                        + PX(i_x,i_xplus).*(...
                                            mass_success_invest+mass_failed_invest+mass_intend_deprec ...
                                        );
            end
        end
    end
    dist_err = norm(Kd_old-Kd_new,Inf);
    Kd_old = Kd_new;
    ttheta_old = sum(sum(Kd_old.*(W>U)));
end
%% Now the individual value function is found, given q forever
A = (W>U);
sum(sum(A))/(size(A,1)*size(A,2));
change = (KK-(1-ddelta)*repmat(K,1,nx)).*(W>U);
ii0(i_q) = sum(sum(change));
kdstar{i_q} = Kd_old;
end;



    
    
    
