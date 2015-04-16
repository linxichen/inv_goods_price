%% do predictability test with AR(1)
%% should do this with annual data later.
%% detrend raw data with linear trend, since HP filter contaminate the data with future info.
X0 = xlsread('ivg raw data.xls','B2:L164');
run tax_series
cpi = X0(:,2);
ivg = X0(:,3);
unc_sale = X0(:,1);
ivg_real = ivg./cpi;
T = length(ivg_real);
% generate seasonal dummy
dum = zeros(T,4);
for t = 1:T
    for ii = 0:3
        if mod(t,4) == ii
            dum(t,ii+1) = 1;
        end
    end
end

% filter the data with a polynomial of time trend, seasonal dummy and AR(1)
% investment goods price
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';ivg_real(1:T-1)']';
yy = ivg_real(2:T);
beta_ivg = (xx'*xx)\(xx'*yy);
eps = yy-xx*beta_ivg;
% uncertainty
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';unc_sale(1:T-1)']';
yy = unc_sale(2:T);
beta_unc = (xx'*xx)\(xx'*yy);
eps_unc = yy-xx*beta_unc;
% capital tax
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';tau_K(1:T-1)']';
yy = tau_K(2:T);
beta_tau_K = (xx'*xx)\(xx'*yy);
eps_tau_K = yy-xx*beta_tau_K;
% FFR
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';FFR(1:T-1)']';
yy = FFR(2:T);
beta_FFR = (xx'*xx)\(xx'*yy);
eps_FFR = yy-xx*beta_FFR;
% MN
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';MN(1:T-1)']';
yy = MN(2:T);
beta_MN = (xx'*xx)\(xx'*yy);
eps_MN = yy-xx*beta_MN;
% OP
xx = [ones(1,T-1);1:T-1;(1:T-1).^2;dum(1:T-1,1:3)';log(OP(1:T-1))']';
yy = OP(2:T);
beta_OP = (xx'*xx)\(xx'*yy);
eps_OP = yy-xx*beta_OP;
%% use uncertainty to predict eps
N = 10;% number of lag 
beta0 = zeros(N,2);
rsq0 = zeros(N,1);
se0 = zeros(N,2);
for j = 1:N
xx = [ones(T-1-j,1),eps_unc(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta0(j,:) = beta_temp;
rsq0(j) = rsq_temp;
se0(j,:) = se_temp;
end

%% use capital tax to predict ivg
beta1 = zeros(N,2);
rsq1 = zeros(N,1);
se1 = zeros(N,2);
for j = 1:N
xx = [ones(T-1-j,1),eps_tau_K(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta1(j,:) = beta_temp;
rsq1(j) = rsq_temp;
se1(j,:) = se_temp;
end

%% use federal fund rate to predict ivg
beta2 = zeros(N,2);
rsq2 = zeros(N,1);
se2 = zeros(N,2);
for j = 1:N
xx = [ones(T-1-j,1),eps_FFR(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta2(j,:) = beta_temp;
rsq2(j) = rsq_temp;
se2(j,:) = se_temp;
end

%% use news of military spending (Ramey Shapiro) to predict ivg
beta3 = zeros(N,2);
rsq3 = zeros(N,1);
se3 = zeros(N,2);
for j = 1:N
xx = [ones(T-1-j,1),eps_MN(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta3(j,:) = beta_temp;
rsq3(j) = rsq_temp;
se3(j,:) = se_temp;
end

%% use real oil cost (Ramey Vine) to predict ivg
beta4 = zeros(N,2);
rsq4 = zeros(N,1);
se4 = zeros(N,2);
for j = 1:N
xx = [ones(T-1-j,1),eps_OP(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta4(j,:) = beta_temp;
rsq4(j) = rsq_temp;
se4(j,:) = se_temp;
end

%% use them all
beta5 = zeros(N,6);
rsq5 = zeros(N,1);
se5 = zeros(N,6);
for j = 1:N
xx = [ones(T-1-j,1),eps_unc(1:end-j),eps_tau_K(1:end-j),eps_FFR(1:end-j),eps_MN(1:end-j),eps_OP(1:end-j)];
yy = eps(1+j:end,1);
beta_temp = (xx'*xx)\(xx'*yy);
rsq_temp = var(xx*beta_temp)/var(yy);
eps_temp = yy-xx*beta_temp;
% Newey West SE, assuming L = 2
se_temp = NeweyWest(eps_temp,xx,2,1);
beta5(j,:) = beta_temp;
rsq5(j) = rsq_temp;
se5(j,:) = se_temp;
end