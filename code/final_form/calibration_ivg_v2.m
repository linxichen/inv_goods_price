% utility function is linear in this version. C=1
% monthly calibration
function [Kmax_lower,Kmin_upper,Kmin_lower]=calibration_ivg_v2(X,PX_low,q_grid)
% input:
% X: prod grid
% PX_low: transition matrix of X
% q_grid: investment goods price grid

% output:
% K_max_lower: lower bound of k_grid(end)
% K_min_upper: upper bound of k_grid(1)
% K_min_lower: lower bound of k_grid(1)

aalpha = .3; % capital share
ddelta = .02;
%C_ss = Y_ss-I_ss; % agg consumption
C_ss = 1; % set it to 1

bbeta = 0.996; % discount rate
nu = 0.6; % labor share
xi_n = 1; % labor disutility

nq = length(q_grid);
K_min1 = zeros(nq,1);K_max2 = zeros(nq,1);K_max3 = zeros(nq,1);K_min4 = zeros(nq,1);
for i = 1:nq
    p_ss = q_grid(i);
    % condition 1
    ggamma = 1/(1-aalpha-nu);
    a1 = (bbeta^(1-nu)*C_ss^(-nu))^ggamma*(nu/xi_n)^(nu*ggamma)*((1-(1-ddelta)*bbeta)*p_ss)^((nu-1)*ggamma);
    b1 = PX_low(end,:)*(X.^ggamma);
    K_min1(i) = a1*b1/(1-ddelta); % k_grid(1)<K_min1
    % condition 2
    a2 = (bbeta^(1-nu)*C_ss^(-nu))^ggamma*(nu/xi_n)^(nu*ggamma)*((1-(1-ddelta)*bbeta)*p_ss)^((nu-1)*ggamma);
    b2 = PX_low(1,:)*(X.^ggamma);
    K_max2(i) = a2*b2/(1-ddelta); % k_grid(end)>K_max2
    % condition 3
    a3 = (bbeta^(1-nu)*C_ss^(-nu))^ggamma*(nu/xi_n)^(nu*ggamma)*((1-(1-ddelta)*bbeta)*p_ss)^((nu-1)*ggamma);
    b3 = PX_low(end,:)*(X.^ggamma);
    K_max3(i) = a3*b3; % k_grid(end)>K_max3
    % condition 4
    a4 = (bbeta^(1-nu)*C_ss^(-nu))^ggamma*(nu/xi_n)^(nu*ggamma)*((1-(1-ddelta)*bbeta)*p_ss)^((nu-1)*ggamma);
    b4 = PX_low(1,:)*(X.^ggamma);
    K_min4(i) = a4*b4; % k_grid(1)>K_min4    
end
Kmin_upper = min(K_min1);
Kmin_lower = max(K_min4);
Kmax_lower = max(max(K_max2),max(K_max3));


