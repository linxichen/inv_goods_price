% this file compute comparative static: compare two permanent regimes (high and low uncertainty). 
% In each regime, step 1: every period fixing price to p0, get value functions
%                 step 2: change price only for the first period, compute
%                 value functions with U and EV got in step 1.
%                 step 3: find the optimal price and corresponding
%                 investment. 
% 2015/2/13

aalpha0 = 1.25; % search efficiency
aalpha1 = .8; % elascticity
bbeta = 0.996;
aalpha = 0.3; % y = a*k^aalpha*l^v
ddelta = 0.02;
rrhox = 0.98;
Ssigmax = [0.04,0.06]; % two regimes: boom and recession
nsx = length(Ssigmax);
nq = 100;
q_star = zeros(nsx,1);% optimal price
inv_star = zeros(nsx,1); % final investment 
resale = 0.2; % resale price 
ttau = .5; % fixed cost

WW = cell(nq,nsx); % value function W
II = cell(nq,nsx); % active zone (extensive)
DD = cell(nq,nsx); % matrix of intensive demand 
ttheta_old = ones(nq,nsx); % tightness ratio
mmu_old = .8*ones(nq,nsx); % finding rate
inv = ones(nq,nsx); % investment
rvn = ones(nq,nsx); % revenues

%% Accuracy control
nk = 120; % number of grid points on capital stock
nx = 20; % number of grid points on idiosyncractic prod.
nK = 11;
m = 3; % support is m s.d. away from mean
tol = 1e-6;
tolmmu = 1e-2;


for isigma = 1:nsx
    ssigmax = Ssigmax(isigma);
    %% Grids
    nq = 100; % grid of price
    q0 = .5; % price in the long run
    q_grid = q0+linspace(-0.1,0.1,nq)';
    [X,PX] = tauchen(nx,0,rrhox,ssigmax,m);
    X = exp(X);% Convert back to level
    uncond_X = PX^3000;
    % Capital stuff
    max_k = 100;
    min_k = 50;
    k_grid = linspace(min_k,max_k,nk)'; % calibrated from ig_calibration_10_4
    noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
    for i_k = 1:nk
        [~,noinvest_ind(i_k)] = min(abs(k_grid-(1-ddelta)*k_grid(i_k)));
    end
    w = 1; % Normalized to one
    inv_mat = repmat(k_grid',nk,1)-(1-ddelta)*repmat(k_grid,1,nk);
    pos_inv = inv_mat>0;
    neg_inv = inv_mat<=0;
    %% Initialize value functions
    W_old = ones(nk,nx); % value of matching with investment goods producer after paying the search cost
    W_new = W_old;
    U_old = ones(nk,nx); % value of not going to search, not related to current q
    U_new = U_old;
    V_old = ones(nk,nx); % maximized value after discrete choice
    V_new = ones(nk,nx); %
    max_movingpart = zeros(nk,nx); % k,kplus
    koptind_active = zeros(nk,nx);
    
    profit = zeros(nk,nx); % Compute operation profit resulting from labor choice
    L = zeros(nk,nx);
    nu = 0.5; % y = x*k^aalpha*L^nu
    for i_k = 1:nk
        for i_x = 1:nx
            L(i_k,i_x) = (X(i_x)*nu*k_grid(i_k)^aalpha/w)^(1/(1-nu));
            profit(i_k,i_x) = X(i_x)*k_grid(i_k)^aalpha*L(i_k,i_x)^nu-L(i_k,i_x)*w;
        end
    end
    diff = 10;
    outer_iter = 0;
    % Start VFI here
    err = 10;
    iter = 0;
    while err > tol
        % ttheta = ttheta_old(i_K,i_ssigmak);
        % mmu = mmu_old(i_K,i_ssigmak);
        mmu = .9;
        EV = V_old(:,:)*PX';
        U_new = profit + bbeta*EV(noinvest_ind,:);
        for i_x = 1:nx
            [max_movingpart(:,i_x),koptind_active(:,i_x)] = max(bbeta*repmat(EV(:,i_x)',nk,1) - q0*inv_mat.*pos_inv - resale*q0*inv_mat.*neg_inv,[],2);
        end
        W_new = profit + mmu*max_movingpart - ttau + bbeta*(1-mmu)*EV(noinvest_ind,:);
        
        %====================================================%
        V_new = max(W_new,U_new);
        
        err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
        V_old = V_new;
        W_old = W_new;
        U_old = U_new;
        iter = iter + 1;
        iter
    end
    EV = V_old(:,:)*PX';
    U = U_old;
    for i = 1:nq
        eer = 1;
        q = q_grid(i);
        mmu_new = mmu_old(i);
        D_new = zeros(nk,nx);
        while eer>tolmmu
            mmu = mmu_new;
            for i_x = 1:nx
                [max_movingpart(:,i_x),koptind_active(:,i_x)] = max(bbeta*repmat(EV(:,i_x)',nk,1) - (q*inv_mat.*pos_inv + resale*q*inv_mat.*neg_inv),[],2);
                D_new(:,i_x) = k_grid(koptind_active(:,i_x))-k_grid(noinvest_ind);
            end
            W_new = profit + mmu*max_movingpart - ttau + bbeta*(1-mmu)*EV(noinvest_ind,:);
            I_new = sum(sum((W_new>U)))/(nk*nx); % fraction being active
            ttheta = 1/I_new;
            mmu_new = 1./(1+ttheta^(-aalpha0))^(1/aalpha0);
            eer = abs(mmu-mmu_new);
        end
        WW{i,isigma} = W_new;
        II{i,isigma} = (W_new>U);
        mmu_old(i,isigma) = mmu_new;
        DD{i,isigma} = D_new;
        inv(i,isigma) = mmu_new*sum(sum(D_new.*(W_new>U)));
        inv_temp = mmu_new*sum(sum(D_new.*((D_new>0)+(D_new<0)*resale).*(W_new>U)));
        rvn(i,isigma) = q*inv(i,isigma);
    end
    [~,optq] = max(rvn(:,isigma)); % best price
    q_star(isigma) = q_grid(optq);
    inv_star(isigma) = inv(optq,isigma); % final investment
end
plot(inv(1:55,1),q_grid(1:55),'LineWidth',2.5);
xlabel('Investment','FontSize', 10) % x-axis label
ylabel('Price','FontSize', 10) % y-axis label
title('Comparing boom and recession', 'FontSize', 12);
hold on
plot(inv_star(1),q_star(1),'kd','LineWidth',4.5);
plot(inv(1:55,2),q_grid(1:55),'r','LineWidth',2.5);
plot(inv_star(2),q_star(2),'yd','LineWidth',4.5);
legend('Demand curve in boom','Equilibrium in boom','Demand curve in recession','Equilibrium in recession')