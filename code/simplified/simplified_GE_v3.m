%% Housekeeping
% new variable: consumption C.
% In this version, producing investment goods does not incur any cost, and
% there's no search cost, we always have C = Y
% log(C) = pphi_CC+pphi_CK*log(K)+pphi_Cssigmax*log(ssigmax)

% 4/8/2015 new version
% investment has cost psi, thus consumption needs to be solved numerically
% 4-11
% deduct cost in firm's revenue function in the simulation part
clear; close all; clc;
bbeta = 0.996;
ttau = 0.1; % search cost
aalpha = 0.3; % y = a*k^aalpha*l^v
v = 0.6;
aalpha0 = 2; % search elasticity
ddelta = 0.02;
pphi = 0.0000001; % price of disinvestment relative to investment
rrhox = 0.98;
ppsi = 0.001; % quadratic cost of investment adjustment

rrhoz = 0.95;
ssigmaz = 0.01;
ssigmax_low = 0.03;
ssigmax_high= 0.06;
Pssigmax = [0.9 0.1; 0.1 0.9];

%% Accuracy control
nk = 30; % number of grid points on capital stock
nfine = 100;
nx = 5; % number of grid points on idiosyncractic prod.
nz = 5; % number of grid points on aggregate productivity
ns = nx*nz*2;
nK = 7;
nq =20;
m = 3; % support is m s.d. away from mean
tol = 1e-2;
maxiter = 2000;
damp = 0.3;

%% Grids
% ssigmagrid = [0.8;1.2];
% Pssigma = 0.5*ones(2,2);
[Z,PZ] = tauchen(nz,0,rrhoz,ssigmaz,m);
Z = exp(Z);
ssigmax = [ssigmax_low,ssigmax_high];
[X,PX_low] = tauchen(nx,0,rrhox,ssigmax_low,m);
PX_high = tauchen_givengrid(0,rrhox,ssigmax_high,X);
X = exp(X);% Convert back to level
P = zeros(2*nx,2*nx);
low = 1;
high = 2;
uncond_X_low = PX_low^3000;
uncond_X_high = PX_high^3000;
for i = 1:ns
    [i_x,i_z,i_ssigmax] = ind2sub([nx nz 2],i);
    for j = 1:ns
        [j_x,j_z,j_ssigmax] = ind2sub([nx nz 2],j);
        P(i,j) = PZ(i_z,j_z)*Pssigmax(i_ssigmax,j_ssigmax)*( (j_ssigmax==low)*PX_low(i_x,j_x)+ (j_ssigmax==high)*PX_high(i_x,j_x) );
    end
end

% Capital stuff
min_k = 10;
k_grid = min_k*ones(nk,1); % calibrated from ig_calibration_10_4
for i_k = 1:nk-1
    k_grid(i_k+1) = (1.02)*k_grid(i_k);
end
fine_grid = linspace(k_grid(1),k_grid(nk),nfine)';
noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
invest_ind = ones(nk,1);
for i_k = 1:nk
    [~,noinvest_ind(i_k)] = min(abs(k_grid-(1-ddelta)*k_grid(i_k)));
    [~,invest_ind(i_k)] = min(abs(k_grid-(1.02)*k_grid(i_k)));
end

noinvest_ind_fine = ones(nfine,1); % for each k, the index of tmr k if no invest
invest_ind_fine = ones(nfine,1);
for i_k = 1:nfine
    [~,noinvest_ind_fine(i_k)] = min(abs(fine_grid-(1-ddelta)*fine_grid(i_k)));
    [~,invest_ind_fine(i_k)] = min(abs(fine_grid-(1.02)*fine_grid(i_k)));
end
w = 1; % Normalized to one
inv_mat = repmat(k_grid',nk,1)-(1-ddelta)*repmat(k_grid,1,nk);
pos_inv = inv_mat>0;
neg_inv = inv_mat<=0;

inv_mat_fine = repmat(fine_grid',nfine,1)-(1-ddelta)*repmat(fine_grid,1,nfine);
pos_inv_fine = inv_mat_fine>0;
neg_inv_fine = inv_mat_fine<=0;

K_grid = linspace(k_grid(1),k_grid(nk),nK)'; % Aggregate capital grid

% Variance grid
qmin = .1;
qmax = .9;
psi = qmin; % cost of investment goods
% a1 = (log(qmax)-log(qmin))/pi;
% b1 = (log(qmax)+log(qmin))/2;
q_grid = linspace(qmin,qmax,nq); % grid for current q
% q_tilde = (tan(log(q_grid)-b1))/a1; % log(q) = atan(q_tilde)*a1+b1; make sure log(q) is always between log(qmin) and log(qmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure q is always between .2 and .6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule of q as function of aggregate states
pphi_qC = log(1);
pphi_qK = 0.0; % aggregate prediction rule for q
pphi_qz = 0.0;
pphi_qssigmax = 0.00;

pphi_KK = 0.99; 
pphi_KC = 0.1; 
pphi_Kz = 0.01;
pphi_Kssigmax = 0.01;% Aggregate Law of motion for aggregate capital

pphi_CK = 0.99; 
pphi_CC = 0.1; 
pphi_Cz = 0.01;
pphi_Cssigmax = 0.01;% Aggregate Law of motion for consumption
%% Initialize value functions
W_old = ones(nk,ns,nK,nq); % value of matching with investment goods producer after paying the search cost
W_new = W_old;
U_old = ones(nk,ns,nK); % value of not going to search, not related to current q
U_new = U_old;
V_old = ones(nk,ns,nK,nq); % maximized value after discrete choice
V_new = ones(nk,ns,nK,nq); %
V_new_fine = zeros(nfine,ns,nK,nq);
% profit_fine = zeros(nfine,ns,nq);
EV_new_fine = zeros(nfine,ns,nK,nq);
W_new_fine = zeros(nfine,ns,nK,nq);
U_new_fine = zeros(nfine,ns,nK);
koptind_active = zeros(nk,ns,nK,nq);
active = zeros(nk,ns,nK,nq);

% Compute operation profit 
profit = zeros(nk,ns); 
L = zeros(nk,ns);
psi_n = .2;
nu = .2;
for i_k = 1:nk
    for i_s = 1:ns
        [i_z,i_x,~] = ind2sub([nz nx 2],i_s);
        profit(i_k,i_s) = (1-nu)*(nu/psi_n)^(nu/(1-nu))*(Z(i_z)*X(i_x)*k_grid(i_k)^aalpha)^(1/(1-nu));
    end
end

% Prepare for Simulation stuff
T = 100;
burnin = 3;
kss = mean(k_grid);
Ksim = kss*ones(1,T);
Csim = Ksim; % production and consumption of final goods
qsim = mean(q_grid)*ones(1,T);
x_cdf_low = cumsum(PX_low,2);
x_cdf_high = cumsum(PX_high,2);
revenue = zeros(nq,1);
demand = zeros(nq,T);
dist_k = ones(nfine,nx,T)/(nfine*nx);
% for i_s = 1:nx
%     dist_k(floor(nfine/2),i_s,1) = uncond_X_low(1,i_s); % Economy starts at low volatility state
% end

% New Discrete Simulation
rng('default');
rng(2015);
markov_shock = rand(1,T);
ssigmaxsim = ones(1,T);
ssigmax_cdf = cumsum(Pssigmax,2);
for t = 1:T-1
    ssigmaxsim(t+1) = find(ssigmax_cdf(ssigmaxsim(t),:) >= markov_shock(t),1,'first');
end

zindsim = ceil(nz/2)*ones(1,T);
zsim = ones(1,T);
z_cdf = cumsum(PZ,2);
z_index_shock = rand(1,T);
for t = 1:T-1
    zindsim(t+1) = find(z_cdf(zindsim(t),:) >= z_index_shock(t),1,'first');
    zsim(t+1) = Z(zindsim(t+1));
end

tic
%% Main Body of KS iter
diff = 10;
outer_iter = 0;

while ((diff > tol) && (outer_iter < maxiter))
    % Inner loop
    err = 10;
    iter = 0;
    while err > tol
        for i_s = 1:ns
            for i_K = 1:nK
                % What is agg state today
                [i_z,i_x,i_ssigmax] = ind2sub([nz nx 2],i_s);                
                [~,i_Kplus] = min(abs(K_grid-exp(pphi_KC+pphi_KK*log(K_grid(i_K))+pphi_Kssigmax*log(i_ssigmax)+pphi_Kz*log(Z(i_z)) )));
                [~,i_qplus] = min(abs(q_grid-exp(pphi_qC+pphi_qK*log(K_grid(i_K))+pphi_qssigmax*log(i_ssigmax)+pphi_qz*log(Z(i_z)) )));
                C = exp(pphi_CC+pphi_CK*log(K_grid(i_K))+pphi_Cssigmax*log(i_ssigmax)+pphi_Cz*log(Z(i_z)));
                                
                % Find expected value across i_splus
                EV_invest = V_old(invest_ind,:,i_Kplus,i_qplus)*P(i_s,:)';
                EV_noinvest = V_old(noinvest_ind,:,i_Kplus,i_qplus)*P(i_s,:)';
                
                % U has no choice, so no need to find max in a loop.
                U_new(:,i_s,i_K) = profit(:,i_s).*C.^(1/(nu-1)) + bbeta*EV_noinvest;
                
                % Find W, which also without choice
                for i_q = 1:nq
                    q = q_grid(i_q);
                    W_new(:,i_s,i_K,i_q) = profit(:,i_s).*C.^(1/(nu-1)) - q*(k_grid(invest_ind)-(1-ddelta)*k_grid)./C + bbeta*EV_invest;
                    V_new(:,i_s,i_K,i_q) = max(W_new(:,i_s,i_K,i_q),U_new(:,i_s,i_K));
                end
            end
        end
        
        err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
        V_old = V_new;
        W_old = W_new;
        U_old = U_new;
        iter = iter + 1;
        if mod(iter,200) == 0
            disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,diff,iter,err);
            disp(disp_text);
        end
    end
    
    active = W_new > repmat(U_new,1,1,1,nq);
    koptind = repmat(noinvest_ind,1,ns,nK,nq).*(1-active) + active.*repmat(invest_ind,1,ns,nK,nq);
    kopt = k_grid(koptind);
    
    % Interpolate on finer grid
    for i_q = 1:nq
        for i_K = 1:nK
            EV_temp = V_new(:,:,i_K,i_q)*P';
            for i_s = 1:ns
                U_new_fine(:,i_s,i_K) = interp1(k_grid,U_new(:,i_s,i_K),fine_grid,'linear')';
                W_new_fine(:,i_s,i_K,i_q) = interp1(k_grid,W_new(:,i_s,i_K,i_q),fine_grid,'linear')';
                % V_new_fine(:,i_s,i_K,i_q) = interp1(k_grid,V_new(:,i_s,i_K,i_q),fine_grid,'linear')';
                % profit_fine(:,i_s) = interp1(k_grid,profit(:,i_s),fine_grid,'linear')';
                % EV_new_fine(:,i_s,i_K,i_q) = interp1(k_grid,EV_temp(:,i_s),fine_grid,'linear')';
            end
        end
    end
    active_fine = W_new_fine > repmat(U_new_fine,1,1,1,nq);
    koptind_fine = repmat(noinvest_ind_fine,1,ns,nK,nq).*(1-active_fine) + active_fine.*repmat(invest_ind_fine,1,ns,nK,nq);
    kopt_fine = fine_grid(koptind_fine);
    
    %% Given individual policies, simulate a large panel to update aggregate law of motion
    
    % I forgot to zero out after each outer iteration.
    dist_k(:,:,2:T) = zeros(nfine,nx,T-1);
    for t = 1:T
        % Find aggregate stuff today
        Ksim(t) = sum(vec(dist_k(:,:,t).*repmat(fine_grid,1,nx)));
        [~,i_K] = min(abs(Ksim(t)-K_grid));
        % Identify which uncertainty
        whichxind = (1:nx)+nx*(ssigmaxsim(t)-1);
        if ssigmaxsim(t) == low
            whichprob = PX_low;
        elseif ssigmaxsim(t) == high
            whichprob = PX_high;
        end
       
        
        % According to policy functions, find the optimal q
        for i_q = 1:nq
            tot_revenue_grid = (q_grid(i_q)-psi)*dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichxind,i_K,i_q));
            tot_revenue_grid(tot_revenue_grid<0) = 0;
            demand_grid = dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichxind,i_K,i_q));
            revenue(i_q,t) = sum(tot_revenue_grid(:));
            demand(i_q,t) = sum(demand_grid(:));
        end
        [~,i_qmax] = max(revenue(:,t));
        qsim(t) = q_grid(i_qmax);
        
        % find consumption
        % Csim(t) = (nu/psi_n)^nu*sum(vec(dist_k(:,:,t).*(fine_grid.^aalpha*X'*zsim(t))))^(1-nu);
        aa1 = demand(i_qmax,t);
        aa2 = (nu/psi_n)^(nu/(1-nu))*sum(vec(dist_k(:,:,t).*(fine_grid.^aalpha*X'*zsim(t))));
        csolve_v3 = @(c) c+psi*aa1-aa2*c^(nu/(nu-1)); % equation 4
        Csim(t) = fsolve(csolve_v3,1);
        
        % Evolution under the argmax q
        if (t<=T-1)
            for i_k = 1:nfine
                for i_x = 1:nx
                    i_s = sub2ind([nx 2],i_x,ssigmaxsim(t));
                    kplus = kopt_fine(i_k,i_s,i_K,i_qmax);
                    
                    % Assign mass to tomorrow's distribution
                    for i_xplus = 1:nx
                        if (kplus>fine_grid(1) && kplus<fine_grid(nfine))
                            lower_ind = find(fine_grid<=kplus,1,'last');
                            upper_ind = lower_ind + 1;
                            denom = fine_grid(upper_ind)-fine_grid(lower_ind);
                            dist_k(lower_ind,i_xplus,t+1) = dist_k(lower_ind,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t)*(kplus-fine_grid(lower_ind))/denom;
                            dist_k(upper_ind,i_xplus,t+1) = dist_k(upper_ind,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t)*(fine_grid(upper_ind)-kplus)/denom;
                        elseif (kplus<=fine_grid(1))
                            dist_k(1,i_xplus,t+1) = dist_k(1,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t);
                        elseif (kplus>=fine_grid(nfine))
                            dist_k(nfine,i_xplus,t+1) = dist_k(nfine,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t);
                        end
                    end
                    
                end
            end
        end
    end
    
    % update kss
    % kss = mean(Ksim(burnin+1:T));
    
    % Regress to get coefficients of K law
    XX = [ones(T-burnin-1,1) log(Ksim(burnin+1:T-1))' log(ssigmax(ssigmaxsim(burnin+1:T-1)))' log(zsim(burnin+1:T-1))'];
    Y = log(Ksim(2+burnin:T)');
    bbeta_K = (XX'*XX)\(XX'*Y);
    e = Y-XX*bbeta_K;
    ytilde = Y-mean(Y);
    Rsq_K = 1-(e'*e)/(ytilde'*ytilde);
    pphi_KC_new = damp*bbeta_K(1)+(1-damp)*pphi_KC; pphi_KK_new = damp*bbeta_K(2)+(1-damp)*pphi_KK;    
    pphi_Kssigmax_new = damp*bbeta_K(3)+(1-damp)*pphi_Kssigmax; pphi_Kz_new = damp*bbeta_K(4)+(1-damp)*pphi_Kz;
    
    % Regress to get q law
    Y = log(qsim(1+burnin:T-1))';
    bbeta_q = (XX'*XX)\(XX'*Y);
    e = Y-XX*bbeta_q;
    ytilde = Y-mean(Y);
    Rsq_q = 1-(e'*e)/(ytilde'*ytilde);
    pphi_qC_new = damp*bbeta_q(1)+(1-damp)*pphi_qC; pphi_qK_new = damp*bbeta_q(2)+(1-damp)*pphi_qK;
    pphi_qssigmax_new = damp*bbeta_q(3)+(1-damp)*pphi_qssigmax; pphi_qz_new = damp*bbeta_q(4)+(1-damp)*pphi_qz;
    % Regress to get q law
    Y = log(Csim(1+burnin:T-1))';
    bbeta_C = (XX'*XX)\(XX'*Y);
    e = Y-XX*bbeta_C;
    ytilde = Y-mean(Y);
    Rsq_C = 1-(e'*e)/(ytilde'*ytilde);
    pphi_CC_new = damp*bbeta_C(1)+(1-damp)*pphi_CC; pphi_CK_new = damp*bbeta_C(2)+(1-damp)*pphi_CK;
    pphi_Cssigmax_new = damp*bbeta_C(3)+(1-damp)*pphi_Cssigmax; pphi_Cz_new = damp*bbeta_C(4)+(1-damp)*pphi_Cz;
    
    % Update mmu_old as well
   diff = norm([pphi_KC,    pphi_KK,    pphi_Kssigmax,    pphi_Kz,    pphi_qC,    pphi_qK,    pphi_qssigmax,    pphi_qz,    pphi_CC,    pphi_CK,    pphi_Cssigmax,    pphi_Cz]-...
                [pphi_KC_new,pphi_KK_new,pphi_Kssigmax_new,pphi_Kz_new,pphi_qC_new,pphi_qK_new,pphi_qssigmax_new,pphi_qz_new,pphi_CC_new,pphi_CK_new,pphi_Cssigmax_new,pphi_Cz_new],Inf);
    
    % Update mmu_old as well
    pphi_KC = pphi_KC_new; pphi_KK = pphi_KK_new; pphi_Kssigmax = pphi_Kssigmax_new; pphi_Kz = pphi_Kz_new;
    pphi_qC = pphi_qC_new; pphi_qK = pphi_qK_new; pphi_qssigmax = pphi_qssigmax_new; pphi_qz = pphi_qz_new;
    pphi_CC = pphi_CC_new; pphi_CK = pphi_CK_new; pphi_Cssigmax = pphi_Cssigmax_new; pphi_Cz = pphi_Cz_new;



    %% Plot something
%     i_x = 4;
%     i_K = 3;
%     i_q = 3;
%     figure
%     plot(k_grid,W_old(:,i_x,i_K,i_q),'-r',k_grid,U_old(:,i_x,i_K),'-b');
%     figure
%     plot(k_grid,V_old(:,i_x,i_K,i_q));
%     figure
%     plot(k_grid,kopt(:,i_x,i_K,i_q));
%     figure
%     plot(q_grid,revenue(:,T),'-.r',q_grid,revenue(:,T-7),'b');
%     figure
%     plot(q_grid,demand(:,T),'-.r',q_grid,demand(:,T-7),'b');
    
    %%
    outer_iter = outer_iter + 1;
    disp_text = sprintf('Rsq_K = %d, Rsq_q = %d, Rsq_C = %d',Rsq_K,Rsq_q,Rsq_C);
    disp(disp_text);
    disp_text = sprintf('log(q) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_qC,pphi_qK,pphi_qssigmax,pphi_qz_new);
    disp(disp_text);
    disp_text = sprintf('log(Kplus) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_KC,pphi_KK,pphi_Kssigmax,pphi_Kz_new);
    disp(disp_text);
    disp_text = sprintf('log(C) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_CC,pphi_CK,pphi_Cssigmax,pphi_Cz_new);
    disp(disp_text);
    disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,diff,iter,err);
    disp(disp_text);
    disp('===============================');
    
end


toc

save main.mat
