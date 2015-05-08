%% Housekeeping
% 05/05/2015 Linxi: Add intensive margin to simplified version. This should
% be the only version we use from now. Heavily comment this one.
% Key assumptions:
% 1. resale price is 0. Means irreversible capital. DONT CHANGE THIS!
% otherwise have to rewrite a lot of codes
% 2. wage = MC = 1
% 3. Partial eqm, risk neutral
% 4. Uses calibration file

%% Housekeeping
clear; close all; clc;
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

%% Accuracy control
nk = 50; % number of grid points on capital stock
nfine = 150; % number of grid points for policy functions and simulation
nx = 7; % number of grid points on idiosyncractic prod.
nz = 3; % number of grid points on aggregate productivity
ns = nx*nz*2; % total dimension of exo state, (idio, agg, ssigmax)
nK = 21; % agg capital state
nq = 15; % number of price points inv producer can set
m = 3; % support is m s.d. away from mean
tol = 1e-3; % when to stop VFI and KS
maxiter = 2000;
damp = 0.3; % help with outer convergence

%% Grids
[Z,PZ] = tauchen(nz,0,rrhoz,ssigmaz,m); Z = exp(Z); % agg TFP
ssigmax = [ssigmax_low,ssigmax_high];
[X,PX_high] = tauchen(nx,0,rrhox,ssigmax_high,m);
%=========================================================================%
% This PX_high(i,j) gives the approximate prob of P(X_t+1(j) | X_t = X(i))
% Follows X_t+1 = rrhox*X(i) + ssigma_high* N(0,1)
%=========================================================================%
PX_low = tauchen_givengrid(0,rrhox,ssigmax_low,X); X = exp(X);
P = zeros(nx*nz*2,nx*nz*2);  % The transition of exo state vec, all independent
low = 1;
high = 2;
uncond_X_low = PX_low^30000;
uncond_X_high = PX_high^30000;
% Below find the transition prob of three indie variable
for i = 1:ns
    [i_x,i_z,i_ssigmax] = ind2sub([nx nz 2],i);
    for j = 1:ns
        [j_x,j_z,j_ssigmax] = ind2sub([nx nz 2],j);
        P(i,j) = PZ(i_z,j_z)*Pssigmax(i_ssigmax,j_ssigmax)*( (j_ssigmax==low)*PX_low(i_x,j_x)+ (j_ssigmax==high)*PX_high(i_x,j_x) );
    end
end

%=========================================================================%
% inv price grid. Since MC is 1, the price must be larger than 1
% From the two period model, q = (xxi - bbeta*(1-ddelta))/(xxi - 1),
% The tail index has to be larger than 1/(1-aalpha) = 1.4286, this implies
% q = 1.05.
%=========================================================================%
qmin = 1.001;
qmax = 1.2;
q_grid = linspace(qmin,qmax,nq); 
% Capital grid
[Kmax_lower,Kmin_upper,Kmin_lower]=calibration_ivg_v2(X,PX_low,q_grid);
min_k = 40;
max_k = Kmax_lower;
k_grid = zeros(nk,1);
for i_k = 1:nk
    k_grid(i_k) = max_k*(1-ddelta)^(nk-i_k);
end
fine_grid = linspace(k_grid(1),k_grid(nk),nfine)';
noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
for i_k = 1:nk
    [~,noinvest_ind(i_k)] = min(abs(k_grid-(1-ddelta)*k_grid(i_k)));
end
noinvest_ind_fine = ones(nfine,1); % for each k, the index of tmr k if no invest
for i_k = 1:nfine
    [~,noinvest_ind_fine(i_k)] = min(abs(fine_grid-(1-ddelta)*fine_grid(i_k)));
end
ppsi_n = 1; % labor distutility
w = ppsi_n; % Normalized to one, w
mc = w; % mc of investment producer is also 1
inv_mat = repmat(k_grid',nk,1)-(1-ddelta)*repmat(k_grid,1,nk);
%=========================================================================%
% This inv_mat(i,j) gives inv amount needed to reach k_t+1 = k(j) when 
% k_t = k(i)
%=========================================================================%
pos_inv = inv_mat>0;
neg_inv = inv_mat<=0;
inv_mat_fine = repmat(fine_grid',nfine,1)-(1-ddelta)*repmat(fine_grid,1,nfine);
pos_inv_fine = inv_mat_fine>0;
neg_inv_fine = inv_mat_fine<=0;

K_grid = linspace(k_grid(1),k_grid(nk),nK)'; % Aggregate capital grid


%% Aggregate rules and its parameters
% Rule of q as function of aggregate states
pphi_qC = log(mean(q_grid)); % constant term
pphi_qK = 0.0; % w.r.t agg K
pphi_qz = 0.0; % w.r.t agg TFP
pphi_qssigmax = 0.00; % w.r.t uncertainty
pphi_q = [pphi_qC,pphi_qK,pphi_qz,pphi_qssigmax];

pphi_KK = 0.99; 
pphi_KC = log(mean(K_grid)); 
pphi_Kz = 0.01;
pphi_Kssigmax = 0.01;% Aggregate Law of motion for aggregate capital
pphi_K = [pphi_KC,pphi_KK,pphi_Kz,pphi_Kssigmax];

pphi_CK = 0.99; 
pphi_CC = 0.1; 
pphi_Cz = 0.01;
pphi_Cssigmax = 0.01;% Aggregate Law of motion for consumption
pphi_C = [pphi_CC,pphi_CK,pphi_Cz,pphi_Cssigmax];
load aggrules.mat;

%% Initialize value functions
if (exist('valuefunctions.mat','file') == 2)
    load valuefunctions.mat
end
if isequal(size(W_old),[nk,ns,nK,nq])
else
    W_old = ones(nk,ns,nK,nq); % value of matching with investment goods producer after paying the search cost
    W_new = W_old;
    U_old = ones(nk,ns,nK); % value of not going to search, not related to current q
    U_new = U_old;
    V_old = ones(nk,ns,nK,nq); % maximized value after discrete choice
    V_new = ones(nk,ns,nK,nq); %
end


% profit_fine = zeros(nfine,ns,nq);
V_new_fine = zeros(nfine,ns,nK,nq);
EV_new_fine = zeros(nfine,ns,nK,nq);
W_new_fine = zeros(nfine,ns,nK,nq);
U_new_fine = zeros(nfine,ns,nK);
koptind_active = zeros(nk,ns,nK,nq);
koptind = zeros(nk,ns,nK,nq);
kopt_fine = zeros(nfine,ns,nK,nq);
active = zeros(nk,ns,nK,nq);

% Compute operation profit at all possible idio state, just Y - w*l
profit = zeros(nk,ns); 
L = zeros(nk,ns);
for i_k = 1:nk
    for i_s = 1:ns
        [i_z,i_x,~] = ind2sub([nz nx 2],i_s);
        L(i_k,i_s) = (w*k_grid(i_k)^(-aalpha)/Z(i_z)/X(i_x)/v)^(1/(v-1));
        profit(i_k,i_s) = Z(i_z)*X(i_x)*k_grid(i_k)^aalpha*L(i_k,i_s)^v - w*L(i_k,i_s);
    end
end

% Compute net profit, i.e. operation profit minus investment cost
netprofit = zeros(nk,nk,ns,nq);
for i_k = 1:nk
    for i_kplus = 1:nk
        for i_s = 1:ns
            for i_q = 1:nq
                netprofit(i_k,i_kplus,i_s,i_q) = profit(i_k,i_s)-q_grid(i_q)*(k_grid(i_kplus)-(1-ddelta)*k_grid(i_k))*(pos_inv(i_k,i_kplus)+neg_inv(i_k,i_kplus)*pphi);
            end
        end
    end
end

% Prepare for Simulation stuff
T = 1000; % How long to simulate
burnin = ceil(0.1*T); % how many periods to discard to get rid of dependence on initial state
kss = mean(k_grid);
Ksim = kss*ones(1,T);
qsim = mean(q_grid)*ones(1,T);
revenue = zeros(nq,1); % revenue of inv producer
demand = zeros(nq,T);
dist_k = zeros(nfine,nx,T); % !! Should change this to Pareto
dist_k(ceil(nfine/2),:,:) = ones(nx,T)/nx; 

% New Discrete Simulation? simulate the path of z and ssigmax
x_cdf_low = cumsum(PX_low,2);
x_cdf_high = cumsum(PX_high,2);
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
outer_diff = 10;
outer_iter = 0;

while ((outer_diff > tol) && (outer_iter < maxiter))
    %============ VFI Begins==============================================%
    % Inner loop
    err = 10;
    iter = 0;
    while err > tol
        for i_s = 1:ns
            [i_z,i_x,i_ssigmax] = ind2sub([nz nx 2],i_s);
            for i_K = 1:nK
                % What is agg state today
                log_aggstate = [1; log(K_grid(i_K)); log(ssigmax(i_ssigmax)); log(Z(i_z))];
                % Forecast future values
                [~,i_Kplus] = min(abs(K_grid-exp(pphi_K*log_aggstate)));
                [~,i_qplus] = min(abs(q_grid-exp(pphi_q*log_aggstate)));
                % C = exp(pphi_CC+pphi_CK*log(K_grid(i_K))+pphi_Cssigmax*log(i_ssigmax)+pphi_Cz*log(Z(i_z)));
                                
                % Find expected value across i_splus
                EV_noinvest = V_old(noinvest_ind,:,i_Kplus,i_qplus)*P(i_s,:)'; % EV_oninvest(i) gives expected V given k_t+1 = K(i), optimally!
                
                % U has no choice, so no need to find max in a loop.
                U_new(:,i_s,i_K) = profit(:,i_s) + bbeta*EV_noinvest;
                
                % Find W
                EV_invest = V_old(:,:,i_Kplus,i_qplus)*P(i_s,:)';
                %=========================================================%
                % This EV_invest(i) is expected value when kplus=K(i)
                % Not optimally! We need to find maximizing kplus
                %=========================================================%
                for i_q = 1:nq
                    candidate = netprofit(:,:,i_s,i_q) + bbeta*repmat(EV_invest',nk,1);
                    %=====================================================%
                    % candidate(i_k,i_kplus) is the thing inside max
                    % operator. expected value tomorrow depends only on
                    % capital tomorrow and aggregate state, therefore I use
                    % repmat here. netprofit doesn't depend on value
                    % function so it is pre-computed outside for loops
                    %=====================================================%
                    [W_new(:,i_s,i_K,i_q),koptind_active(:,i_s,i_K,i_q)] = max(candidate,[],2);
                    V_new(:,i_s,i_K,i_q) = max(W_new(:,i_s,i_K,i_q),U_new(:,i_s,i_K));
                end
            end
        end
        
        err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
        V_old = V_new;
        W_old = W_new;
        U_old = U_new;
        iter = iter + 1;
        if mod(iter,1) == 0
            disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,outer_diff,iter,err);
            disp(disp_text);
        end
    end
    %============ VFI Ends================================================%
    
    active = W_new > repmat(U_new,1,1,1,nq);
    koptind = repmat(noinvest_ind,1,ns,nK,nq).*(1-active) + active.*koptind_active;
    kopt = k_grid(koptind);
    plot(k_grid,kopt(:,sub2ind([nx nz 2],ceil(nx/2),ceil(nz/2),1),ceil(nK/2),ceil(nq/2))-(1-ddelta)*k_grid)
    save('valuefunctions.mat','V_new','W_new','U_new','V_old','W_old','U_old')
    
    % Interpolate on finer grid
    for i_q = 1:nq
        for i_K = 1:nK
            for i_s = 1:ns
                U_new_fine(:,i_s,i_K) = interp1(k_grid,U_new(:,i_s,i_K),fine_grid,'linear')';
                W_new_fine(:,i_s,i_K,i_q) = interp1(k_grid,W_new(:,i_s,i_K,i_q),fine_grid,'linear')';
                kopt_fine(:,i_s,i_K,i_q) = interp1(k_grid,kopt(:,i_s,i_K,i_q),fine_grid,'linear')';
            end
        end
    end    
    active_fine = W_new_fine > repmat(U_new_fine,1,1,1,nq);
    plot(fine_grid,kopt_fine(:,sub2ind([nx nz 2],ceil(nx/2),ceil(nz/2),1),ceil(nK/2),ceil(nq/2))-(1-ddelta)*fine_grid)

    
    %% Given individual policies, simulate a large panel to update aggregate law of motion
    % Zero out after each outer iteration.
    dist_k(:,:,2:T) = zeros(nfine,nx,T-1);
    for t = 1:T
        % Find aggregate stuff today
        Ksim(t) = sum(vec(dist_k(:,:,t).*repmat(fine_grid,1,nx)));
        [~,i_K] = min(abs(Ksim(t)-K_grid));
        % Given uncertainty and agg TFP, find the agg exo state index of each indiviual
        % prod. level. He uses this to forecast expected value function
        whichs = zeros(1,nx);
        for i_x = 1:nx
            whichs(i_x) = sub2ind([nz nx 2],zindsim(t),i_x,ssigmaxsim(t));
        end
        
        % Find transition matrix according to today's state
        if ssigmaxsim(t) == low
            whichprob = PX_low;
        elseif ssigmaxsim(t) == high
            whichprob = PX_high;
        end
       
        % According to policy functions, find the optimal q
        for i_q = 1:nq
            tot_profit_grid = (q_grid(i_q)-mc)*dist_k(:,:,t).*(kopt_fine(:,whichs,i_K,i_q)-(1-ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
            demand_grid = dist_k(:,:,t).*(kopt_fine(:,whichs,i_K,i_q)-(1-ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
            revenue(i_q,t) = sum(tot_profit_grid(:));
            demand(i_q,t) = sum(demand_grid(:));
        end
        [~,i_qmax] = max(revenue(:,t));
        qsim(t) = q_grid(i_qmax);
        
        % Evolution under the argmax q
        if (t<=T-1)
            for i_k = 1:nfine
                for i_x = 1:nx
                    ind_z =find(Z==zsim(t));
                    i_s = sub2ind([nx nz 2],i_x,ind_z,ssigmaxsim(t));
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
    
    %% Regress to get coefficients of K law
    XX = [ones(T-burnin-1,1) log(Ksim(burnin+1:T-1))' log(ssigmax(ssigmaxsim(burnin+1:T-1)))' log(zsim(burnin+1:T-1))'];
    Y = log(Ksim(2+burnin:T)');
    bbeta_K = (XX'*XX)\(XX'*Y);
    e = Y-XX*bbeta_K;
    ytilde = Y-mean(Y);
    Rsq_K = 1-(e'*e)/(ytilde'*ytilde);
    pphi_K_new = damp*bbeta_K' + (1-damp)*pphi_K;
  
    % Regress to get q law
    Y = log(qsim(1+burnin:T-1))';
    bbeta_q = (XX'*XX)\(XX'*Y);
    e = Y-XX*bbeta_q;
    ytilde = Y-mean(Y);
    Rsq_q = 1-(e'*e)/(ytilde'*ytilde);
    pphi_q_new = damp*bbeta_q' + (1-damp)*pphi_q;
  
    % Update mmu_old as well
    outer_diff = norm([pphi_K,pphi_q]-[pphi_K_new,pphi_q_new],Inf);
    
    % Update mmu_old as well
    pphi_K = pphi_K_new;
    pphi_q = pphi_q_new;
   
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
    disp_text = sprintf('Rsq_K = %d, Rsq_q = %d',Rsq_K,Rsq_q);
    disp(disp_text);
    disp_text = sprintf('log(q) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_q);
    disp(disp_text);
    disp_text = sprintf('log(Kplus) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_K);
    disp(disp_text);
    disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,outer_diff,iter,err);
    disp(disp_text);
    disp('===============================');
    save('aggrules.mat','pphi_K','pphi_q');
    save('Rsq.mat','Rsq_K','Rsq_q');
    
end


toc



%% Decompose Demand Curve
% According to policy functions, find the optimal q
% Fix the distribution of firms the last period T
% t = T;
% [~,i_K] = min(abs((Ksim(t)-K_grid))); 
% revenue_lowtfp = zeros(1,nq);
% demand_lowtfp = revenue_lowtfp;
% whichs = zeros(1,nx);
% for i_x = 1:nx
%     whichs(i_x) = sub2ind([nz nx 2],1,i_x,ssigmaxsim(t));
% end
% for i_q = 1:nq
%     tot_profit_grid(:,:,i_q) = (q_grid(i_q)-psi)*dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     % tot_revenue_grid(tot_revenue_grid<0) = 0;
%     demand_grid = dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     revenue_lowtfp(i_q) = sum(vec(tot_profit_grid(:,:,i_q)));
%     demand_lowtfp(i_q) = sum(demand_grid(:));
% end
% figure
% mesh(X,fine_grid,tot_profit_grid(:,:,1))
% 
% t = T;
% [~,i_K] = min(abs((Ksim(t)-K_grid))); 
% revenue_hightfp = zeros(1,nq);
% demand_hightfp = revenue_hightfp;
% whichs = zeros(1,nx);
% for i_x = 1:nx
%     whichs(i_x) = sub2ind([nz nx 2],nz,i_x,ssigmaxsim(t));
% end
% for i_q = 1:nq
%     tot_profit_grid(:,:,i_q) = (q_grid(i_q)-psi)*dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     % tot_revenue_grid(tot_revenue_grid<0) = 0;
%     demand_grid = dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     revenue_hightfp(i_q) = sum(vec(tot_profit_grid(:,:,i_q)));
%     demand_hightfp(i_q) = sum(demand_grid(:));
% end
% 

figure
mesh(X,fine_grid,tot_profit_grid(:,:,1))

plot(q_grid,revenue_lowtfp,'b',q_grid,revenue_hightfp,'r')
save main.mat

checkresults;
