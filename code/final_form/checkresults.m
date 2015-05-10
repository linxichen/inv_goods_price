%%============
% Check thing that should be checked before declare successful convergence.
%==============

%% Distrubtion of firms
if sum((sum(dist_k(:,:,end),2) > 0)) <= 1
    warning('Marginal distribution of capital is degenerate. Check grid for capital and investment good prices.');
else
    disp('================================================================')
    disp('Marginal Distribution of capital is not degenerate. We are good.')
    disp('================================================================')
end

%% Demand is downward sloping?
t = T;
[~,i_K] = min(abs((Ksim(t)-K_grid))); 
profit_check = zeros(1,nq);
demand_check = profit_check;
whichxind = zeros(1,nx);
for i_x = 1:nx
    whichxind(i_x) = sub2ind([nz nx 2],zindsim(t),i_x,ssigmaxsim(t));
end
for i_q = 1:nq
    profit_check(i_q) = sum(vec(tot_profit_grid(:,:,i_q)));
    demand_check(i_q) = sum(demand_grid(:));
end

if sum(diff(demand_check) > 0) ~= 0
    warning('Demand curve is not downward sloping in the end. Check grid for capital and investment good prices.');
else
    disp('================================================================')
    disp('Demand curve seems to be downward sloping.')
    disp('================================================================')
end

%% Is capital grid big enough?
if min(koptind(:)) == 1
    warning('Somebody chooses the smallest capital level in the grid. Decrease min(k_grid).');
else
    disp('================================================================')
    disp('Nobody is choosing the smallest capital level. We are good.')
    disp('================================================================')
end
if max(koptind(:)) == length(k_grid)
    warning('Somebody chooses the largest capital level in the grid. Increase max(k_grid).');
else
    disp('================================================================')
    disp('Nobody is choosing the largest capital level. We are good.')
    disp('================================================================')
end

%% 