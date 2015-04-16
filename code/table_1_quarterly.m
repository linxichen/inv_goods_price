% correlation table
X = xlsread('quarterly investment goods price.xls','B2:I164');
igp = X(:,2:end);% investment goods price measurements and unemployment
N = 10;% number of lag and lead 
[m,n] = size(X);
unc_sale = X(:,1);
table1 = zeros(n-1,2*N+1);
for i = 1:n-1
    table1(i,N+1) = corr(unc_sale,igp(:,i));
    for j = 1:N
        table1(i,N+1+j) = corr(unc_sale(1:end-j),igp(1+j:end,i));% corr(unctainty(t),inv price(t+j))
        table1(i,N+1-j) = corr(unc_sale(1+j:end),igp(1:end-j,i));% corr(unctainty(t),inv price(t-j))
    end
end




