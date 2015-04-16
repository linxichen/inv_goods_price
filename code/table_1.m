% correlation table
X = xlsread('yearly investment goods price.xls','B2:J39');
igp = X(:,3:end);% investment goods price measurements and unemployment
N = 10;% number of lag and lead 
[m,n] = size(X);
unc_tfp = X(:,1);% uncertainty 1
table1 = zeros(n-2,2*N+1);
for i = 1:n-2
    table1(i,N+1) = corr(unc_tfp,igp(:,i));
    for j = 1:N
        table1(i,N+1+j) = corr(unc_tfp(1:end-j),igp(1+j:end,i));% corr(unctainty(t),inv price(t+j))
        table1(i,N+1-j) = corr(unc_tfp(1+j:end),igp(1:end-j,i));% corr(unctainty(t),inv price(t-j))
    end
end
unc_sale = X(:,2);% uncertainty 2
table2 = zeros(n-2,2*N+1);
for i = 1:n-2
    table2(i,N+1) = corr(unc_sale,igp(:,i));
    for j = 1:N
        table2(i,N+1+j) = corr(unc_sale(1:end-j),igp(1+j:end,i));% corr(unctainty(t),inv price(t+j))
        table2(i,N+1-j) = corr(unc_sale(1+j:end),igp(1:end-j,i));% corr(unctainty(t),inv price(t-j))
    end
end

