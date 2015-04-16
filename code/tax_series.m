X = xlsread('tax and income.xls','B2:L176');
CT = X(:,1);
CP = X(:,2);
RI = X(:,3);
PRI = X(:,4);
NI = X(:,5);
FIT = X(:,6);
SIT = X(:,7);
EC = X(:,8);
W = X(:,9);
PT = X(:,10);
CSI = X(:,11);
% capital income
CI = PRI/2+RI+CP+NI;
% personal income tax
tau_P = (FIT+SIT)./(W+PRI/2+CI);
% labor tax
tau_L = (tau_P.*(W+PRI/2)+CSI)./(EC+PRI/2);
% capital tax
tau_K = (tau_P.*CI+CT+PT)./(CI+PT);
