function y = dividelabor(x,aalpha,bbeta,ggamma,ddelta,A)
Nplusone = length(x); % N firms equations plus one equation.
y = zeros(Nplusone,1);
for i_eq = 1:Nplusone-1
    y(i_eq) = x(Nplusone) - ggamma*A(i_eq)^(1/(1-aalpha))*((1/bbeta-1+ddelta)/aalpha)^(aalpha/(aalpha-1))*x(i_eq)^(ggamma*aalpha/(1-aalpha)+ggamma-1);
end
y(Nplusone) = 1 - sum(x(1:Nplusone-1));
end