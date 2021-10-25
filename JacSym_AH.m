%variables:
syms  c1 c2 c3 c4 c5 c6 c7 V real       

%parameters:
syms kT2 k1 k2 K1 K2 c6_star real

x = [c1, c2, c3, c4, c5, c6, c7, V];

par = [kT2, k1, k2, K1, K2, c6_star];

eq = odefun_AH(0,x,par);

symjac = jacobian(eq,x)