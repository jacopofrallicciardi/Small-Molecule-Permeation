function [AH,A] = HHE(c,pH,pKa)

% Known the concentration, pKa and pH of the solution calculate the
% concentration of A+ and AH using Henderson-Hasselbach equation

r = 10^(pH-pKa);

A = c*(r./(1+r));
AH = c*(1./(1+r)); 

end