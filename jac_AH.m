function j = jac_AH(t,x,par)

% Variables
% x(1) = c1 [H_2O]
% x(2) = c2 [Cal/Pyr]
% x(3) = c3 [H_2PO_4-]
% x(4) = c4 [HPO_42-]
% x(5) = c5 [H+]
% x(6) = c6 [AH]
% x(7) = c7 [A+]
% x(8) = V

c1 = x(1);
c2 = x(2);
c3 = x(3);
c4 = x(4);
c5 = x(5);
c6 = x(6);
c7 = x(7);
V = x(8);

% Parameters
kT2 = par(1);
k1 = par(2);
k2 = par(3);
K1 = par(4);
K2 = par(5);
c6_star = par(6);

j = [...
[ 0,                                            0,                                                 0,                                                         0,                                                                      0,                                                         0,                                                         0,                                                                     0]
[ 0, - c2/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V,                                             -c2/V,                                                     -c2/V,                                                                  -c2/V,                                                     -c2/V,                                                     -c2/V,                            (c2*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                        -c3/V, - k1 - c3/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V,                                         (c5*k1)/K1 - c3/V,                                                      (c4*k1)/K1 - c3/V,                                                     -c3/V,                                                     -c3/V,                            (c3*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                        -c4/V,                                         k1 - c4/V, - c4/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V - (c5*k1)/K1,                                                    - c4/V - (c4*k1)/K1,                                                     -c4/V,                                                     -c4/V,                            (c4*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                        -c5/V,                                         k1 - c5/V,                                       - c5/V - (c5*k1)/K1, - c5/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V - (c4*k1)/K1 - (c7*k2)/K2,                                                 k2 - c5/V,                                       - c5/V - (c5*k2)/K2,                            (c5*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                        -c6/V,                                             -c6/V,                                                     -c6/V,                                                      (c7*k2)/K2 - c6/V, - k2 - c6/V - kT2/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V,                                         (c5*k2)/K2 - c6/V, (kT2*(c6 - c6_star))/V^2 + (c6*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                        -c7/V,                                             -c7/V,                                                     -c7/V,                                                    - c7/V - (c7*k2)/K2,                                                 k2 - c7/V, - c7/V - (c2 + c3 + c4 + c5 + c6 + c7 - 1)/V - (c5*k2)/K2,                            (c7*(c2 + c3 + c4 + c5 + c6 + c7 - 1))/V^2]
[ 0,                                            1,                                                 1,                                                         1,                                                                      1,                                                         1,                                                         1,                                                                     0]
 ];
end
 