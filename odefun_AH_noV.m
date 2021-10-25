function dx = odefun_AH_noV(t,x,par)

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
k1 = par(1);
k2 = par(2);
K1 = par(3);
K2 = par(4);
c6_star = par(5);

dx1 = 0;
dx2 = 0;
dx3 = - k1*( c3 - (c4*c5)/K1 );
dx4 = + k1*( c3 - (c4*c5)/K1 );
dx5 = + k1*( c3 - (c4*c5)/K1 ) + k2*( c6 - (c7*c5)/K2 );
dx6 = - k2*( c6 - (c7*c5)/K2 ) + (c6_star - c6);
dx7 = + k2*( c6 - (c7*c5)/K2 ) ;
dx8 = 0; 

dx = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8];
    
end
