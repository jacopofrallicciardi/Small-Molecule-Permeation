function chi2 = ODE_AH_fitfun_Fr(par)

global INPUT INPUT_pR x0 cs_star Mw k1 k2 K1 K2 c6_star

% Experimental relaxation curves
texp = INPUT(:,1);
Fexp = INPUT(:,2);

% Experimental vesicle size distribution
xr = INPUT_pR(:,1)';
pR = INPUT_pR(:,2)';
% Area normalization
pR = pR/sum(pR);

% PARS: 
% 1) K_SV [cm^3/mol]
% 2) Pw [cm/s]
% 3) P2 [cm/s]

% Fitting parameters
K_SV = par(1);
Pw = par(2);
P2 = par(3);

for i = 1 : length(pR)  
    
% Starting radius, surface area and volume of the vesicle
r = xr(i)*10^-7;        %[cm]
S = 4*pi*r^2;    %[cm^2]
V0(i) = (4/3)*pi*r^3; %[cm^3]

x0(8) = V0(i);      %[cm^3]

% SCALING %

% We scale the ODE system to get dimensionless equations
% Concentrations are scaled by the external solute concentration (cs_star)
% The volume is scaled by V0
% The time is scaled by the characteristic time tc that we set equal to the
% volume transport rate kv: tc = kv

% Note that all quantities are dimensionless!!!

% Scaling of the starting conditions 
x0_bar = x0/cs_star;
x0_bar(8) = 1;  % The volume is indeed normalized to V0

% Definition of the transport rate constants
kv(i) = 3*Pw*Mw*cs_star/r;    % Volume transport
kT2(i) = 3*P2/r;              % AH transport

% Scaling of the rate constants with respect to kv 
kT2_bar(i) = kT2(i)/kv(i);
k1_bar(i) = k1/kv(i);
k2_bar(i) = k2/kv(i);

% Dimensionless time
tspan = logspace(-6,4,1000)'; 

% Scaled parameters of the ODE system
parODE = [kT2_bar(i), k1_bar(i), k2_bar(i), K1, K2, c6_star];

%% ODE SOLUTION %%

% Options of the ODE solver.
% We give the jacobian matrix as input
options = odeset('Stats','on','Jacobian',@(t,x)jac_AH(t,x,parODE),'NonNegative',1);

% ODE solution
[t,x] = ode15s(@(t,x) odefun_AH(t,x,parODE), tspan, x0_bar, options);

%% RESCALING %%

% We rescale back to the original dimensions

% Time
td(:,i) = t./kv(i);                 %[s]

% Internal concetrations 
x = x*cs_star; 
c2(:,i) = x(:,2)*10^6;   %[mM]
c5(:,i) = x(:,5)*10^3;   %[M]

end
% I create matrices with the size distribution g(r_0) for integration
pRm = repmat(pR,length(texp),1);
xrm = repmat(xr,length(texp),1);


for j = 1 : length(pR)
    c2q(:,j) = interp1(td(:,j),c2(:,j),texp,'spline');
    c5q(:,j) = interp1(td(:,j),c5(:,j),texp,'spline');
end

%% - AVERAGE F_ratio: CALCEIN READOUT - %%

% Calculation of the F_ratio matrix
for j = 1 : length(pR) 
    F_ratio(:,j) = (1 + K_SV*c2q(1,j)*10^-6) ./ (1 + K_SV*c2q(:,j)*10^-6);
end

% Average <r^3>
r3_av = sum(xr*3.*pR); %[nm^3]

% Calculation of ensemble average F-ratio
Fth = sum(F_ratio.*xrm*3.*pRm,2)/r3_av;

resF = Fexp - Fth;

%%

% Normalized volume and calcein fluorescence ratio
h = figure(1);
clf(h)

axes('position',[0.1 0.2 0.8 0.7])
semilogx(texp,Fexp,'.k',texp,Fth,'-r','LineWidth',2)
ylabel('F(t)/F(0)','fontsize',14, 'FontWeight','bold')
xlim([0.001 600])
allAxes = findall(0,'type','axes');
set(allAxes, 'linewidth', 2)
set(gca,'fontsize',14, 'FontWeight','bold')
set(gca,'XTick',[])

axes('position',[0.1 0.1 0.8 0.1])
semilogx(texp,resF,'k','LineWidth',2)
box on
xlim([0.001 600])
xlabel('time [s]','fontsize',16, 'FontWeight','bold')
ylabel('res.','fontsize',14, 'FontWeight','bold')
set(gca,'fontsize',14, 'FontWeight','bold')
allAxes = findall(0,'type','axes');
set(allAxes, 'linewidth', 2)

chi2 = sum(resF.^2);

