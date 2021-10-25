function chi2 = ODE_AH_fitfun_pH(par)

global INPUT x0 cs_star k1 k2 K1_star K2_star c6_star pKa1 SVR

% Experimental relaxation curves
texp = INPUT(:,1);
pHexp = INPUT(:,2);

% PARS: 
% 1) P2 [cm/s]
% 2) KPi_in [mol/cm^3]

% Fitting parameters
P2 = par(1);
c_KPi_in = par(2);

pH0_in = 6.4;
% Calculation of the KPi species
[KPH,KP] = HHE(c_KPi_in, pH0_in, pKa1);

x0(3) = KPH;      %[mol/cm^3]
x0(4) = KP;       %[mol/cm^3]
x0(5) = (10^-pH0_in)*10^-3; %[mol/cm^3]

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
beta6 = SVR*P2;              % AH transport

% Scaling of the rate constants with respect to kv 
k1_bar = k1/beta6;
k2_bar = k2/beta6;

% Dimensionless time
tspan = (0:1:10000)';

% Scaled parameters of the ODE system
parODE = [k1_bar, k2_bar, K1_star, K2_star, c6_star];

%% ODE SOLUTION %%

% Options of the ODE solver.
% We give the jacobian matrix as input
options = odeset('Stats','on','Jacobian',@(t,x)jac_AH_noV(t,x,parODE),'NonNegative',1);

% ODE solution
[t,x] = ode15s(@(t,x) odefun_AH_noV(t,x,parODE), tspan, x0_bar, options);

%% RESCALING %%

% We rescale back to the original dimensions

% Time
td = t./beta6;                 %[s]

% Internal concetrations 
x = x*cs_star; 
c5 = x(:,5)*10^3;   %[M]

c5q = interp1(td,c5,texp,'spline');

%% - pH-PROBE READOUT - %%

% Calculation of the matrix pH(t) from the computed relaxation curves c5(r,t)
pHth = -log10(c5q);

respH = pHexp - pHth;

%%

% Normalized volume and calcein fluorescence ratio
h = figure(1);
clf(h)

axes('position',[0.15 0.2 0.75 0.75])
semilogx(texp,pHexp,'.k',texp,pHth,'-r','LineWidth',2)
ylabel('pH(t)','fontsize',14, 'FontWeight','bold')
xlim([0 10000])
allAxes = findall(0,'type','axes');
set(allAxes, 'linewidth', 2)
set(gca,'fontsize',14, 'FontWeight','bold')
set(gca,'XTick',[])

axes('position',[0.15 0.15 0.75 0.1])
semilogx(texp,respH,'k','LineWidth',2)
box on
xlim([0 10000])
xlabel('time [s]','fontsize',16, 'FontWeight','bold')
ylabel('res.','fontsize',14, 'FontWeight','bold')
set(gca,'fontsize',14, 'FontWeight','bold')
allAxes = findall(0,'type','axes');
set(allAxes, 'linewidth', 2)

chi2 = sum(respH.^2);

