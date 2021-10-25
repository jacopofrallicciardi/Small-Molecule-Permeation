clear all

global INPUT x0 cs_star k1 k2 K1_star K2_star c6_star pKa1 SVR

% IMPORTANT!!! %
% 1L = 10^3 cm^3
% 1M = mol/L = 10^-3 mol/cm^3
% 1 nm = 10^-7 cm

%% - INPUT OF RELAXATION CURVES - %%

% Input of the first set of data
[filename1, pathname1] = uigetfile('*.dat*','PYRANINE (KCl)', 'MultiSelect','off');

INPUT = importdata(filename1,'\t');

%% PHYSICAL PARAMETERS - BUFFER  %%

% pKa of KPi
pKa1 = 7.21;
% Kpi dissociation equilibrium constant
K1 = 10^-3*10^-pKa1; % [mol/cm^3]

% Starting external pH @ t = t0
pH0_out = 6;
% External proton concentration
H_out = (10^-pH0_out)*10^-3;    %[mol/cm^3]

% External KPi concentration
c_KPi_out = 100E-6; %[mol/cm^3]

%pH0_in = 6.4;
% Internal KPi concentration
c_KPi_in = 100E-6; %[mol/cm^3]
% Calculation of the KPi species
%[KPH,KP] = HHE(c_KPi_in,pH0_in,pKa1);

% Water molar volume
Mw = 18;   % [cm^2/mol]

%% PHYSICAL PARAMETERS - AH  %%

% KPi and AH dissocation rate constant
k1 = 1E6;   %[1/s]
k2 = 1.01E6;   %[1/s]

% pKa of the acid AH
pKa2 = input('Set the pKa of the acid: ');
% Kpi dissociation equilibrium constant
K2 = 10^-3*10^-pKa2; % [mol/cm^3]

% Acid concentration
c_AH = 100; % [mM]
c_AH = c_AH*10^-6; %[mol/cm^3]
% Concentration of AH and A+
[AH,A] = HHE(c_AH,pH0_out,pKa2);
c6_star = AH;

% External concentration of the osmolite (e.g. AH)
osm = 2*c_AH;    %[mol/cm^3]

% External solute concentration
cs_star = (c_KPi_out + osm + H_out);    %[mol/cm^3]

%% PHYSICAL PARAMETERS - VESICLE DIMENSION AND PERMEABILITY %%

% Starting radius, surface area and volume of the vesicle
S = 10^-6;    %[cm^2]
V0 = 10^-10; %[cm^3]
SVR = S/V0;

% Permeabilities
Pw = 0.0;      %[cm/s] H2O
P2 = 0.0003;    %[cm/s] AH

%% STARTING CONDITIONS %%

% Variables
% x(1) = c1 [H_2O]
% x(2) = c2 [Cal/Pyr]
% x(3) = c3 [H_2PO_4-]
% x(4) = c4 [HPO_42-]
% x(5) = c5 [H+]
% x(6) = c6 [AH]
% x(7) = c7 [A+]
% x(8) = V

% Starting conditions
x0(1) = 55E-3;      %[mol/cm^3]
x0(2) = 0E-6;      %[mol/cm^3]
x0(3) = 0;      %[mol/cm^3]
x0(4) = 0;       %[mol/cm^3]
x0(5) = 0;  %[mol/cm^3]
x0(6) = 0;     %[mol/cm^3]
x0(7) = 0;      %[mol/cm^3]
x0(8) = V0;      %[cm^3]

% Scaling of the equilibrium constants
K1_star = K1/cs_star;
K2_star = K2/cs_star;
% Scaling of the external AH concentration
c6_star = c6_star/cs_star;

%% - FIT OF THE DATA - %%

% I perform a fit of the KCl data

% PARS: 
% 1) P2 [cm/s]
% 2) KPi_in [mol/cm^3]

START_GUESS = [P2 c_KPi_in];

[BEST_PARS,errs,chi2,errmatrix] = fminuit('ODE_AH_fitfun_pH',START_GUESS,'runmode');

%%

par = BEST_PARS;
P2 = par(1);
KPi_out = par(2)*10^6;


