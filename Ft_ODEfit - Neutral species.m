clear all

global INPUT INPUT_pR x0 cs_star Mw k1 k2 K1 K2 c6_star

% IMPORTANT!!! %
% 1L = 10^3 cm^3
% 1M = mol/L = 10^-3 mol/cm^3
% 1 nm = 10^-7 cm

%% - INPUT OF RELAXATION CURVES - %%

% Input of the first set of data
[filename1, pathname1] = uigetfile('*.dat*','CALCEIN', 'MultiSelect','off');

INPUT = importdata(filename1,'\t');

%% - INPUT OF SIZE DISTRIBUTION - %%

% Input of the first set of data
[filename2, pathname2] = uigetfile('*.txt*','SIZE DISTRIBUTION', 'MultiSelect','off');

INPUT_pR = importdata(filename2,'\t');

%% PHYSICAL PARAMETERS - BUFFER  %%

% pKa of KPi
pKa1 = 7.21;
% Kpi dissociation equilibrium constant
K1 = 10^-3*10^-pKa1; % [mol/cm^3]

% Starting external and internal pH @ t = t0
pH0 = 7;
% External proton concentration
H_out = (10^-pH0)*10^-3;    %[mol/cm^3]

% Internal KPi concentration
c_KPi = 90E-6; %[mol/cm^3]
% External KPi concentration
c_KPi_out = 100E-6; %[mol/cm^3]
% Calculation of the KPi species
[KPH,KP] = HHE(c_KPi,pH0,pKa1);

% Water molar volume
Mw = 18;   % [cm^2/mol]

% Stern-Volmer dynamic quenching constant
K_SV = 0.1E6; %[cm^3/mol]

%% PHYSICAL PARAMETERS - AH  %%

% KPi and AH dissocation rate constant
k1 = 1E6;   %[1/s]
k2 = 1.01E6;   %[1/s]

% pKa of the acid AH
pKa2 = input('Set the pKa of the acid: ');
% Kpi dissociation equilibrium constant
K2 = 10^-3*10^-pKa2; % [mol/cm^3]

% Acid concentration
c_AH = input('Set the osmolite concentration [mM]: ');
c_AH = c_AH*10^-6; %[mol/cm^3]
% Concentration of AH and A+
[AH,A] = HHE(c_AH,pH0,pKa2);
c6_star = AH;

% External concentration of the osmolite (e.g. AH)
osm = c_AH;    %[mol/cm^3]

% External solute concentration
cs_star = (c_KPi_out + osm + H_out);    %[mol/cm^3]

%% PHYSICAL PARAMETERS - VESICLE DIMENSION AND PERMEABILITY %%

% Starting radius, surface area and volume of the vesicle
%r = 119E-7;        %[cm]
%S = 4*pi*r^2;    %[cm^2]
%V0 = (4/3)*pi*r^3; %[cm^3]

% Permeabilities
Pw = 0.003;      %[cm/s] H2O
P2 = 0.000;    %[cm/s] AH

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
x0(2) = 10E-6;      %[mol/cm^3]
x0(3) = KPH;      %[mol/cm^3]
x0(4) = KP;       %[mol/cm^3]
x0(5) = (10^-pH0)*10^-3; %[mol/cm^3]
x0(6) = 0;     %[mol/cm^3]
x0(7) = 0;      %[mol/cm^3]
%x0(8) = V0;      %[cm^3]

% Scaling of the equilibrium constants
K1 = K1/cs_star;
K2 = K2/cs_star;
% Scaling of the external AH concentration
c6_star = c6_star/cs_star;

%% - FIT OF THE DATA - %%

% I perform a fit of the KCl data

% PARS: 
% 1) K_SV [cm^3/mol]
% 2) Pw [cm/s]
% 3) P2 [cm/s]

START_GUESS = [K_SV Pw P2];

[BEST_PARS,errs,chi2,errmatrix] = fminuit('ODE_AH_fitfun_Fr',START_GUESS,'runmode');


