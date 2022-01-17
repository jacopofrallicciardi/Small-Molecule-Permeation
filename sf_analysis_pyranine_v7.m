%% - HISTORY - %%

% Written by Matteo Gabba %

% The script load and process the pyranine fluorescence kinetics obtained
% from stopped flow experiments on liposomes. 
% The output file can be feeded to the fitting routine.

% INPUT: ascii files obtained from the stopped flow apparatus after
%       elimination of the headers and the absorption data (this can be
%       done with a text editor). The format of the input files is the
%       following: (t,f1,f2,...,fN). I need two files for each excitation
%       wavelength (405 and 453 nm), that are measured with either empty
%       liposomes (blank) or pyranine filled liposomes.
% OUTPUT: the output is a *.dat file with the pH kinetics.

% The pH values are calculated using a measured calibration curve.


%% - 405ex - %%

clear all;

% Load the data and the blank correction for the 405ex. The blank is
% measured with empty vesicles at the same concentration of the sample and
% correct for static scattering
[filename_405, pathname_405] = uigetfile('*.txt*', '405nm', 'MultiSelect','on');
[filename_bg405, pathname_bg405] = uigetfile('*.txt*', 'BG-405nm', 'MultiSelect','off');

% INPUT = importdata(filename_405,'\t');
INPUT = importdata([pathname_405,filename_405],'\t');%HJK edit

% bg_405 = importdata(filename_bg405,'\t');
bg_405 = importdata([pathname_bg405,filename_bg405],'\t');%HJK edit

% Calculate the mean bg and the std from multiple measurements
[r,c] = size(bg_405);
bgm_405 = mean(bg_405(:,2:c),2);
dbgm_405 = std(bg_405(:,2:c),0,2);
    
% Calculate the mean fluorescence signal and the std from multiple measurements
t = INPUT(:,1);      
[r,c] = size(INPUT);
f_405 = mean(INPUT(:,2:c),2);
df_405 = std(INPUT(:,2:c),0,2);

% Substract the bg from the fluorescence signal and propagate the errors
F_405 =  f_405 - bgm_405;
dF_405 = F_405.* sqrt((df_405./f_405).^2 + (dbgm_405./bgm_405).^2);


%% - 453ex - %%

% Load the data and the blank correction for the 405ex. The blank is
% measured with empty vesicles at the same concentration of the sample and
% correct for static scattering
[filename_453, pathname_453] = uigetfile('*.txt*', '453nm', 'MultiSelect','on');
[filename_bg453, pathname_bg453] = uigetfile('*.txt*', 'BG-453nm', 'MultiSelect','off');

% INPUT = importdata(filename_453,'\t');
% bg_453 = importdata(filename_bg453,'\t');

INPUT = importdata([pathname_453,filename_453],'\t'); %HJK edit
bg_453 = importdata([pathname_bg453,filename_bg453],'\t'); %HJK edit

% Calculate the mean bg and the std from multiple measurements
[r,c] = size(bg_453);
bgm_453 = mean(bg_453(:,2:c),2);
dbgm_453 = std(bg_453(:,2:c),0,2);

% Calculate the mean fluorescence signal and the std from multiple measurements    
t = INPUT(:,1);     
[r,c] = size(INPUT);
f_453 = mean(INPUT(:,2:c),2);
df_453 = std(INPUT(:,2:c),0,2);

% Substract the bg from the fluorescence signal and propagate the errors
F_453 =  f_453 - bgm_453;
dF_453 = F_453.* sqrt((df_453./f_453).^2 + (dbgm_453./bgm_453).^2);

% Plot the bg-corrected signals for both 405ex and 453ex
figure
semilogx(t,F_405)
hold on
semilogx(t,F_453)
errorbar(t,F_405,dF_405)
errorbar(t,F_453,dF_453)

%% - RATIO - %%

% Calculate the 453/405-ratio and propagate the errors
r = F_453./F_405;
dr = r.* sqrt( (dF_453./F_453).^2 + (dF_405./F_405).^2 );

% Plot the ratio
figure
semilogx(t,r)
hold on
errorbar(t,r,dr)

% Plot the ratio normalized to t=0
figure 
semilogx(t,r./mean(r(1:5)))

%% - CONVERSION TO pH: 2-EXP FUNCTION - %

% The conversion is based on an empirical conversion function
% experimentally determined with pyranine in KPi (100mM) buffer solutions
% at different pHs.
% This conversion is valid only in the pH interval [6, 7.5]

% Coefficients from calibration
a = 6.633;
b = 0.1152;
c = -1.009;
d = -9.241;

% Conversion function and error propagation. I assume that the coefficients
% have negligible relative errors.
pH = a*exp(b*r) + c*exp(d*r);
dpH = pH .* sqrt( a*exp(2*b*r) .* (a*b*dr).^2 + c*exp(2*d*r) .* (c*d*dr).^2  );

% Plot the pH curve
figure
semilogx(t,pH)
hold on
errorbar(t,pH,dpH)


%% - DELETE POINTS - %% 

% I delete points outside the time resolution range of the stopped flow
% machine.

k = t>0.002;

t = nonzeros(t.*k);
pH = nonzeros(pH.*k);
dpH = nonzeros(dpH.*k);
r = nonzeros(r.*k);
dr = nonzeros(dr.*k);


%% - DATA OUTPUT - %%

% I save the data = [t,pH,dpH,r,dr]
data = [t,pH,abs(dpH),r,dr];
output_name = strrep(filename_405 , '.txt', '_pH+RATIO.dat');
dlmwrite(output_name,data,'delimiter','\t')
