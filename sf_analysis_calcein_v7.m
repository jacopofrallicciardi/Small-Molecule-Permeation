%% - HISTORY - %%

% Written by Matteo Gabba %

% The script load and process the calcein fluorescence kinetics obtained
% from stopped flow experiments on liposomes. 
% The output file can be feeded to the fitting routine.

% INPUT: ascii files obtained from the stopped flow apparatus after
%       elimination of the headers and the absorption data (this can be
%       done with a text editor). The format of the input files is the
%       following: (t,f1,f2,...,fN).
% OUTPUT: the output is a *.dat file with the fluorescence intensity
%       F(t)/F(0) normalized to t(0) = 2 ms that is the dead time of the machine. 

%% - DATA INPUT - %%

clear all;

% Load the data file
[filename, pathname] = uigetfile('*.txt*', 'MultiSelect','off');
%INPUT = importdata(filename,'\t');
INPUT = importdata([pathname,filename],'\t'); %HJK edit to include pathname, then the raw data can be in a subfolder.
[r,c] = size(INPUT);

% Calculate the mean and standard deviation from multiple experiments
t = INPUT(:,1);   
f = INPUT(:,2:c);
F =  mean(f,2);
dF = std(f,0,2);

%% - DISCARD POINTS - %%

% I discard the points lying outside the time resolution range of the
% stopped flow apparatus (t<2ms). Indeed, the dead time of the mixing
% device is 2ms with the used mixing chamber.

k = t>0.002;
ts = nonzeros(t.*k);
F = nonzeros(F.*k);
dF = nonzeros(dF.*k);

%% - DATA MODIFICATION - %%

%I calculate F(0) as the mean value of the first N-points
N = 20;
F0 = mean(F(1:N));
dF0 = std(F(1:N));

% Normalize the data to t=0. 
Fn = F/F0;
dFn = Fn.* sqrt((dF./F).^2 + (dF0/F0).^2);

% Plot the data
figure
semilogx(t,f(:,1)/mean(f(1:N,1)),t,f(:,2)/mean(f(1:N,2)),t,f(:,3)/mean(f(1:N,3)))
hold on
semilogx(ts,F/mean(F(1:N)),'-k','LineWidth',2)
errorbar(ts,Fn,dFn)

%% - DATA OUTPUT - %%

% I save the normalized data
data = [ts,Fn,dFn];
output_name = strrep(filename , '.txt', '_Fn_t.dat');
dlmwrite(output_name,data,'delimiter','\t')

