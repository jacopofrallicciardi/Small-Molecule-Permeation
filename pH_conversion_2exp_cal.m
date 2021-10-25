clear all;

[filename, pathname] = uigetfile('*.txt*', 'MultiSelect','off');

INPUT = importdata(filename,'\t');

t = (1:length(INPUT.data))*6;
f390 = INPUT.data(:,4);
f470 = INPUT.data(:,5);

r = f390./f470;

% pH-conversion with the bi-exponential calibration curve (see 23 Aug 2017)
a = 5.33;
b = 0.1507;
c = -5.195; 
d = -5.109;

pH = a*exp(b*r) + c*exp(d*r);

%figure
plot(t,pH,'.')

%% - DATA OUTPUT - %%

% I save the data = [t,pH,r]
data = [t',pH,r];
output_name = strrep(filename , '.txt', '_pH+RATIO.dat');
dlmwrite(output_name,data,'delimiter','\t')

