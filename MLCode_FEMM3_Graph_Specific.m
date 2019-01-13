%This program isolates the plots desired based on z-value and frequency;
%Use MLCode_FEMM3.m jointly to get data for using this code
%This program was made by Aydin Wells 5/18/18

%Note: In each cell there are 48 matrices every combination of 4 is a freq.

%1-4: 0.12
%5-8: 0.25
%9-12: 0.5
%13-16: 1
%17-20: 2
%21-24: 5
%25-28: 10
%29-32: 20
%33-36: 30
%37-40: 60
%41-44: 120
%45-48: 600
%PhiValues=((FImaginary./FReal)*(180/pi)) %calculation for phase shift

cell=43;

VCreal_tot_squ=sum(sum(realVC{cell}(1:11,1:29)));
VCimag_tot_squ=sum(sum(imagVC{cell}(1:11,1:29))); 
Ereal_tot_squ=sum(sum(realmagE{cell}(1:11,1:29)));
Eimag_tot_squ=sum(sum(imagmagE{cell}(1:11,1:29)));
Freal_tot_squ=sum(sum(realF{cell}(1:11,1:29)));
Fimag_tot_squ=sum(sum(imagF{cell}(1:11,1:29)));

VCreal_tot_vol=pi*VCreal_tot_squ*(1.74^2)*0.66;
VCimag_tot_vol=pi*VCimag_tot_squ*(1.74^2)*0.66;
Ereal_tot_vol=pi*Ereal_tot_squ*(1.74^2)*0.66;
Eimag_tot_vol=pi*Eimag_tot_squ*(1.74^2)*0.66;
Freal_tot_vol=pi*Freal_tot_squ*(1.74^2)*0.66;
Fimag_tot_vol=pi*Fimag_tot_squ*(1.74^2)*0.66;

P=(Fimag_tot_vol/Freal_tot_vol);

subplot (3,1,1)
heatmap(realVC{cell}(1:11,1:29))
colormap('jet')
title('RealVC')
grid off

subplot (3,1,2)
heatmap(realmagE{cell}(1:11,1:29))
colormap('jet')
title('RealMagE')
grid off

subplot (3,1,3)
heatmap(realF{cell}(1:11,1:29))
colormap('jet')
title('RealF')
grid off

disp(VCreal_tot_vol)
disp(VCimag_tot_vol)
disp(Ereal_tot_vol)
disp(Eimag_tot_vol)
disp(Freal_tot_vol)
disp(Fimag_tot_vol)
disp(P)
