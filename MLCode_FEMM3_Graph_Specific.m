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

cell=8;

subplot (3,1,1)
heatmap(realVC{cell}(1:12,1:30))
colormap('jet')
grid off

subplot (3,1,3)
[x,y]=meshgrid(0:1:29,0:1:29);
quiver(x,y,realER{(cell)},realEZ{(cell)})
axis([-1 30 -1 12]);

subplot (3,1,2)
heatmap(realmagE{cell}(1:12,1:30))
colormap('jet')
grid off