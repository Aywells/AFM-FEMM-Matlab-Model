%This program creates, analises and files FEMM simulation data using a power law (Surface)
%This program was started by Aydin Wells 7/2/2018
%Edit appropriately for what data to collect and what figures to scan for

clc
clear
clf
format long

f=[0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1,5,10,50,100,500,1000,5000,10000,20000,30000]; %set for frequencies
cond=[4.35967E-09,5.36739E-09,6.06164E-09,7.06553E-09,8.69869E-09,1.07093E-08,1.20946E-08,1.40976E-08,1.73562E-08,2.81284E-08,3.46301E-08,5.61235E-08,6.90962E-08,1.11981E-07,1.37865E-07,2.23432E-07,2.75077E-07,3.38659E-07,3.82464E-07]; %set for conductivities

addpath('C:\femm42\mfiles');
savepath;

    freqset=1;
while (freqset<=19) %frequency and conductivity loop (using power law)
    
    freq=f(freqset);
    sig=cond(freqset);

    movement=4; 
while (movement<=7) %z-distance movement loop
    openfemm
    opendocument(sprintf('3.17.17_1_%d.fec',movement)); %femm prototype loaded
    
    ci_probdef('meters','axi',freq,0.00000001,0,30); %set femm settings
    ci_modifymaterial('Clay',1,sig)%r-conductivity
    ci_modifymaterial('Clay',2,sig)%z-conductivity (change to 0 if anisotropic)
    
    if movement==4 %movement to z value loop
           dist=1.25;
        elseif movement==5
           dist=1.50; 
        elseif movement==6
           dist=1.75; 
        elseif movement==7
           dist=2.00; 
    end
    
    ci_analyze; %analysis solution with settings
    ci_loadsolution;
    co_zoom(-.1,-.1,1.1,1);
    co_seteditmode('contour');

         co_selectpoint(0.03,0.781); %scan of the portion above surface
         co_selectpoint(1.63,0.781);
         co_makeplot(5,150,(sprintf('z = %d_at_%d_top.txt',dist,freq)),0);
         co_seteditmode('point');
         co_seteditmode('contour');
         
         co_selectpoint(0.03,0.715); %scan of the portion below surface
         co_selectpoint(1.63,0.715);
         co_makeplot(5,150,(sprintf('z = %d_at_%d_bot.txt',dist,freq)),0);
         co_seteditmode('point');
         co_seteditmode('contour');
         
    movement=movement+1;
    closefemm
end
freqset=freqset+1;
end
