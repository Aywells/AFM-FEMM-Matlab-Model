%This program creates, analises and files FEMM simulation data (Volumetric)
%This program was started by Aydin Wells 1/10/18

clc
clear
clf

addpath('C:\femm42\mfiles');
savepath;

    freqset=1;
while (freqset<=12) %frequency loop
    if freqset==1
        freq=0.12;
    elseif freqset==2
        freq=0.25;
    elseif freqset==3
        freq=0.5;
    elseif freqset==4
        freq=1;
    elseif freqset==5
        freq=2;
    elseif freqset==6
        freq=5;
    elseif freqset==7
        freq=10;
    elseif freqset==8
        freq=20;
    elseif freqset==9
        freq=30;
    elseif freqset==10
        freq=60;
    elseif freqset==11
        freq=120;
    elseif freqset==12
        freq=600;
    end

    movement=0; 
while (movement<=3) %z-distance movement loop
    openfemm
    opendocument(sprintf('3.17.17_2_%d.fec',movement)); %femm prototype loaded
    
    ci_probdef('meters','axi',freq,0.00000001,0,30); %set femm settings
    ci_modifymaterial('Clay',1,1e-8)%r-conductivity
    ci_modifymaterial('Clay',2,0)%z-conductivityn 
    
    if movement==0 %movement to z value loop
           dist=0.25;
        elseif movement==1
           dist=0.50; 
        elseif movement==2
           dist=0.75; 
        elseif movement==3
           dist=1.00; 
    end
    
    ci_analyze; %analysis solution with settings
    ci_loadsolution;
    co_zoom(-.1,-.1,1.1,1);
    co_seteditmode('contour');
    hori1=0;
    vert1=0;
        
        while (hori1<=12 && hori2<=12) %horizontal volumetric scan loop
         co_selectpoint(0.01,0.06*hori1);
         co_selectpoint(1.91,0.06*hori1);
         co_makeplot(5,320,(sprintf('z = %d_hori%d_at_%d.txt',dist,hori1,freq)),1);
         co_seteditmode('point');
         co_seteditmode('contour');
         hori1=hori1+1;
        end
    
        while (vert1<=32) %vertical volumetric scan loop
         co_selectpoint(0.06*vert1,0.71);
         co_selectpoint(0.06*vert1,0.01);
         co_makeplot(5,120,(sprintf('z = %d_vert%d_at_%d.txt',dist,vert1,freq)),1);
         co_seteditmode('point');
         co_seteditmode('contour');
         vert1=vert1+1;
        end
    
    movement=movement+1;
    closefemm
end
freqset=freqset+1;
end