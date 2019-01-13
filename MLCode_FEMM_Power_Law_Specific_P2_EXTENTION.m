%This program is used for calculation of FEMM Data for Smectite Clay
    %simulation (Surface)with an included power law for conductivity\
    %(NOTE: THIS PROGRAM WILL ONLY ANALYSE Z = 1.25 TO 2.00; USE FILE:P2
    %(REGULAR) FOR Z = 0.25 TO 1.75)
    %This program was started by Aydin Wells 7/25/2018

%Notes:Program assumes there are 4 "z" distances under a specific frequency (8 files)

clc
clear
clf
format long

freq=[0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1,5,10,50,100,500,1000,5000,10000,20000,30000]; %set for frequencies
fix=9
f=freq(fix); %change for specific frequency desired
chart=1;

mtx125a=load((sprintf('z = 1.250000e+00_at_%d_top.txt',f))); %data associated with specific __tropy, permiability and frequency (150 values)
mtx125b=load((sprintf('z = 1.250000e+00_at_%d_bot.txt',f)));
mtx150a=load((sprintf('z = 1.500000e+00_at_%d_top.txt',f)));
mtx150b=load((sprintf('z = 1.500000e+00_at_%d_bot.txt',f)));
mtx175a=load((sprintf('z = 1.750000e+00_at_%d_top.txt',f)));
mtx175b=load((sprintf('z = 1.750000e+00_at_%d_bot.txt',f)));
mtx200a=load((sprintf('z = 2_at_%d_top.txt',f)));
mtx200b=load((sprintf('z = 2_at_%d_bot.txt',f)));

r=mtx125a(:,1); %general length matrix
count=0;

    while count<4 %count used for 4 sections of z-value distances
        
        if count==0
    x1=1.125;
    x2=1.25;
    x3=1.375;
    mtxa=mtx125a;
    mtxb=mtx125b;
        elseif count==1
    x1=1.375;
    x2=1.50;
    x3=1.625;
    mtxa=mtx150a;
    mtxb=mtx150b;
        elseif count==2
    x1=1.625;
    x2=1.75;
    x3=1.875;
    mtxa=mtx175a;
    mtxb=mtx175b;
        elseif count==3
    x1=1.875;
    x2=2.00;
    x3=2.125;
    mtxa=mtx200a;
    mtxb=mtx200b;
        end
        
dif3=(mtxa(:,3)-mtxb(:,3))*(8.854*10^-12); %surface charge calculation (real)
dif4=(mtxa(:,4)-mtxb(:,4))*(8.854*10^-12); %surface charge calculation (imag)
a1=((dif3.*r)./((((r.^2)+(x1.^2)).^(3/2))));
a2=((dif3.*r)./((((r.^2)+(x2.^2)).^(3/2))));
a3=((dif3.*r)./((((r.^2)+(x3.^2)).^(3/2))));
b1=((dif4.*r)./((((r.^2)+(x1.^2)).^(3/2))));
b2=((dif4.*r)./((((r.^2)+(x2.^2)).^(3/2))));
b3=((dif4.*r)./((((r.^2)+(x3.^2)).^(3/2))));
Efield(1,1)=sum((r(2:150,:)-r(1:149,:)).*(a1(2:150,:)+a1(1:149,:))/2);
Efield(2,1)=sum((r(2:150,:)-r(1:149,:)).*(a2(2:150,:)+a2(1:149,:))/2);
Efield(3,1)=sum((r(2:150,:)-r(1:149,:)).*(a3(2:150,:)+a3(1:149,:))/2);
Efield(1,2)=sum((r(2:150,:)-r(1:149,:)).*(b1(2:150,:)+b1(1:149,:))/2);
Efield(2,2)=sum((r(2:150,:)-r(1:149,:)).*(b2(2:150,:)+b2(1:149,:))/2);
Efield(3,2)=sum((r(2:150,:)-r(1:149,:)).*(b3(2:150,:)+b3(1:149,:))/2);
FinalMatrix(1,1)=Efield(1,1)*(x1/(2*8.854*10^-12)); %E-field calculation (real)
FinalMatrix(2,1)=Efield(2,1)*(x2/(2*8.854*10^-12));
FinalMatrix(3,1)=Efield(3,1)*(x3/(2*8.854*10^-12));
FinalMatrix(1,2)=Efield(1,2)*(x1/(2*8.854*10^-12)); %E-field calculation (imag)
FinalMatrix(2,2)=Efield(2,2)*(x2/(2*8.854*10^-12));
FinalMatrix(3,2)=Efield(3,2)*(x3/(2*8.854*10^-12));

        if count==0
    SurfChrgMatrix125=[dif3,dif4];
    SeriesIntegralMatrix125=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix125=Efield;
    FinalMatrix125=FinalMatrix;
        elseif count==1
    SurfChrgMatrix150=[dif3,dif4];
    SeriesIntegralMatrix150=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix150=Efield;
    FinalMatrix150=FinalMatrix;
        elseif count==2
    SurfChrgMatrix175=[dif3,dif4];
    SeriesIntegralMatrix175=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix175=Efield;
    FinalMatrix175=FinalMatrix;
        elseif count==3
    SurfChrgMatrix200=[dif3,dif4];
    SeriesIntegralMatrix200=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix200=Efield;
    FinalMatrix200=FinalMatrix;
        end
        
count=count+1;
    end
    
SurfChrgMatrix=[SurfChrgMatrix125,SurfChrgMatrix150,SurfChrgMatrix175,SurfChrgMatrix200];
    %1-150:1-2 are 125 real-imag. real-imag. ...
    %1-150:3-4 are 150 real-imag. real-imag. ...
    %1-150:5-6 are 175 real-imag. real-imag. ...
    %1-150:7-8 are 200 real-imag. real-imag. ...
SeriesIntegralMatrix=[SeriesIntegralMatrix125,SeriesIntegralMatrix150,SeriesIntegralMatrix175,SeriesIntegralMatrix200];
    %1-150:1-6 are 125 real-imag. real-imag. real-imag. ...
    %1-150:7-12 are 150 real-imag. real-imag. real-imag. ...
    %1-150:13-18 are 175 real-imag. real-imag. real-imag. ...
    %1-150:19-24 are 200 real-imag. real-imag. real-imag. ...
IntegralMatrix=[IntegralMatrix125,IntegralMatrix150,IntegralMatrix175,IntegralMatrix200];
    %1-3:1-2 are 125 real-imag. ...
    %1-6:3-4 are 150 real-imag. ...
    %1-6:5-6 are 175 real-imag. ...
    %1-6:7-8 are 200 real-imag. ...
FinalMatrix=[FinalMatrix125,FinalMatrix150,FinalMatrix175,FinalMatrix200];
    %1-6:1-2 are 125 real-imag. ...
    %1-6:3-4 are 150 real-imag. ...
    %1-6:5-6 are 175 real-imag. ...
    %1-6:7-8 are 200 real-imag. ...
FinalMatrixIm=[FinalMatrix(:,2);FinalMatrix(:,4);FinalMatrix(:,6);FinalMatrix(:,8)];

disp('The resulting data (above) is listed in the following order for each property, respectively:');
disp('    z=1.25');
disp('    z=1.50');
disp('    z=1.75');
disp('    z=2.00');

ElectricalFieldReal=[FinalMatrix(2,1);FinalMatrix(2,3);FinalMatrix(2,5);FinalMatrix(2,7)]
ElectricalFieldImaginary=[FinalMatrix(2,2);FinalMatrix(2,4);FinalMatrix(2,6);FinalMatrix(2,8)]

dx1=diff(1.125:0.125:1.375);
dy1=diff(FinalMatrix(:,1)');
dx2=diff(1.375:0.125:1.625);
dy2=diff(FinalMatrix(:,3)');
dx3=diff(1.625:0.125:1.875);
dy3=diff(FinalMatrix(:,5)');
dx4=diff(1.875:0.125:2.125);
dy4=diff(FinalMatrix(:,7)');
deriv125re=abs(dy1/dx1); %force derivative calculation (real)
deriv150re=abs(dy2/dx2);
deriv175re=abs(dy3/dx3);
deriv200re=abs(dy4/dx4);
DerivativeElectricalFieldReal=[deriv125re;deriv150re;deriv175re;deriv200re]

dx5=diff(1.125:0.125:1.375);
dy5=diff(FinalMatrix(:,2)');
dx6=diff(1.375:0.125:1.625);
dy6=diff(FinalMatrix(:,4)');
dx7=diff(1.625:0.125:1.875);
dy7=diff(FinalMatrix(:,6)');
dx8=diff(1.875:0.125:2.125);
dy8=diff(FinalMatrix(:,8)');
deriv125im=(dy5/dx5); %force derivative calculation (imag)
deriv150im=(dy6/dx6);
deriv175im=(dy7/dx7);
deriv200im=(dy8/dx8);
DerivativeElectricalFieldImaginary=[deriv125im;deriv150im;deriv175im;deriv200im]

if chart == 1

subplot(2,2,1); %plot Surface Charge, z=1.25
hold on
plot(r,SurfChrgMatrix(:,1),'-b')
plot(r,SurfChrgMatrix(:,2),'-b')
title('"rho" vs. r (z=1.25)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,2); %plot Surface Charge, z=1.50
hold on
plot(r,SurfChrgMatrix(:,3),'-r')
plot(r,SurfChrgMatrix(:,4),'-r')
title('"rho" vs. r (z=1.50)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,3); %plot Surface Charge, z=1.75
hold on
plot(r,SurfChrgMatrix(:,5),'-c')
plot(r,SurfChrgMatrix(:,6),'-c')
title('"rho" vs. r (z=1.75)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,4); %plot Surface Charge, z=2.00
hold on
plot(r,SurfChrgMatrix(:,7),'-g')
plot(r,SurfChrgMatrix(:,8),'-g')
title('"rho" vs. r (z=2.00)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

saveas(gcf,(sprintf('SC_freq=%d.jpg',f)));

elseif chart == 2

subplot(2,2,1) %plot E-field on all z-value distances
hold on
plot((1.125:0.125:1.375),FinalMatrix(:,1),'-bo')
plot((1.375:0.125:1.625),FinalMatrix(:,3),'-bo')
plot((1.625:0.125:1.875),FinalMatrix(:,5),'-bo')
plot((1.875:0.125:2.125),FinalMatrix(:,7),'-bo')
plot([(1.125:0.125:1.375)';(1.375:0.125:1.625)';(1.625:0.125:1.875)';(1.875:0.125:2.125)'],FinalMatrixIm,'-ro')
title('E(z) (all)')
xlabel('z')
ylabel('E')
grid on
hold off    

subplot(2,2,2) %plot E-field on main z-value distances
hold on
f=fit((1.25:0.25:2)',ElectricalFieldReal,'exp2');
plot(f,(1.25:0.25:2),ElectricalFieldReal,'ob');
plot((1.25:0.25:2),ElectricalFieldImaginary,'-ro')
f2=fit((1.25:0.25:2)',abs(ElectricalFieldReal),'exp2');
plot(f2,(1.25:0.25:2),abs(ElectricalFieldReal),'oc')
title('E(z) (main)')
xlabel('z')
ylabel('E')
grid on
legend off
hold off
%{
subplot(4,3,1); %plot Series Integral, z=1.25
hold on
plot(r,SeriesIntegralMatrix125(:,3),'-b')
plot(r,SeriesIntegralMatrix125(:,4),'-b')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=1.25)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,2); %plot Series Integral, z=1.50
hold on
plot(r,SeriesIntegralMatrix150(:,3),'-r')
plot(r,SeriesIntegralMatrix150(:,4),'-r')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=1.50)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,7); %plot Series Integral, z=1.75
hold on
plot(r,SeriesIntegralMatrix175(:,3),'-c')
plot(r,SeriesIntegralMatrix175(:,4),'-c')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=1.75)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,8); %plot Series Integral, z=2.00
hold on
plot(r,SeriesIntegralMatrix200(:,3),'-g')
plot(r,SeriesIntegralMatrix200(:,4),'-g')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=2.00)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off
%}
subplot(2,2,3); %plot Force Derivative both real and imag.
hold on
f3=fit((1.25:0.25:2)',DerivativeElectricalFieldReal,'exp2');
f4=fit((1.25:0.25:2)',DerivativeElectricalFieldImaginary,'exp2');
plot(f3,(1.25:0.25:2),DerivativeElectricalFieldReal,'oc')
plot(f4,(1.25:0.25:2),DerivativeElectricalFieldImaginary,'or')
title('|dE/dz| vs. z')
xlabel('z')
ylabel('|dE/dz|')
grid on
legend off
hold off

subplot(2,2,4); %plot Force Derivative Ratio
hold on
f5=fit((1.25:0.25:2)',(DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal),'exp2');
plot(f5,(1.25:0.25:2)',(DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal),'ob')
title('|(dIm[E]/dz)/(dRe[E]/dz)| vs. z')
xlabel('z')
ylabel('|(dIm[E]/dz)/(dRe[E]/dz)|')
grid on
legend off
hold off

f=freq(fix);
saveas(gcf,(sprintf('E_F_etc_freq=%d.jpg',f)));
end

PhiValues=((DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal)*(180/pi)) %calculation for phase shift