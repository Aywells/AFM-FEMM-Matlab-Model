%This program is used for calculation of FEMM Data for Smectite Clay
    %simulation (Surface)with an included power law for conductivity
    %This program was started by Aydin Wells 7/2/2018

%Notes:Program assumes there are 4 "z" distances under a specific frequency (8 files)

clc
clear
clf
format long

freq=[0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1,5,10,50,100,500,1000,5000,10000,20000,30000]; %set for frequencies
fix=9
f=freq(fix); %change for specific frequency desired
chart=1;

mtx025a=load((sprintf('z = 2.500000e-01_at_%d_top.txt',f))); %data associated with specific __tropy, permiability and frequency (150 values)
mtx025b=load((sprintf('z = 2.500000e-01_at_%d_bot.txt',f)));
mtx050a=load((sprintf('z = 5.000000e-01_at_%d_top.txt',f)));
mtx050b=load((sprintf('z = 5.000000e-01_at_%d_bot.txt',f)));
mtx075a=load((sprintf('z = 7.500000e-01_at_%d_top.txt',f)));
mtx075b=load((sprintf('z = 7.500000e-01_at_%d_bot.txt',f)));
mtx100a=load((sprintf('z = 1_at_%d_top.txt',f)));
mtx100b=load((sprintf('z = 1_at_%d_bot.txt',f)));

r=mtx025a(:,1); %general length matrix
count=0;

    while count<4 %count used for 4 sections of z-value distances
        
        if count==0
    x1=0.125;
    x2=0.25;
    x3=0.375;
    mtxa=mtx025a;
    mtxb=mtx025b;
        elseif count==1
    x1=0.375;
    x2=0.50;
    x3=0.625;
    mtxa=mtx050a;
    mtxb=mtx050b;
        elseif count==2
    x1=0.625;
    x2=0.75;
    x3=0.875;
    mtxa=mtx075a;
    mtxb=mtx075b;
        elseif count==3
    x1=0.875;
    x2=1.00;
    x3=1.125;
    mtxa=mtx100a;
    mtxb=mtx100b;
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
    SurfChrgMatrix025=[dif3,dif4];
    SeriesIntegralMatrix025=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix025=Efield;
    FinalMatrix025=FinalMatrix;
        elseif count==1
    SurfChrgMatrix050=[dif3,dif4];
    SeriesIntegralMatrix050=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix050=Efield;
    FinalMatrix050=FinalMatrix;
        elseif count==2
    SurfChrgMatrix075=[dif3,dif4];
    SeriesIntegralMatrix075=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix075=Efield;
    FinalMatrix075=FinalMatrix;
        elseif count==3
    SurfChrgMatrix100=[dif3,dif4];
    SeriesIntegralMatrix100=[a1,b1,a2,b2,a3,b3];
    IntegralMatrix100=Efield;
    FinalMatrix100=FinalMatrix;
        end
        
count=count+1;
    end
    
SurfChrgMatrix=[SurfChrgMatrix025,SurfChrgMatrix050,SurfChrgMatrix075,SurfChrgMatrix100];
    %1-150:1-2 are 025 real-imag. real-imag. ...
    %1-150:3-4 are 050 real-imag. real-imag. ...
    %1-150:5-6 are 075 real-imag. real-imag. ...
    %1-150:7-8 are 100 real-imag. real-imag. ...
SeriesIntegralMatrix=[SeriesIntegralMatrix025,SeriesIntegralMatrix050,SeriesIntegralMatrix075,SeriesIntegralMatrix100];
    %1-150:1-6 are 025 real-imag. real-imag. real-imag. ...
    %1-150:7-12 are 050 real-imag. real-imag. real-imag. ...
    %1-150:13-18 are 075 real-imag. real-imag. real-imag. ...
    %1-150:19-24 are 100 real-imag. real-imag. real-imag. ...
IntegralMatrix=[IntegralMatrix025,IntegralMatrix050,IntegralMatrix075,IntegralMatrix100];
    %1-3:1-2 are 025 real-imag. ...
    %1-6:3-4 are 050 real-imag. ...
    %1-6:5-6 are 075 real-imag. ...
    %1-6:7-8 are 100 real-imag. ...
FinalMatrix=[FinalMatrix025,FinalMatrix050,FinalMatrix075,FinalMatrix100];
    %1-6:1-2 are 025 real-imag. ...
    %1-6:3-4 are 050 real-imag. ...
    %1-6:5-6 are 075 real-imag. ...
    %1-6:7-8 are 100 real-imag. ...
FinalMatrixIm=[FinalMatrix(:,2);FinalMatrix(:,4);FinalMatrix(:,6);FinalMatrix(:,8)];

disp('The resulting data (above) is listed in the following order for each property, respectively:');
disp('    z=0.25');
disp('    z=0.50');
disp('    z=0.75');
disp('    z=1.00');

ElectricalFieldReal=[FinalMatrix(2,1);FinalMatrix(2,3);FinalMatrix(2,5);FinalMatrix(2,7)]
ElectricalFieldImaginary=[FinalMatrix(2,2);FinalMatrix(2,4);FinalMatrix(2,6);FinalMatrix(2,8)]

dx1=diff(0.125:0.125:0.375);
dy1=diff(FinalMatrix(:,1)');
dx2=diff(0.375:0.125:0.625);
dy2=diff(FinalMatrix(:,3)');
dx3=diff(0.625:0.125:0.875);
dy3=diff(FinalMatrix(:,5)');
dx4=diff(0.875:0.125:1.125);
dy4=diff(FinalMatrix(:,7)');
deriv125re=abs(dy1/dx1); %force derivative calculation (real)
deriv150re=abs(dy2/dx2);
deriv175re=abs(dy3/dx3);
deriv100re=abs(dy4/dx4);
DerivativeElectricalFieldReal=[deriv125re;deriv150re;deriv175re;deriv100re]

dx5=diff(0.125:0.125:0.375);
dy5=diff(FinalMatrix(:,2)');
dx6=diff(0.375:0.125:0.625);
dy6=diff(FinalMatrix(:,4)');
dx7=diff(0.625:0.125:0.875);
dy7=diff(FinalMatrix(:,6)');
dx8=diff(0.875:0.125:1.125);
dy8=diff(FinalMatrix(:,8)');
deriv125im=(dy5/dx5); %force derivative calculation (imag)
deriv150im=(dy6/dx6);
deriv175im=(dy7/dx7);
deriv100im=(dy8/dx8);
DerivativeElectricalFieldImaginary=[deriv125im;deriv150im;deriv175im;deriv100im]

if chart == 1

subplot(2,2,1) %plot E-field on all z-value distances
hold on
plot((0.125:0.125:0.375),FinalMatrix(:,1),'-bo')
plot((0.375:0.125:0.625),FinalMatrix(:,3),'-bo')
plot((0.625:0.125:0.875),FinalMatrix(:,5),'-bo')
plot((0.875:0.125:1.125),FinalMatrix(:,7),'-bo')
plot([(0.125:0.125:0.375)';(0.375:0.125:0.625)';(0.625:0.125:0.875)';(0.875:0.125:1.125)'],FinalMatrixIm,'-ro')
title('E(z) (all)')
xlabel('z')
ylabel('E')
grid on
hold off

subplot(2,2,2) %plot E-field on main z-value distances
hold on
f=fit((0.25:0.25:1)',ElectricalFieldReal,'exp2');
plot(f,(0.25:0.25:1),ElectricalFieldReal,'ob');
plot((0.25:0.25:1),ElectricalFieldImaginary,'-ro')
f2=fit((0.25:0.25:1)',abs(ElectricalFieldReal),'exp2');
plot(f2,(0.25:0.25:1),abs(ElectricalFieldReal),'oc')
title('E(z) (main)')
xlabel('z')
ylabel('E')
grid on
legend off
hold off
%{
subplot(4,3,1); %plot Series Integral, z=0.25
hold on
plot(r,SeriesIntegralMatrix025(:,3),'-b')
plot(r,SeriesIntegralMatrix025(:,4),'-b')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.25)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,2); %plot Series Integral, z=0.50
hold on
plot(r,SeriesIntegralMatrix050(:,3),'-r')
plot(r,SeriesIntegralMatrix050(:,4),'-r')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.50)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,7); %plot Series Integral, z=0.75
hold on
plot(r,SeriesIntegralMatrix075(:,3),'-c')
plot(r,SeriesIntegralMatrix075(:,4),'-c')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.75)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,8); %plot Series Integral, z=1.00
hold on
plot(r,SeriesIntegralMatrix100(:,3),'-g')
plot(r,SeriesIntegralMatrix100(:,4),'-g')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=1.00)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off
%}
subplot(2,2,3); %plot Force Derivative both real and imag.
hold on
f3=fit((0.25:0.25:1)',DerivativeElectricalFieldReal,'exp2');
f4=fit((0.25:0.25:1)',DerivativeElectricalFieldImaginary,'exp2');
plot(f3,(0.25:0.25:1),DerivativeElectricalFieldReal,'oc')
plot(f4,(0.25:0.25:1),DerivativeElectricalFieldImaginary,'or')
title('|dE/dz| vs. z')
xlabel('z')
ylabel('|dE/dz|')
grid on
legend off
hold off

subplot(2,2,4); %plot Force Derivative Ratio
hold on
f5=fit((0.25:0.25:1)',(DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal),'exp2');
plot(f5,(0.25:0.25:1)',(DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal),'ob')
title('|(dIm[E]/dz)/(dRe[E]/dz)| vs. z')
xlabel('z')
ylabel('|(dIm[E]/dz)/(dRe[E]/dz)|')
grid on
legend off
hold off

f=freq(fix);
saveas(gcf,(sprintf('E_F_etc_freq=%d.jpg',f)));

elseif chart == 2

subplot(2,2,1); %plot Surface Charge, z=0.25
hold on
plot(r,SurfChrgMatrix(:,1),'-b')
plot(r,SurfChrgMatrix(:,2),'-b')
title('"rho" vs. r (z=0.25)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,2); %plot Surface Charge, z=0.50
hold on
plot(r,SurfChrgMatrix(:,3),'-r')
plot(r,SurfChrgMatrix(:,4),'-r')
title('"rho" vs. r (z=0.50)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,3); %plot Surface Charge, z=0.75
hold on
plot(r,SurfChrgMatrix(:,5),'-c')
plot(r,SurfChrgMatrix(:,6),'-c')
title('"rho" vs. r (z=0.75)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(2,2,4); %plot Surface Charge, z=1.00
hold on
plot(r,SurfChrgMatrix(:,7),'-g')
plot(r,SurfChrgMatrix(:,8),'-g')
title('"rho" vs. r (z=1.00)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

saveas(gcf,(sprintf('SC_freq=%d.jpg',f)));

end

PhiValues=((DerivativeElectricalFieldImaginary./DerivativeElectricalFieldReal)*(180/pi)) %calculation for phase shift
