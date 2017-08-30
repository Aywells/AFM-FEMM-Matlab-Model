%This program is used for calculation of FEMM Data
%This program was made by Aydin Wells 8/10/17

    %Notes:Program assumes there are 4 "z" distances (8 files)
clc
clear
format long

mtx025a=load('z = 0.25 above.txt');
mtx025b=load('z = 0.25 below.txt');
mtx050a=load('z = 0.50 above.txt');
mtx050b=load('z = 0.50 below.txt');
mtx075a=load('z = 0.75 above.txt');
mtx075b=load('z = 0.75 below.txt');
mtx100a=load('z = 1.00 above.txt');
mtx100b=load('z = 1.00 below.txt');
    %Put files to load into folder with program to work

r=mtx025a(:,1);
count=0;

    while count<4
        
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
        
dif3=(mtxa(:,3)-mtxb(:,3))*(8.854*10^-12);
dif4=(mtxa(:,4)-mtxb(:,4))*(8.854*10^-12);
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
FinalMatrix(1,1)=Efield(1,1)*(x1/(2*8.854*10^-12));
FinalMatrix(2,1)=Efield(2,1)*(x2/(2*8.854*10^-12));
FinalMatrix(3,1)=Efield(3,1)*(x3/(2*8.854*10^-12));
FinalMatrix(1,2)=Efield(1,2)*(x1/(2*8.854*10^-12));
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
    %1-3:1-2 are 025 real-imag.
    %1-6:3-4 are 050 real-imag.
    %1-6:5-6 are 075 real-imag.
    %1-6:7-8 are 100 real-imag.
FinalMatrix=[FinalMatrix025,FinalMatrix050,FinalMatrix075,FinalMatrix100];
    %1-6:1-2 are 025 real-imag.
    %1-6:3-4 are 050 real-imag.
    %1-6:5-6 are 075 real-imag.
    %1-6:7-8 are 100 real-imag.
    
FinalMatrixIm=[FinalMatrix(:,2);FinalMatrix(:,4);FinalMatrix(:,6);FinalMatrix(:,8)];

subplot(4,3,3)
hold on
plot((0.125:0.125:0.375),FinalMatrix(:,1),'-b.')
plot((0.375:0.125:0.625),FinalMatrix(:,3),'-b.')
plot((0.625:0.125:0.875),FinalMatrix(:,5),'-b.')
plot((0.875:0.125:1.125),FinalMatrix(:,7),'-b.')
plot([(0.125:0.125:0.375)';(0.375:0.125:0.625)';(0.625:0.125:0.875)';(0.875:0.125:1.125)'],FinalMatrixIm,'-r.')
title('E(z) (all)')
xlabel('z')
ylabel('E')
grid on
hold off

FinalMatrixRemain=[FinalMatrix(2,1);FinalMatrix(2,3);FinalMatrix(2,5);FinalMatrix(2,7)];
FinalMatrixImmain=[FinalMatrix(2,2);FinalMatrix(2,4);FinalMatrix(2,6);FinalMatrix(2,8)];

subplot(4,3,6)
hold on
f=fit((0.25:0.25:1)',FinalMatrixRemain,'exp1');
plot(f,(0.25:0.25:1),FinalMatrixRemain,'.b');
plot((0.25:0.25:1),FinalMatrixImmain,'-r.')
f2=fit((0.25:0.25:1)',abs(FinalMatrixRemain),'exp1');
plot(f2,(0.25:0.25:1),abs(FinalMatrixRemain),'.c')
title('E(z) (main)')
xlabel('z')
ylabel('E')
grid on
hold off

%{
subplot(4,3,9);
hold on
f=fit((0.25:0.25:1)',dydx(FinalMatrixRemain),'exp1');
plot(f,(0.25:0.25:1),dydx(FinalMatrixRemain),'.b');
plot((0.25:0.25:1),dydx(FinalMatrixImmain),'-r.')
title('dE/dz vs. z')
xlabel('z')
ylabel('dE/dz')
grid on
hold off
%}

subplot(4,3,1);
hold on
plot(r,SeriesIntegralMatrix025(:,3),'-b')
plot(r,SeriesIntegralMatrix025(:,4),'-b')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.25)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,2);
hold on
plot(r,SeriesIntegralMatrix050(:,3),'-r')
plot(r,SeriesIntegralMatrix050(:,4),'-r')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.50)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,7);
hold on
plot(r,SeriesIntegralMatrix075(:,3),'-c')
plot(r,SeriesIntegralMatrix075(:,4),'-c')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=0.75)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,8);
hold on
plot(r,SeriesIntegralMatrix100(:,3),'-g')
plot(r,SeriesIntegralMatrix100(:,4),'-g')
title('("rho"*r)/(r^2 + z^2)^3/2 vs. r (z=1.00)')
xlabel('r')
ylabel('("rho"*r)/(r^2 + z^2)^3/2')
grid on
hold off

subplot(4,3,4);
hold on
plot(r,SurfChrgMatrix(:,1),'-b')
plot(r,SurfChrgMatrix(:,2),'-b')
title('"rho" vs. r (z=0.25)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(4,3,5);
hold on
plot(r,SurfChrgMatrix(:,3),'-r')
plot(r,SurfChrgMatrix(:,4),'-r')
title('"rho" vs. r (z=0.50)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(4,3,10);
hold on
plot(r,SurfChrgMatrix(:,5),'-c')
plot(r,SurfChrgMatrix(:,6),'-c')
title('"rho" vs. r (z=0.75)')
xlabel('r')
ylabel('"rho"')
grid on
hold off

subplot(4,3,11);
hold on
plot(r,SurfChrgMatrix(:,7),'-g')
plot(r,SurfChrgMatrix(:,8),'-g')
title('"rho" vs. r (z=1.00)')
xlabel('r')
ylabel('"rho"')
grid on
hold off