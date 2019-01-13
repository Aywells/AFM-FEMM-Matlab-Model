%With the help of MLCode_FEMM2.m, this program calculates FEMM Data for Smectite Clay simulation (Volumetric)
%This program was made by Aydin Wells 3/10/18

clc
clear
clf
format long

z=0; %settings for starting arrays
level=0;
level2=0;
freq=0;
index=1;
column=3;

nanocharge=zeros(31,31); %size of full particle figure for volumetric charge density
nanoelectricR=zeros(31,31); %size of full particle figure for r-direction E-field
nanoelectricZ=zeros(31,31); %size of full particle figure for z-direction E-field
nanoforce=zeros(31,31); %size of full particle figure for F-field
nanoelectricMAG=zeros(31,31); %size of full particle figure for total E-field magnetud

realVC=cell(1,48); %cell for real volumetric charge values
imagVC=cell(1,48); %cell for imaginary volumetric charge values

realER=cell(1,48); %cell for real r-direction E-field values
imagER=cell(1,48); %cell for imag r-direction E-field values
realEZ=cell(1,48); %cell for real z-direction E-field values
imagEZ=cell(1,48); %cell for imag z-direction E-field values

realmagE=cell(1,48); %cell for magnetudes of total real Electric fields
imagmagE=cell(1,48); %cell for magnetudes of total imaginary Electric fields

realF=cell(1,48); %cell for real F-field values
imagF=cell(1,48); %cell for imag F-field values

while column<5 %loop to cover both real and imaginary computations
    
for x4=1:12%frequency loop 
        if x4==1
        freq=0.12;
    elseif x4==2
        freq=0.25;
    elseif x4==3
        freq=0.5;
    elseif x4==4
        freq=1;
    elseif x4==5
        freq=2;
    elseif x4==6
        freq=5;
    elseif x4==7
        freq=10;
    elseif x4==8
        freq=20;
    elseif x4==9
        freq=30;
    elseif x4==10
        freq=60;
    elseif x4==11
        freq=120;
    elseif x4==12
        freq=600;
         end

for S1=1:4 %z-value loop
   if S1==1
     z=2.500000e-01;
     Z1=0.125;
     Z2=0.25;
     Z3=0.375;
     dx=diff(0.125:0.125:0.375);
   elseif S1==2
     z=5.000000e-01;
     Z1=0.1375;
     Z2=0.50;
     Z3=0.625;
     dx=diff(0.375:0.125:0.625);
   elseif S1==3
     z=7.500000e-01;
     Z1=0.625;
     Z2=0.75;
     Z3=0.875;
     dx=diff(0.625:0.125:0.875);
   else
     z=1;
     Z1=0.875;
     Z2=1.00;
     Z3=1.125;
     dx=diff(0.875:0.125:1.125);
   end
for P2=0:30 %horizontal volumetric scan loop
   level=P2;
for P3=0:10 %vertical volumetric scan loop
   level2=P3;

mtx025a=load(sprintf('z = %d_vert%d_at_%d.txt',z,level,freq)); %loading a the full partical scan size (10+2 x 28+2) for E-field calculations
mtx025b=load(sprintf('z = %d_vert%d_at_%d.txt',z,level+1,freq));
mtx025c=load(sprintf('z = %d_vert%d_at_%d.txt',z,level+2,freq));

mtx025d=load(sprintf('z = %d_hori%d_at_%d.txt',z,level2,freq));
mtx025e=load(sprintf('z = %d_hori%d_at_%d.txt',z,level2+1,freq));
mtx025f=load(sprintf('z = %d_hori%d_at_%d.txt',z,level2+2,freq));

ra=mtx025b((1+(level2*10)):(10+(level2*10)),1); %selecting scan location and area from full particle scan
rb=mtx025e((1+(level*10)):(10+(level*10)),1);

mtxa=mtx025a((1+(level2*10)):(10+(level2*10)),column);
mtxb=mtx025b((1+(level2*10)):(10+(level2*10)),column);
mtxc=mtx025c((1+(level2*10)):(10+(level2*10)),column);

mtxd=mtx025d((1+(level*10)):(10+(level*10)),column);
mtxe=mtx025e((1+(level*10)):(10+(level*10)),column);
mtxf=mtx025f((1+(level*10)):(10+(level*10)),column); 

dif1b=(mtxb-mtxa)*(8.854*10^-12); %computation for charge density in z-direction
a1b=((dif1b.*ra)./((((ra.^2)+(z.^2)).^(3/2))));
Efield1b=sum((ra(2:10,:)-ra(1:9,:)).*(a1b(2:10,:)+a1b(1:9,:))/2);
Ed1b=Efield1b*(z/(2*8.854*10^-12)); %computation for E-field in r-direction

dif2b=(mtxc-mtxb)*(8.854*10^-12); %computation for charge density in r-direction
a2b=((dif2b.*ra)./((((ra.^2)+(z.^2)).^(3/2))));
Efield2b=sum((ra(2:10,:)-ra(1:9,:)).*(a2b(2:10,:)+a2b(1:9,:))/2);
Ed2b=Efield2b*(z/(2*8.854*10^-12)); %computation for E-field in z-direction

dif3b=(mtxe-mtxd)*(8.854*10^-12); %computation for charge density in z-direction
a3b=((dif3b.*rb)./((((rb.^2)+(z.^2)).^(3/2))));
Efield3b=sum((rb(2:10,:)-rb(1:9,:)).*(a3b(2:10,:)+a3b(1:9,:))/2);
Er1b=Efield3b*(z/(2*8.854*10^-12)); %computation for E-field in r-direction

dif4b=(mtxf-mtxe)*(8.854*10^-12); %computation for charge density in z-direction
a4b=((dif4b.*rb)./((((rb.^2)+(z.^2)).^(3/2))));
Efield4b=sum((rb(2:10,:)-rb(1:9,:)).*(a4b(2:10,:)+a4b(1:9,:))/2);
Er2b=Efield4b*(z/(2*8.854*10^-12)); %computation for E-field in r-direction
%-----------------------------------------------------------------------------------
dif1a=(mtxb-mtxa)*(8.854*10^-12); %computation for charge density in r-direction
a1a=((dif1a.*ra)./((((ra.^2)+(Z1.^2)).^(3/2))));
Efield1a=sum((ra(2:10,:)-ra(1:9,:)).*(a1a(2:10,:)+a1a(1:9,:))/2);
Ed1a=Efield1a*(Z1/(2*8.854*10^-12)); %computation for E-field in r-direction

dif2a=(mtxc-mtxb)*(8.854*10^-12); %computation for charge density in r-direction
a2a=((dif2a.*ra)./((((ra.^2)+(Z1.^2)).^(3/2))));
Efield2a=sum((ra(2:10,:)-ra(1:9,:)).*(a2a(2:10,:)+a2a(1:9,:))/2);
Ed2a=Efield2a*(Z1/(2*8.854*10^-12)); %computation for E-field in z-direction

dif3a=(mtxe-mtxd)*(8.854*10^-12); %computation for charge density in z-direction
a3a=((dif3a.*rb)./((((rb.^2)+(Z1.^2)).^(3/2))));
Efield3a=sum((rb(2:10,:)-rb(1:9,:)).*(a3a(2:10,:)+a3a(1:9,:))/2);
Er1a=Efield3a*(Z1/(2*8.854*10^-12)); %computation for E-field in r-direction

dif4a=(mtxf-mtxe)*(8.854*10^-12); %computation for charge density in z-direction
a4a=((dif4a.*rb)./((((rb.^2)+(Z1.^2)).^(3/2))));
Efield4a=sum((rb(2:10,:)-rb(1:9,:)).*(a4a(2:10,:)+a4a(1:9,:))/2);
Er2a=Efield4a*(Z1/(2*8.854*10^-12)); %computation for E-field in r-direction
%-----------------------------------------------------------------------------------
dif1c=(mtxb-mtxa)*(8.854*10^-12); %computation for charge density in z-direction
a1c=((dif1c.*ra)./((((ra.^2)+(Z3.^2)).^(3/2))));
Efield1c=sum((ra(2:10,:)-ra(1:9,:)).*(a1c(2:10,:)+a1c(1:9,:))/2);
Ed1c=Efield1c*(Z3/(2*8.854*10^-12)); %computation for E-field in r-direction

dif2c=(mtxc-mtxb)*(8.854*10^-12); %computation for charge density in r-direction
a2c=((dif2c.*ra)./((((ra.^2)+(Z3.^2)).^(3/2))));
Efield2c=sum((ra(2:10,:)-ra(1:9,:)).*(a2c(2:10,:)+a2c(1:9,:))/2);
Ed2c=Efield2c*(Z3/(2*8.854*10^-12)); %computation for E-field in z-direction

dif3c=(mtxe-mtxd)*(8.854*10^-12); %computation for charge density in z-direction
a3c=((dif3c.*rb)./((((rb.^2)+(Z3.^2)).^(3/2))));
Efield3c=sum((rb(2:10,:)-rb(1:9,:)).*(a3c(2:10,:)+a3c(1:9,:))/2);
Er1c=Efield3c*(Z3/(2*8.854*10^-12)); %computation for E-field in r-direction

dif4c=(mtxf-mtxe)*(8.854*10^-12); %computation for charge density in z-direction
a4c=((dif4c.*rb)./((((rb.^2)+(Z3.^2)).^(3/2))));
Efield4c=sum((rb(2:10,:)-rb(1:9,:)).*(a4c(2:10,:)+a4c(1:9,:))/2);
Er2c=Efield4c*(Z3/(2*8.854*10^-12)); %computation for E-field in r-direction

sigma=((8.854*10^-12)/0.06)*((Ed1b-Ed2b)+(Er1b-Er2b)); %final computation for volumetric charge in (1x1) square area
%nanocharge
Emagb=(sqrt((Er2b-Er1b)^2 + (Ed2b-Ed1b)^2)); %computation for E-field magnetude
Emaga=(sqrt((Er2a-Er1a)^2 + (Ed2a-Ed1a)^2));
Emagc=(sqrt((Er2c-Er1c)^2 + (Ed2c-Ed1c)^2));

NforF=[Emaga,Emagb,Emagc];

F=abs(diff(NforF)/dx); %computation for electrical force

nanocharge(abs(level2-11),level+1)=sigma; %organization in matrix for volumetric charge values
nanoelectricR(abs(level2-11),level+1)=(Er1b-Er2b);
nanoelectricZ(abs(level2-11),level+1)=(Ed1b-Ed2b);
nanoelectricMAG(abs(level2-11),level+1)=Emagb;
nanoforce(abs(level2-11),level+1)=F;
end
end

if column==3
realVC{(index)}=nanocharge; %cell assignment for volumetric charge data (real)
realER{(index)}=nanoelectricR; %cell assignment for r-direction E-field data (real)
realEZ{(index)}=nanoelectricZ; %cell assignment for z-direction E-field data (real)
realF{(index)}=nanoforce; %cell assignment for F-field data (real)
realmagE{(index)}=nanoelectricMAG; %cell assignment for total E-field data (real)

elseif column==4
imagVC{(index)}=nanocharge; %cell assignment for volumetric charge data (imag)
imagER{(index)}=nanoelectricR; %cell assignment for r-direction E-field data (imag)
imagEZ{(index)}=nanoelectricZ; %cell assignment for z-direction E-field data (imag)
imagF{(index)}=nanoforce; %cell assignment for F-field data (imag)
imagmagE{(index)}=nanoelectricMAG; %cell assignment for total E-field data (imag)
end

index=index+1;
end
end

column=column+1;
z=0; %resetting settings for starting arrays, and repeat
level=0;
level2=0;
freq=0;
index=1;
end
 
