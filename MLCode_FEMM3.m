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

nanocharge=zeros(30,30); %size of full particle figure for volumetric charge density
nanoelectricR=zeros(30,30); %size of full particle figure for r-direction E-field
nanoelectricZ=zeros(30,30); %size of full particle figure for z-direction E-field
nanoforce=zeros(30,30); %size of full particle figure for F-field
nanoelectricMAG=zeros(30,30); %size of full particle figure for total E-field magnetude
%nanophi=zeros(30,30); %size of full particle figure for phase

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
%phi=cell(1,48); %cell for phase values

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

for x1=1:4 %z-value loop
   if x1==1
    z=2.500000e-01;
   elseif x1==2
     z=5.000000e-01;
   elseif x1==3
     z=7.500000e-01;
   else
     z=1;
   end

for x2=0:29 %horizontal volumetric scan loop
   level=x2;
for x3=0:11 %vertical volumetric scan loop
   level2=x3;

mtx025a=load(sprintf('z = %d_vert%d_at_%d.txt',z,level,freq)); %loading a the full partical scan size (11+1 x 29+1) for E-field calculations
mtx025b=load(sprintf('z = %d_vert%d_at_%d.txt',z,level+1,freq));
mtx025c=load(sprintf('z = %d_hori%d_at_%d.txt',z,level2,freq));
mtx025d=load(sprintf('z = %d_hori%d_at_%d.txt',z,level2+1,freq));

ra=mtx025a((1+(level2*10)):(10+(level2*10)),1); %selecting scan location and area from full particle scan
rb=mtx025c((1+(level*10)):(10+(level*10)),1);
mtxa=mtx025a((1+(level2*10)):(10+(level2*10)),column);
mtxb=mtx025b((1+(level2*10)):(10+(level2*10)),column);
mtxc=mtx025c((1+(level*10)):(10+(level*10)),column);
mtxd=mtx025d((1+(level*10)):(10+(level*10)),column);
    
dif3=(mtxb-mtxa)*(8.854*10^-12); %computation for charge density in z-direction
a2=((dif3.*ra)./((((ra.^2)+(z.^2)).^(3/2))));
Efield=sum((ra(2:10,:)-ra(1:9,:)).*(a2(2:10,:)+a2(1:9,:))/2);
Ed=Efield*(z/(2*8.854*10^-12)); %computation for E-field in z-direction

dif4=(mtxd-mtxc)*(8.854*10^-12); %computation for charge density in r-direction
a3=((dif4.*rb)./((((rb.^2)+(z.^2)).^(3/2))));
Efield2=sum((rb(2:10,:)-rb(1:9,:)).*(a3(2:10,:)+a3(1:9,:))/2);
Er=Efield2*(z/(2*8.854*10^-12)); %computation for E-field in z-direction

Emag=(sqrt(Er^2 + Ed^2));

sigma=((8.854*10^-12)/1.92)*((Ed)+(Er)); %final computation for volumetric charge in (1x1) square area
%nanocharge
pointy=(0:1:30);
pointx=(0:1:30);

F=Ed/Er;

%P=1;

nanocharge(abs(level2-12),level+1)=sigma; %organization in matrix for volumetric charge values
nanoelectricR(abs(level2-12),level+1)=Er;
nanoelectricZ(abs(level2-12),level+1)=Ed;
nanoelectricMAG(abs(level2-12),level+1)=Emag;
nanoforce(abs(level2-12),level+1)=F;
%nanophi(abs(level2-12),level+1)=P;
end
end

%phi{(index)}=nanophi; %cell assignment for phase shift data

if column==3
realVC{(index)}=nanocharge; %cell assignment for volumetric charge data (real)
realER{(index)}=nanoelectricR; %cell assignment for r-direction E-field data (real)
realEZ{(index)}=nanoelectricZ; %cell assignment for z-direction E-field data (real)
realF{(index)}=nanoforce; %cell assignment for F-field data (real)
realmagE{(index)}=nanoelectricMAG; %cell assignment for total E-field data (real)

else
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