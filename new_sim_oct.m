clc
clear
close all
% Input=======================
area_domain =400; %% dimensionless area
mesh=(0:0.0001:1).^4; %% mesh points 
lambda = 0.01;  %% membrane tension at the end \lambda_0
R0 = 15;  %% scale for length
gamma =20;   %% exponent of tanh function
k0 =42; %bending rigidity
dk =1;  % ration of bending rigidity with protein coat to the bare membrane
dkd=0; % fraction change of deviatoric bending rigidity with protein coat to the bare membrane 
P = 0; % pressure
dkG = 0; % fraction change in Gaussian bending rigidity with protein coat to the bare membrane 
%rIn=0;
%rOut=0;
%step=60;
%r1=linspace(1,5900/(2*pi*R0^2),step/3);
%r2=linspace(5900/(2*pi*R0^2),490000/(2*pi*R0^2),step/3);
%r3=492000/(2*pi*R0^2)*ones(1,step/3);
%rF=[r1 r2 r3];

%f0=0;
acoat=10;  %2*600000/(2*pi*R0^2); %% area of protein coat
%c1=zeros(1,2*step/3);
%c2=1*linspace(0,0.004,step/3);
%C0=0.02; %0.0265;  %0.65*0.0075;       %[c1 c2];
%C0=0.65*0.0075;
%D0=0.00;   %0.75*0.0075;
%C0=0.02-D0;
D1=0.0; %0.028;  %% deviatoric curvature at collar
acn=acoat/10;  %% collar area
fc=0.0*7.95/(acn*2*pi*R0^2);  % collar force  %%
%z1=linspace(0,1200,step/3);
%z2=1200* ones(1,2*step/3);
%zpRng=[z1 z2];
%lamRng=lambda*(5:-1:1);



initSol =endoInit(area_domain,mesh, lambda, k0,R0);      %% initialization

%for ki=1:5
DRng=linspace(0,0.035,31);      %% deviatoric curvature
CRng=linspace(0.01,0.035,31);  %% spontaneous curvature

for i=1:length(DRng)
for j=1:length(CRng)

    %D0=4*42*DRng(i)/(k0+4*42);
    %C0=0.035; %-D0;
    %C0=0.02-DRng(i);

try
[t, Sol]= tubedBAR_solver_new_par(area_domain,mesh, lambda, fc, acn, acoat,k0, dk,dkd,dkG, P, gamma, CRng(j), DRng(i), D1, R0, initSol);
catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', CRng(j),DRng(i)))
        Sol=initSol;
        


end

initSOl=Sol;
  %  fc=fc*1.1
%Nn=size(Sol);
%N=Nn(2);
%xx(N)=Sol(1,1);
%yy(N)=Sol(2,1);

%for i=1:N-1
   %xx(N-i)=-Sol(1,i+1);
   %xx(N+i)=Sol(1,i+1);
  % yy(N-i)=Sol(2,i+1);
 %  yy(N+i)=Sol(2,i+1);
%end

%plot(xx,yy)

%pause(1)
end
end