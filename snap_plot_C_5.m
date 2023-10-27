clc;
clear;
%lam=0.001*(1:1:10);
lam2=linspace(0.01,0.035,31);
Nn=size(lam2);
N=Nn(2);
tt=zeros(N:4);
b=(1:8:100);


%lam3=0.01*(5:-0.1:1);
%Mn=size(lam3);
%M=Mn(2);
%ttb=zeros(M:4);
lam3=linspace(0,0.035,31);
D0=lam3(16);

for i=1:N
    %tt(i,:)=load(['/Users/ariji/Documents/MATLAB/fig5/tension/energy' sprintf('%f',lam(i)) '.txt']); 
   try
    tt21(i,:)=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/CD_c/energy' sprintf('%f',D0) 'a' sprintf('%f',lam2(i)) '.txt']); 
    
    catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', D0,lam2(i)))
        tt21(i,:)=tt21(i-1,:);
        


   end


   try
    tt(i,:)=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k_15/collar/energy' sprintf('%f',D0) 'a' sprintf('%f',lam2(i)) '.txt']); 
    
    catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', D0,lam2(i)))
        tt(i,:)=tt(i-1,:);
        


   end

   
   try
    tt42(i,:)=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k_20/collar/energy' sprintf('%f',D0) 'a' sprintf('%f',lam2(i)) '.txt']); 
    
    catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', D0,lam2(i)))
        tt42(i,:)=tt42(i-1,:);
        


   end

   
    try
    tt84(i,:)=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k25/k_15/collar/energy' sprintf('%f',D0) 'a' sprintf('%f',lam2(i)) '.txt']); 
    
    catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', D0,lam2(i)))
        tt84(i,:)=tt84(i-1,:);
        


    end

   
     try
    tt168(i,:)=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k30/k_20/collar/energy' sprintf('%f',D0) 'a' sprintf('%f',lam2(i)) '.txt']); 
    
    catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.3f D= %0.3f', D0,lam2(i)))
        tt168(i,:)=tt168(i-1,:);
        


   end

   
   
    %tt(i,:)=load(['/Users/ariji/Documents/MATLAB/fig5/k_1/energy' sprintf('%f',lam2(i)) '.txt']);  
    %tt84(i,:)=load(['/Users/ariji/Documents/MATLAB/fig5/k_L05/energy' sprintf('%f',lam2(i)) '.txt']);  
    %tt168(i,:)=load(['/Users/ariji/Documents/MATLAB/fig5/k_L1/energy' sprintf('%f',lam2(i)) '.txt']); 
end




%L=tt(:,1);
%W=tt(:,2);
%z=-tt(:,4);


figure
aa=plot(lam2,tt21(:,3),'--','color',[0.6350, 0.0780, 0.1840],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);
hold on;
bb=plot(lam2,tt(:,3),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);
cc=plot(lam2,tt42(:,3),':','color',[0.1, 0.55, 0.1],'linewidth',1.8);
dd=plot(lam2,tt84(:,3),'--','color',[0.3, 0.1, 0.6],'linewidth',2); %,'MarkerIndices',b,'Markersize',5);
ee=plot(lam2,tt168(:,3),'-p','color',[0.25, 0.25, 0.25],'linewidth',2,'MarkerIndices',b,'Markersize',3);
hold off;
%bb=plot(L(1:28),W(1:28),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);
%bb2=plot(L(33:end),W(33:end),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);

%bc=plot(tt84(:,1),tt84(:,2),'-','color',[0.1, 0.55, 0.1],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);%'color',[0.1, 0.4, 0.1],

%cc=plot(tt168(:,1),tt168(:,2),'-','color',[0.25, 0.25, 0.25],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);%'color',[0.1, 0.4, 0.1],

%legend([aa,bb,cc],'$\kappa=10~\mathrm{pN/nm}$','$\kappa=15~\mathrm{pN/nm}$','$\kappa=20~\mathrm{pN/nm}$');
axis([0.01 0.035  50 600]);
%hold off
xlabel('$C_0$ ($\mathrm{~nm}^{-1}$)', 'FontSize', 20);
ylabel('$E$ ($\mathrm{~pN-nm}$)', 'FontSize', 20);
set(gca,'FontSize',12);


figure
aa1=plot(lam2,-tt21(:,5),'--','color',[0.6350, 0.0780, 0.1840],'linewidth',1.8); % ,'MarkerIndices',b,'Markersize',5);
hold on;
bb1=plot(lam2,-tt(:,5),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); % ,'MarkerIndices',b,'Markersize',5);
cc1=plot(lam2,-tt42(:,5),':','color',[0.1, 0.55, 0.1],'linewidth',1.8);
dd1=plot(lam2,-tt84(:,5),'--','color',[0.3, 0.1, 0.6],'linewidth',2); %,'MarkerIndices',b,'Markersize',5);
ee1=plot(lam2,-tt168(:,5),'-p','color',[0.25, 0.25, 0.25],'linewidth',2,'MarkerIndices',b,'Markersize',3);
%bb1=plot(L(1:28),z(1:28),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);
%bb12=plot(L(33:end),z(33:end),'-','color',[0.3, 0.1, 0.6],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);
%bc1=plot(tt84(:,1),-tt84(:,4),'-','color',[0.1, 0.55, 0.1],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);%'color',[0.1, 0.4, 0.1],
%cc1=plot(tt168(:,1),-tt168(:,4),'-','color',[0.25, 0.25, 0.25],'linewidth',1.8); %,'MarkerIndices',b,'Markersize',5);%'color',[0.1, 0.4, 0.1],
legend([aa1,bb1,cc1],'$\kappa=10~\mathrm{pN/nm}$','$\kappa=15~\mathrm{pN/nm}$','$\kappa=20~\mathrm{pN/nm}$');
%legend([aa1,bb1,bc1,cc1],'$\lambda=0.001$','$\lambda=0.01$','$\lambda=0.05$','$\lambda=0.1$');
%daspect([1 1 1])
axis([0.01 0.035  0 250]);
legend([aa1,bb1,cc1,dd1,ee1],'$\kappa=10~\mathrm{pN/nm}$','$\kappa=15~\mathrm{pN/nm}$','$\kappa=20~\mathrm{pN/nm}$','$\kappa=25~\mathrm{pN/nm}$','$\kappa=30~\mathrm{pN/nm}$');%legend([aa1,bb1,bc1,cc1],'$\lambda=0.001$','$\lambda=0.01$','$\lambda=0.05$','$\lambda=0.1$');
hold off
xlabel('$C_0$ ($\mathrm{~nm}^{-1}$)', 'FontSize', 20);
ylabel('$l$ ($\mathrm{~nm}$)', 'FontSize', 20);
set(gca,'FontSize',12);


