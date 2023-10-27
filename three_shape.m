clc;
clear;
clf;
DRng=linspace(0,0.035,31);
CRng=linspace(0.01,0.035,31);

k=16;  %D
l=16;   %C
D0=DRng(k);
C0=CRng(l);
D1=DRng(k);
C1=CRng(l);
D2=DRng(k);
C2=CRng(l);

tube_21=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/CD/shape' sprintf('%f',D0) 'a' sprintf('%f',C0) '.txt']); 
tube_42=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k_15/shape' sprintf('%f',D1) 'a' sprintf('%f',C1) '.txt']); 
tube_84=load(['/Users/ariji/Documents/MATLAB/Itay/codes/files/gr_e5/k_20/shape' sprintf('%f',D2) 'a' sprintf('%f',C2) '.txt']); 
%tube_31=load('C:\Users\ariji\Documents\MATLAB\figure1\pannel_a\tubed_a15x42.txt');
%tube_36=load('C:\Users\ariji\Documents\MATLAB\figure1\pannel_a\tubed_a20x42.txt');
%tube_38=load('C:\Users\ariji\Documents\MATLAB\figure1\pannel_a\tubed_a25x42.txt');
%tube_39=load('C:\Users\ariji\Documents\MATLAB\figure1\pannel_a\tubed_a30x42.txt');
%tube_41=load('tube_41.txt');
%tube_51=load('tube_51.txt');
%tube_84=load('tube_84.txt');
%b=[1,25,51,80,120,160,200,240,280,320,360];
%a=[1,63,127,200,300,400,500,600,700,800,900];
R0=15;
figure
%aa=plot(R0*tube_21(:,1),R0*tube_21(:,2),':','color',[0.9290, 0.6940, 0.1250],'linewidth',4); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0]
cca=tube_21(:,3)/max(tube_21(:,3));


Nn=size(tube_21);
N=Nn(1);
kk=0;
xxa=R0*tube_21(:,1);
yya=R0*tube_21(:,2);
for i=1:N
    if(round(cca(i))>0.5)
        if(round(cca(i))>round(cca(i-1)))
          ksa=i;  
        end
        kea=i;
        %kk=kk+1;
        %ccx(kk)=xxa(i);
        %ccy(kk)=yya(i);
    end  
end
aac=plot(xxa(ksa:kea),yya(ksa:kea),'--','color',[0, 0.4470, 0.8410],'linewidth',1.8);  %[0, 0.4470, 0.8410][0.5, 0, 0.5]
hold on
aa1=plot(xxa(1:ksa),yya(1:ksa),'--','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
aa2=plot(xxa(kea:N),yya(kea:N),'--','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
axis([-200 200 -350 50]);
%hold off
xlabel('$r$ ($\mathrm{~nm}$)', 'FontSize', 20);
ylabel('$z $ ($\mathrm{~nm}$)', 'FontSize', 20);
set(gca,'FontSize',12);

%bb=plot(R0*tube_42(:,1),R0*tube_42(:,2),'-','color',[0.9290, 0.6940, 0.1250],'linewidth',4); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0]
ccb=tube_42(:,3)/max(tube_42(:,3));


Nnb=size(tube_42);
Nb=Nnb(1);
kkb=0;
xxb=R0*tube_42(:,1);
yyb=R0*tube_42(:,2);
for i=1:Nb
      if(round(ccb(i))>0.5)
        if(round(ccb(i))>round(ccb(i-1)))
          ksb=i;  
        end
        keb=i;
      end
end
%figure
bbc=plot(xxb(ksb:keb),yyb(ksb:keb),'-','color',[0, 0.4470, 0.8410],'linewidth',1.8);  %[0, 0.4470, 0.8410][0.5, 0, 0.5]
%hold on
bb1=plot(xxb(1:ksb),yyb(1:ksb),'-','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
bb2=plot(xxb(keb:Nb),yyb(keb:Nb),'-','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
%hold off
%legend([aa1,bb1],'$\lambda_0=0.0$','$\lambda_0=0.023$');
%legend([aa1,bb1],'$f_c=0$','$f_c=8~\mathrm{pN/nm}$');
daspect([1 1 1])
axis([-200 200 -350 50]);
%title('$\kappa=168 ~ \mathrm{pN \cdot nm}$')
%hold off
xlabel('$r$ ($\mathrm{~nm}$)', 'FontSize', 20);
ylabel('$z $ ($\mathrm{~nm}$)', 'FontSize', 20);
set(gca,'FontSize',12);


ccc=tube_84(:,3)/max(tube_84(:,3));

Nnc=size(tube_84);
Nc=Nnc(1);
kkc=0;
xxc=R0*tube_84(:,1);
yyc=R0*tube_84(:,2);
for i=1:Nc
      if(round(ccc(i))>0.5)
        if(round(ccc(i))>round(ccc(i-1)))
          ksc=i;  
        end
        kec=i;
      end
end
%figure
ccc=plot(xxc(ksc:kec),yyc(ksc:kec),':','color',[0, 0.4470, 0.8410],'linewidth',1.8);  %[0, 0.4470, 0.8410][0.5, 0, 0.5]
%hold on
cc1=plot(xxc(1:ksc),yyc(1:ksc),':','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
cc2=plot(xxc(kec:Nc),yyc(kec:Nc),':','color',[0.9290, 0.6940, 0.1250],'linewidth',1.8); %[0.9290, 0.6940, 0.1250] [0, 0.5, 0];
%hold off
%legend([aa1,bb1],'$\lambda_0=0.0$','$\lambda_0=0.023$');
%legend([aa1,bb1],'$f_c=0$','$f_c=8~\mathrm{pN/nm}$');
daspect([1 1 1])
axis([-200 200 -350 50]);
%title('$\kappa=168 ~ \mathrm{pN \cdot nm}$')
hold off
xlabel('$r$ ($\mathrm{~nm}$)', 'FontSize', 20);
ylabel('$z $ ($\mathrm{~nm}$)', 'FontSize', 20);
set(gca,'FontSize',12);


%legend([aa1,bb1,cc1],'$\kappa=10 \mathrm{~k_B T}$','$\kappa=15 \mathrm{~k_B T}$','$\kappa=20 \mathrm{~k_B T}$');