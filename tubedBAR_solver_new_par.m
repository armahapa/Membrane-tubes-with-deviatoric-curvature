%% Solve the membrane shape with spontaneous mean and deviatoric curvature 

%%

% Inputs:
%   alpha - ND membrane area
%   mesh - meshing for the domain, 
%   lambda - membrane tension at the boundary (pN/nm)
%   alpha0 - coat area of protein (ND)
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   dkd - ration of kappa_d with protein to kappa_d with bare membrane -1
%   dkG - ration of kappa_Gaussian with protein to kappa_Gaussian with bare membrane -1
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the solution

% Outputs:
%   t - area mesh points
%   Sol - solution array


function [t, Sol] = tubedBAR_solver_new_par(alpha, mesh, lambda, fc, acn, alpha0, k0, dk, dkd, dkG, P, gamma, C0, D0, D1, R0, initSol)

t=alpha*mesh;   % area mesh points

% declare and assign global variables to be used in nested functions
global a acap f gt a0 ac iSol lam c0 d0 d1 g delta deltad delkG p pw wd

a = alpha;
a0 = alpha0;
ac=acn/2;
gt = t;
lam = lambda*R0^2/k0;   % dimensionless membrane tension
c0 = R0*C0;             % dimensionless preferred curvature
d0=R0*D0;
d1=R0*D1;
g = gamma;
delta = dk;
deltad=dkd;
delkG = dkG;
p = P*R0^3/k0;          % dimensionless pressure
f=fc*R0^3/k0;
pw = 1*R0^3/k0;         % dimensionless pressure from rigid wall
wd = 10/R0;             % dimensionless wall "distance"
iSol = initSol;
acap=1/((c0+d0))^2;
% initial guess structure
solinit = bvpinit(t,@mat4init);

% solver options; increasing maximum number of mesh points solver will use
options = bvpset('NMax', 100*length(t), 'RelTol', 1e-3);

% solve the boundary value problem
sol = bvp4c(@mat4ode,@mat4bc,solinit,options);

% evaluate the solution on the original mesh points
Sol = deval(sol,t);

% plot the resultant profile of the membrane
coatArea = [0 alpha0];
xLim = [-sqrt(10*alpha0)*R0 sqrt(10*alpha0)*R0];
%xLim = [-350 350];
yLim = [-2*(c0+d0)*a0*R0 R0];
%yLim = [-100 300];

file_name=['shape' sprintf('%f',D0) 'a' sprintf('%f',C0)  '.txt'];
fileID4=fopen(file_name,'w+');

file_name=['energy' sprintf('%f',D0) 'a' sprintf('%f',C0)  '.txt'];
fileID5=fopen(file_name,'w+');


C = 0.5*c0*(tanh(g*(t+ a0)) - tanh(g*(t - a0)));

c = 0.5*c0*(tanh(g*(t+ a0)) - tanh(g*(t - a0)));
m= 0.5*d0*(tanh(g*(t +a0)) - tanh(g*(t - a0)))-0.5*d0*(tanh(g*(t+ acap)) - tanh(g*(t - acap)))+0.5*d1*(tanh(g*(t-(a0-ac))) - tanh(g*(t - (a0+ac)))); 

fcl=0.5*f*(tanh(g*(t-(a0-ac))) - tanh(g*(t - (a0+ac)))); 
b = 1 + 0.5*(delta-1)*(tanh(g*(t + a0)) - tanh(g*(t - a0)));

bd= 1+0.5*(deltad)*(tanh(g*(t + a0)) - tanh(g*(t - a0)));

W=b.*(Sol(4,:)-c).^2+bd.*(sin(Sol(3,:))./(Sol(1,:))-Sol(4,:)-m).^2;
Nn=size(Sol);
N=Nn(2);
xx=zeros(2*N-1);
yy=zeros(2*N-1);
xx(N)=Sol(1,1);
yy(N)=Sol(2,1);
cc(N)=C(1);
ww(N)=W(1);
HH(N)=Sol(4,1);
llam(N)=Sol(6,1);

W1=0;
W2=0;

for i=1:N-1
   xx(N-i)=-Sol(1,i+1);
   xx(N+i)=Sol(1,i+1);
   yy(N-i)=Sol(2,i+1);
   yy(N+i)=Sol(2,i+1);
   cc(N-i)=C(i+1);
   cc(N+i)=C(i+1);
   ww(N-i)=W(i+1);
   ww(N+i)=W(i+1);
   HH(N-i)=Sol(4,i+1);
   HH(N+i)=Sol(4,i+1);
   llam(N-i)=Sol(6,i+1);
   llam(N+i)=Sol(6,i+1);
   W1=W1+W(i)*(t(i+1)-t(i));  %bending
   W2=W2+(W(i)+Sol(6,i))*(t(i+1)-t(i));  %bending plus tension
end

for i=1:2*N-1
    
   fprintf(fileID4,'%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f',xx(i),yy(i),cc(i),ww(i),HH(i),llam(i)); 
   fprintf(fileID4,'\n');
   
  % writing files for x,y,C,W, H, lam 
end


  % fprintf(fileID4,'\n');

fclose(fileID4);


fprintf(fileID5,'%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f ',lambda, k0, W1*k0*2*pi,W2*k0*2*pi,R0*yy(N)); 
 %fprintf(fileID5,'\n');
fclose(fileID5);
% writing files for lambda_0,k_0, bending energy density, total energy
% density, length of tube
plotTitle = sprintf(''); %('k = %0.3f, D = %0.4f', k0, C0);
plotMemProfileArea(Sol, t, R0, coatArea, [], [], xLim, yLim, plotTitle, 0);


%%Define the variables
% X(1)=x;
% X(2)=y;
% X(3)=phi;
% X(4)=H;
% X(5)=L;

%%%%%the differential equations
%------------------------------------
function dXdt = mat4ode(t, X)
%parameters
global a acap f c0 d0 d1 a0 ac g delta deltad delkG p pw wd

% spontaneous curvature
%c = 0.5*c0*(1 - tanh(g*(t - a0)));

C = 0.5*c0*(tanh(g*(t+ a0)) - tanh(g*(t - a0)))+0.5*d0*(tanh(g*(t+ acap)) - tanh(g*(t - acap)));
M = 0.5*d0*(tanh(g*(t +a0)) - tanh(g*(t - a0)))-0.5*d0*(tanh(g*(t+ acap)) - tanh(g*(t - acap)))+0.5*d1*(tanh(g*(t-(a0+ac))) - tanh(g*(t - (a0+2*ac)))) ; 
fcl=0.5*f*(tanh(g*(t-(1.5*a0-ac))) - tanh(g*(t - (1.5*a0+ac)))); 
% derivative of spontaneous curvature
%dc = 0.5*c0*g*(tanh(g*(t - a0))^2 - 1);
dC = 0.5*c0*g*(tanh(g*(t - a0))^2 - tanh(g*(t + a0))^2)+0.5*d0*g*(tanh(g*(t - acap))^2 - tanh(g*(t + acap))^2);
dM = 0.5*d0*g*(tanh(g*(t - a0))^2 - tanh(g*(t + a0))^2)-0.5*d0*g*(tanh(g*(t - acap))^2 - tanh(g*(t + acap))^2)+0.5*d1*g*(tanh(g*(t - (a0+ac)))^2 - tanh(g*(t - (a0+2*ac)))^2);

% bending modulus
%b = 1 + 0.5*(delta-1)*(1 - tanh(g*(t - a0)));
b = 1 + 0.5*(delta-1)*(tanh(g*(t + a0)) - tanh(g*(t - a0)));

bd= 1+0.5*(deltad)*(tanh(g*(t + a0)) - tanh(g*(t - a0)));
% derivative of bending modulus
%db = 0.5*(delta-1)*g*(tanh(g*(t - a0))^2 - 1);
db = 0.5*(delta-1)*g*(tanh(g*(t - a0))^2 - tanh(g*(t + a0))^2);

dbd = 0.5*(deltad)*g*(tanh(g*(t - a0))^2 - tanh(g*(t + a0))^2);
% Derivative of Gaussian modulus
%dkg = 0.5*delkG*g*(tanh(g*(t - a0))^2 - 1);
dkg = 0.5*delkG*g*(tanh(g*(t - a0))^2 - tanh(g*(t + a0))^2);
% Second derivative of Gaussian modulus
%ddkg = delkG*g^2*tanh(g*(t - a0))*(1 - tanh(g*(t - a0))^2);
ddkg = delkG*g^2*(tanh(g*(t - a0))*(1 - tanh(g*(t - a0))^2)-tanh(g*(t + a0))*(1 - tanh(g*(t + a0))^2));

% opposing pressure from rigid wall
%dp = -p*(tanh(g*(X(2)))) - pw/2*(1 + tanh(g*(X(2) - wd)));

% no wall to oppose pressure
dp = p;

k2=1;
k1=1;
kg=-0;

% ODEs - see Hassinger et al, 2016
dXdt = [cos(X(3))/X(1)
        sin(X(3))/X(1)
        (2*X(1)*X(4)-sin(X(3)))/X(1)^2
        X(5)/((k1*b+k2*bd)*X(1)^2)-k1*db*(X(4)-C)/(k1*b+k2*bd)+k2*dbd*(sin(X(3))/X(1)-X(4)-M)/(k1*b+k2*bd)+(k1*b*dC-k2*bd*dM)/(k1*b+k2*bd)+k2*bd/(k1*b+k2*bd)*(cos(X(3))/X(1)^2)*2*(X(4)-sin(X(3))/X(1))-(dkg)*(kg/(k1*b+k2*bd))*(sin(X(3))/X(1))+(cos(X(3))/X(1)^2)*(2*k2*bd*(sin(X(3))/X(1)-X(4)-M)/(k1*b+k2*bd)) %     (1+k2*bd/(2*k1*b))*X(5)/X(1)^2-k2*bd/(k1*b)*X(4)*cos((X(3)))/X(1)^2+k2*bd/(k1*b)*sin(X(3))*cos(X(3))/X(1)^3+k2*bd/(2*k1*b)*dM/X(1)-(k2*bd/(2*k1*b))*cos(X(3))*(sin(X(3))/X(1)-X(4))/(X(1))^2-(db*2*k1+dbd*k2)/(b*2*k1)*(X(4)-C)-(k2*dbd/(2*k1*b))*(sin(X(3))/X(1)-X(4)-M)-(dkg/b)*(kg/(2*k1))*(sin(X(3))/X(1))+(1+k2*bd/(2*k1*b))*dC
        (fcl*sin(X(3))-2*k1*b*(X(4)-C)*(X(4)^2+(X(4)-sin(X(3))/X(1))^2)+2*X(4)*(k1*b*(X(4)-C)^2+k2*bd*(sin(X(3))/X(1)-X(4)-M)^2+X(6)-2*k2*bd*(sin(X(3))/X(1)-X(4)-M)*(sin(X(3))/X(1)-X(4)))-(k1*db+k2*dbd)*X(5))/(k1*b+k2*bd)%            -((db*2*k1+dbd*k2)/(b*2*k1+bd*k2))*X(5)-2*(X(4))*(X(4)^2+(X(4)-sin(X(3))/X(1))^2)+2*X(4)*((X(4)-C)^2-k2*bd/(2*k1*b+k2*bd)*(sin(X(3))/X(1)-X(4)-M)^2+X(6)+2*k2*bd/(2*k1*b+k2*bd)*(sin(X(3))/X(1)-X(4))*(sin(X(3))/X(1)-X(4)-M)) %- 0   %-2*cos(X(3))/X(1)*(k2/(2*k1)*X(5)/X(1)-k2/k1*X(4)*cos(X(3))/X(1)+k2/k1*sin(X(3))*cos(X(3))/X(1)^2+k2/(2*k1)*dM);
         fcl*cos(X(3))/X(1)+2*k1*b*(X(4)-C)*dC+2*k2*bd*(sin(X(3))/X(1)-X(4)-M)*dM-k1*db*(X(4)-C)^2-k2*dbd*(sin(X(3))/X(1)-X(4)-M)^2-dkg*kg*(X(4)^2-(sin(X(3))/X(1)-X(4))^2) %+dbd*k2/(2*k1*b+k2*bd)*(sin(X(3))/X(1)-X(4)-M)^2-2*k2*bd/(2*k1*b+k2*bd)*(sin(X(3))/X(1)-X(4)-M)*dM-(db*2*k1+dbd*k2)/(b*2*k1+bd*k2)*(X(4)-C)^2-dkg*kg/(2*k1+k2)*(X(4)^2-(sin(X(3))/X(1)-X(4))^2)+2*(X(4)-C)*dC
        ];% +db*k2/(2*k1+k2)*(sin(X(3))/X(1)-X(4)-M)^2-db*(X(4))^2-dkg*kg*(X(4)-(sin(X(3))/X(1)-X(4))^2)  ];
            

%-------------------------boundary conditions-------------

function res = mat4bc(Xa,Xb) 
global lam
ds = 1e-4;  % small offset to prevent division by 0
   
    % boundary conditions - see Hassinger et al, 2016
    res = [ Xa(1) - ds
            Xb(2) - 0         
            Xa(3)             
            Xb(3)
            Xa(5)
            Xb(6) - lam
            ];
        
        
%-----------------------------------Initial guesses------------


function Xinit = mat4init(t)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Xinit = iSol(:,find(gt==t));