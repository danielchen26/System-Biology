%file band_detector.m
%model of synthetic band detector system
%adapted from Basu et al. (2005) Nature 434 177-193.
%Figure 7.30


function band_detector

clear all

%declare model parameters
global aG
global bG
global KL
global aL1
global aL2
global KC
global bL
global aC
global bC
global k1
global k2
global KR
global RT
global A

%assign parameter values
aG = 2; % muM/min
 bG = 0.07; % /min 
 KL= 0.8; % muM
 aL1= 1; % muM/min
 aL2= 1; % muM/min
 KC= 0.008; % muM
 bL = 0.02; % /min 
 aC= 1; % muM/min
 bC = 0.07 ;% /min 
 k1 = 0.5; % /muM^3 /min
 k2 = 0.02; % /min 
 KR= 0.01; % muM
 RT = 0.5 % muM
 A=0.1;

%set simulation parameters
ODEFUN=@bandddt;
options=odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Refine', 1);
Tend = 300;


%set initial condition, state vector is [G L C R]
x0=[0 0 0 0]';
 
%set grid-size 
 N=90
 Amesh=zeros(1,N);
 gmesh=zeros(1,N);
 
 %over the grid of inputs, simulate and record steady-state values
 for i=1:N+1
     i
     A=10^(-4 + 4*(i-1)/N);
     [t,s]=ode15s(ODEFUN, [0,Tend], x0, options);
     Amesh(i) = A;
     gmesh(i)=s(length(t),1);
     lmesh(i)=s(length(t),2);
     cmesh(i)=s(length(t),3);
     rmesh(i)=s(length(t),4);
 end
 
 %produce figure 7.30
 figure(2)
 set(gca,'fontsize',14)
 semilogx(Amesh, gmesh, 'k', Amesh, lmesh, 'k--',Amesh, cmesh, 'k-.', 'Linewidth', 3)
 xlabel('Extracellular AHL concentration (\mu M)')
 ylabel('Intracellular concentration (\mu M)')
 legend('GFP', 'LacI', 'cI')

 

end

%dynamics
function dS=bandddt(t,s)

global aG
global bG
global KL
global aL1
global aL2
global KC
global bL
global aC
global bC
global k1
global k2
global KR
global RT
global A


G=s(1);
L=s(2);
C=s(3);
R=s(4);



 Gddt = -bG*G + aG*(1/(1+(L/KL)^2));
 Lddt = -bL*L + aL1/(1+(C/KC)^2) + aL2*R/(KR+R);
 Cddt = -bC*C + aC*R/(KR+R);
 Rddt = -k2*R + k1*(RT-2*R)^2*A^2; 

    

dS =[Gddt,  Lddt, Cddt, Rddt]';

end

