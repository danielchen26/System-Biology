%file quorum_sensing.m
%model of quorum sensing mechanism of Vibrio fischeri
%adapted from James et al. (2000) J. of Molecular Biology 296 pp. 1127-1137.
%Figure 7.25


function quorum_sensing

clear all

%decalre model parameters
global k0
global n
global b
global a0
global RT
global k1
global k2
global nofeedflag
global diff
global popsize
global KM
global a


%assign parameter values
k1=0.5; % /muM^3 /min
k2=0.02; % /min 
b=0.07 ;% /min 
KM=0.01; % muM
a=10; % muM/min
diff=1000;
popsize=1000;
RT=0.5;
k0=0.0008;
n=0.6;
a0=0.05

%set flag for full model simualtion
nofeedflag=0;

%set simulation parameters
ODEFUN=@qsddt;
options=odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Refine', 1);
Tend=600;

%set initial condition: state vector is [A I Rstar A_ext]
x0=[0 0 0 0]';

%set mesh size (N=240 for figure 7.25: takes a few minutes to run)
N=24;
popmesh=zeros(1,N);
Ifeedmesh=zeros(1,N);
Inofeedmesh=zeros(1,N);
 
%run simulations over the mesh of popsize values, 
 for i=1:N+1
     popsize=10^(0 + 4*(i-1)/N);
     popsize=5000*(i-1)/N;
     nofeedflag=0;
     [t,s]=ode15s(ODEFUN, [0,Tend], x0, options);
     popmesh(i) = popsize;
     Ifeedmesh(i)=s(length(t),2); %record steady state I value
     nofeedflag=1;
     [t,s]=ode15s(ODEFUN, [0,Tend], x0, options);
     Inofeedmesh(i)=s(length(t),2);
 end
 
 %produce figure 7.25
 figure(4)
 set(gca,'fontsize',14)
 plot(popmesh, Ifeedmesh, 'k', popmesh, Inofeedmesh,'k--','Linewidth', 3)
 xlabel('Population density (arbitrary units)')
 ylabel('Intracellular LuxI concentration (arbitrary units)')
 legend('Original model', 'Modified model (no positive feedback)')
 


end

function dS=qsddt(t,s)

global k0
global n
global b
global a0
global RT
global k1
global k2
global nofeedflag
global diff
global popsize
global KM
global a


A=s(1); 
I=s(2);
Rstar=s(3);
Aout=s(4)

%original or modified dynamics for A
if nofeedflag==1
Addt= k0*15 -n*(A-Aout) - 2*k1*A^2*(RT-2*Rstar)^2 + 2*k2*Rstar;
else
Addt= k0*I -n*(A-Aout) - 2*k1*A^2*(RT-2*Rstar)^2 + 2*k2*Rstar;
end

Rstarddt= k1*A^2*(RT-2*Rstar)^2 - k2*Rstar;  
Iddt= -b*I + a0 + a*Rstar/(KM+Rstar);
Aoutddt=popsize*n*(A-Aout)-diff*Aout;

dS =[Addt, Iddt, Rstarddt Aoutddt]'; 

end

