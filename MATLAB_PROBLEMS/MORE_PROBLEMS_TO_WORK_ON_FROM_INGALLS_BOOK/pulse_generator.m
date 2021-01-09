%file pulse_generator.m
%model of synthetic pulse generating system
%adapted from Basu et al. (2004) PNAS 101 6355-6360
%Figure 7.28

function pulse_generator

clear all

%decalre model parameters
global aG
global bG
global KC
global aC
global bC
global k1
global k2
global KR
global RT
global A

%assign parameter values
aG = 80; % muM/min
bG = 0.07; % /min 
 KC= 0.008; % muM
 aC= 0.5; % muM/min
 bC = 0.3; % /min 
 k1 = 0.5; % /muM^3 /min
 k2 = 0.02; % /min 
 KR= 0.02; % muM
 RT = 0.5; % muM


%set simulation parameters
ODEFUN=@pulseddt;
options=odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Refine', 1);
Tend = 200;

%set initial condition: state vector is [G C R]
x0=[0 0 0]';
 
%run simulation
[t,s]=ode15s(ODEFUN, [-10,Tend], x0, options);

%generate figure 7.28
 figure(1)
 set(gca,'fontsize',14)
 plot(t,s(:,1), 'k', t,s(:,2), 'k--',t,s(:,3), 'k-.', 'Linewidth', 3)
 legend('GFP', 'cI', 'LuxR:AHL complex')
 axis([-5 50 0 2])
 xlabel('Time (minutes)')
 ylabel('Intracellular concentration (arbitrary units)')


end

%dynamics
function dS=pulseddt(t,s)

global aG
global bG
global KC
global aC
global bC
global k1
global k2
global KR
global RT
global A


%figure in text:
if t < 0
    A=0;
else
    A=10;
end

G=s(1);
C=s(2);
R=s(3);



 Gddt = - bG*G + aG*((R/KR)/(1+(R/KR) + (C/KC)^2 + (R/KR)*(C/KC)^2));
 Cddt = - bC*C + aC*R/(KR+R);
 Rddt = - k2*R + k1*[RT-2*R]^2*A^2; 

    

dS =[Gddt, Cddt, Rddt]';

end

