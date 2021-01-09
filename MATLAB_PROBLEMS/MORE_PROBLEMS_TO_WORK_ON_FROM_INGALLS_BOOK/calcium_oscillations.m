% file calcium_oscillations.m
% model of calcium-induced calcium release in hepatocytes
% adapted from Othmer and Tang (1997) in 'Case studies in mathematical
% modelling', Othmer et al. eds, Prentice Hall
% Figure 6.18

function calcium_oscillations

%declare parameters
global gamma0 ;
global gamma1;
global Cs;
global p1 ;
global p2;
global k1;
global km1;
global k2;
global km2 ;
global k3;
global km3 ;
global I;
global vr;

%mode=1 for Figure 6.18A (steady input I)
%mode=2 for Figure 6.18B (time-varying input I)
global mode

%assign parameter values
gamma0=0.1;
gamma1=20.5; 
p1=8.5;
p2=0.065;
k1=12;
k2=15;
k3=1.8;
km1=8;
km2=1.65;
km3=0.21;
Cs=8.37; 
vr=0.185;
I=1;


ODEFUN=@dSdt;

% mode 1: constant input 
mode=1;

%set initial conditions.  State is S = [ C R RI RIC RICC]
S0=[0 1 0 0 0];

%run simualtion
[t,S]=ode45(ODEFUN, [0,25], S0);

%generate Figure 6.18A
figure(1)
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,4), 'k:', t, S(:,5), 'k--', 'Linewidth', 2)
legend('cytosolic [Ca^{2+}] (\mu M)', 'fraction of open channels (RIC+)', 'fraction of closed channel (RIC+C-)')
xlabel('Time (s)')
ylabel('Abdundance')
str1(1) = {'A'};
text(-3,2.5,str1, 'Fontsize', 40)


% mode 2: varying input
mode=2

%set initial condition
S0=[0 1 0 0 0];

%run simulation
[t,S]=ode45(ODEFUN, [0,120], S0);

%generate Figure 6.18B
figure(2)
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', 'Linewidth', 2)
xlabel('Time (s)')
ylabel('cytosolic [Ca^{2+}] (\mu M)')
axis([0 120 0 2.5])
str1(1) = {'B'};
text(-15,2.5,str1, 'Fontsize', 40)


end

%dynamics
function dS = dSdt(t,S)

global gamma0 ;
global gamma1;
global Cs;
global p1 ;
global p2;
global k1;
global km1;
global k2;
global km2 ;
global k3;
global km3 ;
global vr;
global I;
global mode

%if mode=2, set the input parameter according to a time-varying profile:

if mode==2
    
   I=0;

   if (t>20)
     I=0.7;
   end

   if (t>60)
     I=1.2;
   end

   if (t>90)
     I=4;
   end
end

%dynamics
dS=[vr*(gamma0+gamma1*S(4))*(Cs-S(1)) - (p1*S(1)^4)/(p2^4+S(1)^4);-k1*I*S(2)+km1*S(3); -(km1+k2*S(1))*S(3)+ k1*I*S(2) + km2*S(4); -(km2+k3*S(1))*S(4) + k2*S(1)*S(3) + km3*S(5); k3*S(1)*S(4) - km3*S(5)];


end


