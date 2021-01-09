%file apoptosis.m
%Model of apoptosis signalling pathway
%modified from Eissing et al. (2004) Journal of Biological Chemistry 279, pp.
%36892-36897
%Figure 6.16

function apoptosis

%declare parameters
global k1 ;
global k2;
global k3;
global k4;
global k5;
global k6;
global k7;
global k8 ;
global k9 ;
global k10 ;
global k11 ;
global k12 ;
global k13;
global k14 ;
global k15 ;
global k16;
global k17;
global k18;
global k19;
global input;

%assign paramter values
k1=507; 
k2=3.9e-3;
k3=1e-5;  
k4=5.8e-3; 
k5=5e-4; 
k6=0.21; 
k7=81.9; 
k8=3.9e-3; 
k9=5.8e-6; 
k10=5.8e-3;
k11=5e-4; 
k12=0.21; 
k13=40; 
k14=1e-3;
k15=464; 
k16=1.16e-2; 
k17=3e-4; 
k18=1.16e-2;
k19=1.73e-2; 

input=0;



%assign initial condition.  These values were taken from a previous
%simulation run to steadt state with input=0
%state vector S=[C8 C8star C3 C3star BAR IAP C8starBAR C3starIAP]
 S0=1e5*[ 1.30 0 0.21 0 0.4 0.4 0 0];
 
%run simulation
ODEFUN=@dSdt;
[t,S]=ode45(ODEFUN, [0,1800], S0);

%generate figure 6.16
figure(1)
set(gca,'fontsize',14)
plot(t, S(:,2), 'k--', t, S(:,4), 'k', 'Linewidth', 3)
legend('C8*', 'C3*')
axis([0 1600 0 100000])
xlabel('Time (min)')
ylabel('Concentration (molecules per cell)')


end

%dynamics
function dS = dSdt(t,S)

global k1 ;
global k2;
global k3;
global k4;
global k5;
global k6;
global k7;
global k8 ;
global k9 ;
global k10 ;
global k11 ;
global k12 ;
global k13;
global k14 ;
global k15 ;
global k16 ;
global k17 ;
global k18 ;
global k19;

global input;


%assign time-varying profile for parameter 'input'
input=0;

   if t>100 && t<100+1100
       input=200;
   end

%assign state variables
C8=S(1);
C8star=S(2);
C3=S(3);
C3star=S(4);
BAR=S(5);
IAP=S(6);
C8starBAR=S(7);
C3starIAP=S(8);

%assign dynamics
dS=[k1 - k2*C8 - k3*(C3star+input)*C8;
    k3*(C3star+input)*C8- k4*C8star - k5*C8star*BAR + k6*C8starBAR;
    k7-k8*C3 - k9*C8star*C3;
    k9*C8star*C3 - k10*C3star-k11*C3star*IAP+k12*C3starIAP;
    k13 - k5*C8star*BAR + k6*C8starBAR - k14*BAR;
    k15 - k11*C3star*IAP+k12*C3starIAP - (k16+k17*C3star)*IAP;
    k5*C8star*BAR - k6*C8starBAR - k18*C8starBAR;
    k11*C3star*IAP - k12*C3starIAP - k19*C3starIAP];


end

