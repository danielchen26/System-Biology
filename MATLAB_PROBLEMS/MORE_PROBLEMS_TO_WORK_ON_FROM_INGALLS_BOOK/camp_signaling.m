% file camp_signaling.m
% model of cAMP signaling pathway in Dictyostelium
% described in Maeda et al (2004) Science 304 pp. 875-878. 
% problem 6.8.15

function camp_signaling

%declare model parameters
global k1;
global k2;
global k3;
global k4;
global k5;
global k6;
global k7;
global k8;
global k9;
global k10;
global k11;
global k12;
global k13;
global k14;


 %assign parameter values
 k1=2.0;
 k2=0.9;
 k3=2.5;
 k4=1.5;
 k5=0.6;
 k6=0.8;
 k7=1.0;
 k8=1.3;
 k9=0.3;
 k10=0.8;
 k11=0.7;
 k12=4.9;
 k13=23.0;
 k14=4.5;

 


ODEFUN=@campddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);
Tend=100;


%state vector:
% ACA=S(1);
% PKA=S(2);
% ERK2=S(3);
% REGA=S(4);
% CAMPi=S(5);
% CAMPe=S(6);
% CAR1=S(7);

[t,S]=ode45(ODEFUN, [0,Tend], [1 1 1 1 1 1 1]);


figure(1)
plot(t, S(:,5), t, S(:,2),t, S(:,3), 'k', 'LineWidth',3)
%axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')
legend('CAMPi', 'PKA', 'ERK2')

end


%dynamics 
function dS = campddt(t,S)

global k1;
global k2;
global k3;
global k4;
global k5;
global k6;
global k7;
global k8;
global k9;
global k10;
global k11;
global k12;
global k13;
global k14;

%define state variables
ACA=S(1);
PKA=S(2);
ERK2=S(3);
REGA=S(4);
CAMPi=S(5);
CAMPe=S(6);
CAR1=S(7);

%dynamics
dACAdt=k1*CAR1 - k2*ACA*PKA;
dPKAdt=k3*CAMPi-k4*PKA;
dERK2dt=k5*CAR1-k6*ERK2*PKA;
dREGAdt=k7-k8*REGA*ERK2;
dCAMPidt=k9*ACA-k10*REGA*CAMPi;
dCAMPedt=k11*ACA-k12*CAMPe;
dCAR1dt=k13*CAMPe-k14*CAR1;

dS=[dACAdt dPKAdt dERK2dt dREGAdt dCAMPidt dCAMPedt dCAR1dt]';

    
end


