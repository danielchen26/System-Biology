%file passive_membrane.m
%model of the response of a passive membrane to perturbations in voltage
%Figure 8.3

function passive_membrane

%declare parameters
global C;
global g;
global V_Nernst_K;

%assign parameter values
C= 0.98; %microfarad/cm^2 
g=0.0144; % millisiemens/ cm^2
V_Nernst_K = -93.6; %millivolts 

%set simulation parameters
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);
ODEFUN=@passive_membrane_ddt;
Tend=800;

%set initial condition (resting voltage)
S0=[-93.6];

%run multiple simulations, from perturbed points
[t1,S1]=ode15s(ODEFUN, [0 500], S0, OPTIONS);
[t2, S2]=ode15s(ODEFUN, [t1(length(t1)) Tend+ t1(length(t1))], [-63.6], OPTIONS);
[t3, S3]=ode15s(ODEFUN, [t2(length(t2)) Tend+ t2(length(t2))], [-123.6], OPTIONS);
[t4, S4]=ode15s(ODEFUN, [t3(length(t3)) Tend+ t3(length(t3))], [-78.6], OPTIONS);
[t5, S5]=ode15s(ODEFUN, [t4(length(t4)) Tend+200+ t4(length(t4))], [-108.6], OPTIONS);

%produce Figure 8.3
figure(1)
set(gca, 'fontsize', 14)
plot([t1; t2; t3; t4; t5], [S1;S2;S3;S4;S5], 'k', 'Linewidth', 3)
xlabel('Time (msec)', 'fontsize',12)
ylabel('Membrane Voltage (mV)', 'fontsize',12)
axis([0 3800 -130 -50])

end

%dynamics
function dSdt = passive_membrane_ddt(t,S)

global C;
global g;
global V_Nernst_K;

V=S(1);

ddt_V = (1/C)*(-g*(V-V_Nernst_K));

dSdt=[ddt_V];

end
 
