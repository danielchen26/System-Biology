%file gprotein_pathway.m
%Model of g-protein signalling pathway
%adapted from Yi et al. 2003, PNAS 100, pp. 10764-10769 
%Figure 6.5

function gprotein_pathway

%declare parameters
global krl;
global krlm;
global krs;
global krd0;
global krd1;
global kg1; 
global gt;
global rt;
global kga;
global kgd1;
global kgd0;
global lt;
global per;
global input_dist;
global dB_gain
global phase_difference
global dose_response_trigger
global krs
global krd0
global krd1

%trigger=0 to simulate time-series with varying ligand (Figure 6.3A)
%trigger set to 1 below to simulate multiple steady states for
%dose-response curve (Figure 6.5B)
dose_response_trigger=0;

%parameter assignment
%kinetic parameters
krl=2*1e6;
krlm=1e-2;
kg1=1; 
kga=1e-5;
kgd1=0.11;
%concentration totals
gt=1e4;
rt=4000;
lt=0; 


%set simulation time
Tend=1200;

%set initial condition.  The state vector is S=[RL, G, Ga]
S0=[0,10000,0];


%set simulation parameters
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);
ODEFUN=@gprotddt;

%run simulation
[t,S]=ode15s(ODEFUN, [0 Tend], S0, OPTIONS);

%generate figure 6.5A
figure(1)
set(gca,'fontsize',14)
plot(t, S(:,3), 'k', t, S(:,1), 'k--', 'Linewidth', 3)
xlabel('Time (s)')
ylabel('Abundance (molecules per cell)')
legend('Active G\alpha-GTP (Ga)', 'Receptor-ligand complex (RL)')
axis([0 Tend 0 900])
str1(1) = {'A'};
text(-150,900,str1, 'Fontsize', 40)


%generate dose response curve

%turn off time-varying lingand profile
dose_response_trigger=1;


%set mesh-size
N=30;

for j=1:N+1

ltspan(j)=20*((j-1)/N)*1e-9;
lt=ltspan(j);

%set simulation length
Tend=3000;

%run simulation
[t,S]=ode15s(ODEFUN, [0,Tend], S0, OPTIONS);

%assign steady-state values
Garesp(j) = S(length(t),3);
RLresp(j) = S(length(t),1);
Gresp(j) = S(length(t),2);


end

%generate Figure 6.3B
figure(2)
set(gca,'fontsize',14)
plot(ltspan*1e9, Garesp, 'k', ltspan*1e9, RLresp, 'k-.', 'Linewidth', 3)
xlabel('Ligand concentration (nM)')
ylabel('Steady-state abundance (molecules per cell)')
legend('Active G\alpha-GTP (Ga)', 'Receptor-ligand complex (RL)')
str1(1) = {'B'};
text(-2.8,3500,str1, 'Fontsize', 40)


end


%dynamics
function dSdt = gprotddt(t,S)


global lt;
global rt;
global gt;
global per;
global input_dist;
global krl;
global krlm;
global kg1; 
global kga;
global kgd1;

global dose_response_trigger

%if the trigger is zero, set a time-varying ligand profile
%(otherwise, the ligand profile will be passed from the main function)
if dose_response_trigger < 1
    
    lt=10^(-9);
    
    if t < 100
        lt=0;
    end
    if t> 700
         lt=0;
    end

end


%locally define state variables:
rl=S(1);
g=S(2);
ga=S(3);

%conservations:

r=rt-rl;
gbg=gt-g;
gd=gt-g-ga;



%concentration kinetics:

ddt_rl= krl*lt*r - krlm*rl;
ddt_g= -kga*rl*g + kg1*gd*gbg;
ddt_ga= kga*rl*g - kgd1*ga;

dSdt=[ddt_rl; ddt_g; ddt_ga];

end