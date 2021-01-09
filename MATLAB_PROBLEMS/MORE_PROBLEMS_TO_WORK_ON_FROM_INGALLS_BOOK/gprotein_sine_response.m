%simple model of g-protein signalling from Yi PNAS 2003
%Model of g-protein signalling pathway
%adapted from Yi et al. 2003, PNAS 100, pp. 10764-10769 
%Figure 6.19

function gprotein_sine_response

%parameter declaration
global krl;
global krlm;
global kg1; 
global gt;
global rt;
global kga;
global kgd1;
global lt;
global per;
global amp;

%parameter assignment

%kinetic parameters
krl=2*1e6;
krlm=1e-2;
kg1=1; 
kga=1e-5;
kgd1=0.11; 

%fixed concentrations
gt=1e4;
rt=4000;

%input signal parameters
per=200/6.283;
amp=5;

%nominal input level (Molar)
lt=10^(-9);


%simulation time 
Tend=1000;

%initial condition. State S=[RL, G, Ga]
S0=[0,10000,0];

%run simulation
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);
ODEFUN=@gproteinddt;
[t,S]=ode15s(ODEFUN, [0 Tend], S0, OPTIONS);

%define time-varying input signal (for plotting)
ltcurrent=lt + (lt/amp)*cos(t/per);

%generate figure 6.19
plot(t, 1e9*ltcurrent, 'k--' ,'Linewidth', 3)
xlabel('Time (sec)','fontsize',14)
ylabel('Ligand (input) (nM)','fontsize',14)
axis([0 1000 0 2.2])
set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
hold on
plot([0,1],[-1,-1], 'k--', 'Linewidth', 3)
plot(t, S(:,3), 'k', 'Linewidth', 3)
ylabel('G_a-GTP (output) (molecules per cell)','fontsize',14)
set(gca,'fontsize',14)
legend('Ligand (input)', 'G_a-GTP (response)')
axis([0 1000 0 700])
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])


end



%dynamics
function dSdt = gproteinddt(t,S)


global lt;
global rt;
global gt;
global per;
global amp;
global krl;
global krlm;
global kg1; 
global kga;
global kgd1;


%define state variables:
rl=S(1);
g=S(2);
ga=S(3);

%conservations:
r=rt-rl;
gbg=gt-g;
gd=gt-g-ga;

%determine current value of input signal (as a function of time t)
ltcurrent=lt + (lt/amp)*cos(t/per);

%dynamics
ddt_rl= krl*ltcurrent*r - krlm*rl;
ddt_g= -kga*rl*g + kg1*gd*gbg;
ddt_ga= kga*rl*g - kgd1*ga;

dSdt=[ddt_rl; ddt_g; ddt_ga];

end
