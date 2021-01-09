%file two_component_pathway.m
%bacterial two-component signalling pathway
%generates figure 6.3

function two_component_signalling

%declare parameters
global k1 ;
global km1;
global k2;
global k3;
global LT;
global RT;
global PT;
global dose_response_trigger

%trigger=0 to simulate time-series with varying ligand (Figure 6.3A)
%trigger set to 1 below to simulate multiple steady states for
%dose-response curve (Figure 6.3B)
dose_response_trigger=0;

%assign parameter values
%kinetic parameters
k1 =5;
km1 = 1;
k2=6;
k3=2;

%concentration totals
RT=2;
LT=3;
PT=8;

%assign simulation time
Tend=10;

%assign initial conditions.  The state vector S = [RL, P]
S0=[0, PT];

%run simulation
ODEFUN=@dSdt;
[t,S]=ode45(ODEFUN, [0,Tend], S0);

%determine conserved spceces for plotting
R=RT-S(:,1);
Pstar=PT-S(:,2);

%generate figure
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, Pstar, 'k', t,S(:,1), 'k-.', 'LineWidth',3)
%plot stair-case curve for ligand concentration
plot([0,0.99999,1,2.99999,3,5], [0,0,3,3,0,0], 'k:', 'LineWidth',3)
legend('Active Response Protein (P*)', 'Receptor-Ligand Complex (RL)', 'Total Ligand (L_T)')

xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')
axis([0 10 0 9])
str1(1) = {'A'};
text(-1,9,str1, 'Fontsize', 40)


%generate dose response curve

%turn off time-varying lingand profile
dose_response_trigger=1;

%generate mesh 
N=41;
LTmesh=zeros(1,N);
Pstarmesh=zeros(1,N);
RLmesh=zeros(1,N);

%run multiple simulations to steady state
Tend=12;
for j=1:N

LT=1*(j-1)/(N-1);
[t,S]=ode45(ODEFUN, [0,Tend], S0);

Pstar=PT-S(:,2);
LTmesh(j)=LT;
Pstarmesh(j)=Pstar(length(t));
RLmesh(j) = S(length(t),1);

end

%generate figure
figure(2)
hold on
set(gca,'fontsize',14)
plot(LTmesh, Pstarmesh, 'k', LTmesh, RLmesh, 'k-.', 'LineWidth',3)
xlabel('Total Ligand Concentration (L_T) (arbitrary units)')
ylabel('Steady-state concentration (arbitrary units)')
legend('Active Response Protein (P*)', 'Receptor-Ligand Complex (RL)')
axis([0 1 0 9])
str1(1) = {'B'};
text(-.1,9,str1, 'Fontsize', 40)

end

%dynamics for bacterial two-component pathway
function dS = dSdt(t,S)

global k1;
global km1;
global k2;
global k3;
global LT;
global RT;
global PT;
global dose_response_trigger

%if the trigger is zero, set a time-varying ligand profile
%(otherwise, the ligand profile will be passed from the main function)
if dose_response_trigger <1
LT=0;
if t> 1
if t < 3
LT=3;
end
end
end

%declare the state variables
RL=S(1);
P=S(2);
%determine conserved species
R=RT-RL;
Pstar=PT-P;
%assign dynamics
RLdot=k1*R*LT - km1*RL;
Pdot=-k2*RL*P+k3*Pstar;

dS=[RLdot; Pdot];

end

