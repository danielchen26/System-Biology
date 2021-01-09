

function Problem5_6_1

%declare model parameters
global e0;
global e1;
global e2; 
global S0; 
global km1;
global km2;
global km3;
global km4; 
global km5;
global v1; 
global v2; 
global v3; 
global v4; 
global v5;

e0=1;
e1=1.5;
e2=2; 
S0=1; 
km1=1;
km2=1;
km3=1;
km4=1; 
km5=1;
v1=1; 
v2=1; 
v3=2; 
v4=1; 
v5=0.5;


%set final time for simulation
Tend=10;
%set initial condition 
Sinit=[0 0];
%set right-hand-side for original model
ODEFUN=@dSdt;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], Sinit);

%generate plot
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), t, S(:,2), 'Linewidth', 3)
xlabel('Time (arbitrary units)')
end

function dS = dSdt(t,S)

global e0;
global e1;
global e2; 
global S0; 
global km1;
global km2;
global km3;
global km4; 
global km5;
global v1; 
global v2; 
global v3; 
global v4; 
global v5;

s1=S(1);
s2=S(2);

%model equations
ds1 = e0*(v1*S0-v2*s1)/(1+S0/km1+s1/km2) - e1*(v3*s1-v4*s2)/(1+s1/km3+s2/km4);
ds2 = e1*(v3*s1-v4*s2)/(1+s1/km3+s2/km4) - e2*(v5*s2)/(1+s2/km5);

dS=[ds1 ds2]';

end