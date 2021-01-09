

function Problem5_6_1

%declare model parameters
global v;
global k;
global k1; 
global k2; 
global k3;
global q;

v=10;
k=1;
k1=1; 
k2=1; 
k3=1;
q=1;



%set final time for simulation
Tend=40;
%set initial condition 
Sinit=[0 0 0];
%set right-hand-side for original model
ODEFUN=@dSdt;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], Sinit);

%generate plot
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), t, S(:,2), t, S(:,3), 'Linewidth', 3)
xlabel('Time (arbitrary units)')
end

function dS = dSdt(t,S)

global v;
global k;
global k1; 
global k2; 
global k3;
global q;

s1=S(1);
s2=S(2);
s3=S(3)

%model equations
ds1 = v/(1+(s3/k)^q) - k1*s1;
ds2 = k1*s1 -	k2*s2;
ds3 = k2*s2 -	k3*s3;


dS=[ds1 ds2 ds3]';

end