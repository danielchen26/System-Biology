%file chemotaxis.m
%Model of E. coli chemotaxis signalling pathway
%Figure 6.14

function chemotaxis_model

%declare parameters
global k1 ;
global k2;
global k3;
global km1 ;
global km2;
global km3;
global k5;
global k4;
global km5 ;
global km4;
global KM1;
global KM2;
global cheR;
global L;

%assign parameter values
k1=200;
k2=1;
k3=1; 
km1=1;
km2=1;
km3=1; 
k5=0.05;
km5=0.005;
k4=1;
km4=1;
KM1=1;
KM2=1;
cheR=1;
L=20;

%set simlatuion options
Tend=50;
OPTIONS = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
ODEFUN=@dSdt;


%set up ligand profile for plotting
lt=[0 10 10 30 30 50];
l=0.01+[0.002 0.002 0.004 0.004 0.008 0.008];


%assign initial condition species vector S=[Am AmL A AL B BP]
%these values were determined by running a previous simulation to steady
%state with L=20.
S0=[0.0360    1.5593    0.0595    0.3504    0.7356    0.2644];
[t,S]=ode45(ODEFUN, [0,Tend], S0, OPTIONS);



%generate figure
figure(1)
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', lt,l, 'k--', 'Linewidth', 3)
xlabel('Time (arbitrary units)')
ylabel('Concentration of active CheA [Am] (arbitrary units)')
axis([0 50 0.01 0.04]);
h1=gca;
 

end

%dynamics
function dS = dSdt(t,S)


global k1 ;
global k2;
global k3;
global km1 ;
global km2;
global km3;
global k5;
global k4;
global km5 ;
global km4;
global KM1;
global KM2;
global cheR;
global L;

%assign state variables
Am=S(1);
AmL=S(2);
A=S(3);
AL=S(4);
B=S(5);
BP=S(6);

%set up time-varying ligan profile
if (t>10)
L=40;
end
if (t>30)
L=80;
end


%dynamics
 dS=[km1*cheR - (k1*BP*Am)/(KM1 + Am) - k3*Am*L + km3*AmL ;
    km2*cheR - (k2*BP*AmL)/(KM2 + AmL) + k3*Am*L - km3*AmL;
    -km1*cheR + (k1*BP*Am)/(KM1 + Am) - k4*A*L + km4*AL;
    -km2*cheR + (k2*BP*AmL)/(KM2 + AmL) + k4*A*L - km4*AL;
    -(k5*Am*B) +  (km5*BP);
     (k5*Am*B) -  (km5*BP)];
 
end

