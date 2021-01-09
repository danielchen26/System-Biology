%simple model of g-protein signalling from Yi PNAS 2003
%Model of g-protein signalling pathway
%adapted from Yi et al. 2003, PNAS 100, pp. 10764-10769 
%Figure 6.20

function gprot_sine_resp

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
amp=5;

%nominal input level (Molar)
lt=10^(-9);

OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);


%Figure 6.20C 

%initial condition. State S=[RL, G, Ga]
S0=[664,9428,572];

per=2.5e1/(6.283);

%simulation time
Tend=60*per;

[t1,S1]=ode15s(@gproteinddt, [0 Tend], S0, OPTIONS);


graph_end=Tend/2;

%produce inpiut signal (for plotting)
LTseries=(lt) + (lt/amp)*cos(t1/per);

%generate Figure 6.20C
figure(1)
plot(t1(floor(1*length(t1)/2):length(t1))- t1(floor(1*length(t1)/2)), 1e9*(LTseries(floor(1*length(t1)/2):length(t1))), 'k--','Linewidth', 3)
xlabel('Time (sec)','fontsize',18)
ylabel('Ligand input (nM)','fontsize',18)
axis([0 graph_end 0.5 1.5])
set(gca,'fontsize',18)
str1(1) = {'C'};
text(-15,1.5,str1, 'Fontsize', 40)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
hold on
plot([0,1],[-1,-1], 'k--', 'Linewidth', 3)
plot(t1(floor(1*length(t1)/2):length(t1))- t1(floor(1*length(t1)/2)), S1(floor(1*length(t1)/2):length(t1),3), 'k', 'Linewidth', 3)
ylabel('G_\alpha-GTP response  (molecule per cell)','fontsize',18)
axis([0 graph_end 469 671])
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])





%Figure 6.20B

per=1e3/(6.283);
Tend=60*per;
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);
[t2,S2]=ode15s(@gproteinddt, [0 Tend], S0, OPTIONS);


graph_end=Tend/2;
LTseries=(lt) + (lt/amp)*cos(t2/per);

figure(2)

plot(t2(floor(1*length(t2)/2):length(t2))- t2(floor(1*length(t2)/2)), 1e9*(LTseries(floor(1*length(t2)/2):length(t2))), 'k--', 'Linewidth', 3)
xlabel('Time (sec)','fontsize',18)
ylabel('Ligand input (nM)','fontsize',18)
axis([0 graph_end 0.5 1.5])
set(gca,'fontsize',18)
str1(1) = {'B'};
text(-600,1.5,str1, 'Fontsize', 40)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
hold on
plot([0,1],[-1,-1], 'k--', 'Linewidth', 3)
plot(t2(floor(1*length(t2)/2):length(t2))- t2(floor(1*length(t2)/2)), S2(floor(1*length(t2)/2):length(t2),3), 'k', 'Linewidth', 3)
ylabel('G_\alpha-GTP response  (molecule per cell)','fontsize',18)
axis([0 graph_end 469 671])
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])


%Figure 6.20A

per=4e4/(6.283);
Tend=60*per;
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-9, 'refine',5);
[t3,S3]=ode15s(@gproteinddt, [0 Tend], S0, OPTIONS);


graph_end=Tend/2;
LTseries=(lt) + (lt/amp)*cos(t3/per);

figure(3)
plot(t3(floor(1*length(t3)/2):length(t3))- t3(floor(1*length(t3)/2)), 1e9*(LTseries(floor(1*length(t3)/2):length(t3))), 'k--', 'Linewidth', 3)
xlabel('Time (sec)','fontsize',18)
ylabel('Ligand input (nM)','fontsize',18)
axis([0 graph_end 0.5 1.5])
set(gca,'fontsize',18)
h1=gca;
str1(1) = {'A'};
text(-25000,1.5,str1, 'Fontsize', 40)
h2 = axes('Position',get(h1,'Position'));
hold on
plot([0,1],[-1,-1], 'k--', 'Linewidth', 3)
plot(t3(floor(1*length(t3)/2):length(t3))- t3(floor(1*length(t3)/2)), S3(floor(1*length(t3)/2):length(t3),3), 'k', 'Linewidth', 3)
ylabel('G_\alpha-GTP response (molecule per cell)','fontsize',18)
%legend('L', 'Ga')
axis([0 graph_end 469 671])
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


