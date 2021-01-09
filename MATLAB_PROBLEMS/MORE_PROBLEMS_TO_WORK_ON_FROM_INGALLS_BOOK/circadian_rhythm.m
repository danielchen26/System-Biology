%file circadian_rhythm.m
%Goldbeter model of circadian rhythm in Drosophila
%reviewed in Goldbeter (1996) Biochemical oscillations and cellular
%rhythms, Cambridge University Press
%Figure 7.19


function circadian_rhythm

%declare model parameters
global vs ;
global vm;
global vd;
global ks ;
global kt1;
global kt2;
global v1;
global v2;
global v3 ;
global v4;
global k1 ;
global k2 ;
global k3 ;
global k4 ;
global ki ;
global km1;
global kd;
global n;

%assign parameter values
vs=0.76;
vm=0.65;
vd=0.95;
ks=0.38;
kt1=1.9;
kt2=1.3;
v1=3.2;
v2=1.58;
v3=5;
v4=2.5;
k1=1;
k2=1;
k3=2;
k4=2;
ki=1;
km1=0.5;
kd=0.2;
n=4;





%set simulation parameters
ODEFUN=@circadianddt;
Tend=200;

%set initial condition: state = [M, P0, P1, P2, PN]
S0=[1 1 0 0 0 ];

%run simulation (run from before time zero to the limit cycle by zero
[t,S]=ode45(ODEFUN, [-50,Tend], S0);

%generate figure 7.19A
figure(1)
set(gca, 'fontsize', 14)
plot(t, S(:,1), 'k', t, S(:,2)+ S(:,3) + S(:,4) + S(:,5), 'k--', t, S(:,5), 'k-.', 'LineWidth', 3)
axis([0 72 0 6])
legend('mRNA (M)', 'total PER (P_T)', 'nuclear PER (P_N)')
str1(1) = {'A'};
text(-8,6,str1, 'Fontsize', 40)
xlabel('Time (h)')
ylabel('Concentration (\muM)')
    
    
%generate continuation diagram Figure 7.17B

Tend=1000;
%choose a mesh of vd values
N=15;

peakepsilon=0.05;
vd_values=zeros(N,1);
period=zeros(N,1);

%over the mesh, determine the period of oscillation
for i=1:N
    vd=0.5 + (2.55 - 0.5)*(i-1)/(N-1);
    vd_values(i)=vd;
    [t,s]=ode45(ODEFUN, [-500,Tend], S0);
    
figure(99)
%determine peak
inds=peakdetect(s(:,1));
maxpeak=max(s(floor(0.8*length(t)):length(t),1));


%determine peak-to-peak time interval
flag1=0; 
n1=0;
while flag1==0
    index_of_last_peak=inds(length(inds)-n1)
    if maxpeak-s(index_of_last_peak,1)<peakepsilon
        flag1=1;
    else
        n1=n1+1;
    end
end
    
flag2=0;
n2=n1+1;
while flag2==0
    index_of_secondlastpeak=inds(length(inds)-n2)
    if maxpeak-s(index_of_secondlastpeak,1)<peakepsilon
        flag2=1;
    else
        n2=n2+1;
    end
end    
    
period(i)=t(index_of_last_peak)-t(index_of_secondlastpeak);

end

%produce figure 7.19B
figure(2)
set(gca, 'fontsize', 14)
plot(vd_values, period, 'k', 'Linewidth', 3)
xlabel('maximal degradation rate of PER (v_d) (\muM/h)')
ylabel('Period of oscillations (h)')
axis([0.4 2.7 10 70])
str1(1) = {'B'};
text(0.1,70,str1, 'Fontsize', 40)
   
end

%model dynamics
function dS = circadianddt(t,S)

global vs ;
global vm;
global vd;
global ks ;
global kt1;
global kt2;
global v1;
global v2;
global v3 ;
global v4;
global k1 ;
global k2 ;
global k3 ;
global k4 ;
global ki ;
global km1;
global kd;
global n;


M=S(1);
P0=S(2);
P1=S(3);
P2=S(4);
PN=S(5);


dS=[vs/(1+(PN/ki)^n) - (vm*M)/(km1+M);ks*M - (v1*P0)/(k1+P0) + (v2*P1)/(k2 + P1);(v1*P0)/(k1+P0) - (v2*P1)/(k2 + P1) - (v3*P1)/(k3+P1) + (v4*P2)/(k4+P2);(v3*P1)/(k3+P1) - (v4*P2)/(k4+P2)  - kt1*P2 + kt2*PN - (vd*P2)/(kd+P2);kt1*P2 - kt2*PN
];

end