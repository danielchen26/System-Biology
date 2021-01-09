%file Goodwin.m
%Goodwin oscillator model
%from Goodwin (1965) Adv. in. Enz. Reg. 3 pp. 425-428.
%Figure 7.17

function Goodwin

clear all

%declare model parameters
global a1
global kappa1
global k1
global b1
global alpha1
global beta1
global gamma1
global delta1
global n

%assign parameter values
a1=360;
kappa1=43;
k1=1;
b1=1;
alpha1=1;
beta1=0.6;
gamma1=1;
delta1=0.8;
n=12


%set simulation parameters
Tend=35
ODEFUN=@goodwinddt;
options=odeset('Refine', 6);

%set initial condition: state = [X, Y, Z]
x0=[0 0 0]';
%run simulation
[t,s]=ode45(ODEFUN, [0,Tend], x0, options);

%generate figure 7.17A
figure(1)
set(gca,'fontsize',14)
plot(t,s(:,1), 'k',t,s(:,2), 'k--',t, s(:,3), 'k:','Linewidth', 3)
 xlabel('Time (arbitrary units)')
 ylabel('Concentration (arbitrary units)')
 legend('mRNA (X)', 'enzyme (Y)', 'metabolite (Z)')
str1(1) = {'A'};
text(-3.5,6,str1, 'Fontsize', 40)

%run simulations for figure 7.17B

 [t1,s1] = ode15s(ODEFUN,[0,Tend], [0,0,0]);
 [t2,s2] = ode15s(ODEFUN,[0,Tend], [0,2.5,0]);
 [t3,s3] = ode15s(ODEFUN,[0,Tend], [0,5,0]);

%generate figure 7.17B
figure(2)
set(gca, 'fontsize', 14)
 plot3(s1(:,1), s1(:,2), s1(:,3),'k--', 'Linewidth', 2)
 hold on
 plot3(s2(:,1), s2(:,2), s2(:,3),'k--', 'Linewidth', 2)
 plot3(s3(:,1), s3(:,2), s3(:,3),'k--', 'Linewidth', 2)
 xlabel('[X] (arbitrary units)')
 ylabel('[Y] (arbitrary units)')
 zlabel('[Z] (arbitrary units)')
 str1(1) = {'B'};
text(-3.5, -3.5,5,str1, 'Fontsize', 40)
 
 


end


%dynamics of Goodwin model
function dS=goodwinddt(t,x)

global a1
global kappa1
global k1
global b1
global alpha1
global beta1
global gamma1
global delta1
global n


X=x(1);
Y=x(2);
Z=x(3);

 

dS =[a1/(kappa1+k1*Z^n) - b1*X,  alpha1*X - beta1*Y, gamma1*Y - delta1*Z]';

end

