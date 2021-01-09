%file asymmetric_network.m
%Model of asymmetric network from Figure 4.1. This code generates phase plane in Figures 4.2, 4.3, 4.4A, 4.5 and 
%continuation diagram in Figure 4.18

function asymmetric_network

%declare model parameters
global k1;
global k2;
global k3;
global k4;
global k5;
global n;


%assign values to model parameters
k1=20;
k2=5;
k3=5;
k4=5;
k5=2;
n=4;


%declare dynamics and set simulation options
ODEFUN=@asyddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);

%set simulation end time
Tend=1.5;

%simulate a set of trajectories from several initial conditions
[t1,S1]=ode45(ODEFUN, [0,Tend], [0,0], options);
[t2,S2]=ode45(ODEFUN, [0,Tend], [0.5,0.6], options);
[t3,S3]=ode45(ODEFUN, [0,Tend], [0.17,1.1], options);
[t4,S4]=ode45(ODEFUN, [0,Tend], [0.25,1.9], options);
[t5,S5]=ode45(ODEFUN, [0,Tend], [1.85,1.7], options);

%generate Figure 4.2A
figure(1)
set(gca,'fontsize',14)
plot(t1, S1(:,1), 'k', t1, S1(:,2), 'k--', 'LineWidth',3)
axis([0 1.5 0 2])
legend('S_1', 'S_2')
xlabel('Time')
    ylabel('Concentration')
         str1(1) = {'A'};
text(-.18,2,str1, 'Fontsize', 40)

%generate Figure 4.2B
figure(2)
set(gca,'fontsize',14)
plot(S1(:,1),S1(:,2),'k','LineWidth',3)
axis([0 2 0 2])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')
        str1(1) = {'B'};
text(-.25,2,str1, 'Fontsize', 40)

%generate Figure 4.3B
figure(3)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2),'k','LineWidth',3)
plot(S2(:,1),S2(:,2),'k','LineWidth',3)
plot(S3(:,1),S3(:,2),'k','LineWidth',3)
plot(S4(:,1),S4(:,2),'k','LineWidth',3)
plot(S5(:,1),S5(:,2),'k','LineWidth',3)
hold off         
str1(1) = {'B'};
text(-.2,2,str1, 'Fontsize', 40)

axis([0 2 0 2])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')

%generate Figure 4.3A
figure(4)
set(gca,'fontsize',14)
hold on
plot(t1, S1(:,1),  'k',t1, S1(:,2), 'k','LineWidth',3)
plot(t2, S2(:,1), 'kx', t2, S2(:,2), 'kx', 'LineWidth',3)
plot(t3, S3(:,1), 'k.-', t3, S3(:,2), 'k.-', 'LineWidth',3)
plot(t4, S4(:,1), 'k--', t4, S4(:,2), 'k--', 'LineWidth',3)
plot(t5, S5(:,1), 'k:', t5, S5(:,2), 'k:', 'LineWidth',3)
hold off
axis([0 1.5 0 2])
xlabel('Time')
ylabel('Concentration')
str1(1) = {'A'};
text(-.18,2,str1, 'Fontsize', 40)

    
%generate Figure 4.4A
figure(5)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2),'k','LineWidth',3)
plot(S2(:,1),S2(:,2),'k','LineWidth',3)
plot(S3(:,1),S3(:,2),'k','LineWidth',3)
plot(S4(:,1),S4(:,2),'k','LineWidth',3)
plot(S5(:,1),S5(:,2),'k','LineWidth',3)
%generate direction field
[xx,yy]=meshgrid(0:0.1:2, 0:0.1:2)
xdot=k1./(1+yy.^n) - k3*xx - k5*xx;
ydot=k2 - k4*yy + k5*xx;
L = sqrt(xdot.^2 + ydot.^2); % vector lengths
quiver(xx,yy,xdot./L,ydot./L,0.5, 'Color', 'black');
axis([0 2 0 2])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')
str1(1) = {'A'};
text(-.22,2,str1, 'Fontsize', 40)


%generate Figure 4.5A
figure(6)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2),'k','LineWidth',3)
plot(S2(:,1),S2(:,2),'k','LineWidth',3)
plot(S3(:,1),S3(:,2),'k','LineWidth',3)
plot(S4(:,1),S4(:,2),'k','LineWidth',3)
plot(S5(:,1),S5(:,2),'k','LineWidth',3)
%s1 nullcline
ns12=linspace(0,2,100);
ns11=k1./((k3+k5)*(1.+ns12.^n));
ns21=linspace(0,2,100);
ns22=(k2+k5*ns21)/k4;
%plot nullclines
plot(ns11, ns12, 'k--', 'LineWidth',2)
plot(ns21, ns22, 'k--', 'LineWidth',2)
axis([0 2 0 2])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')
str1(1) = {'A'};
text(-.22,2,str1, 'Fontsize', 40)



%generate Figure4.5B
figure(7)
set(gca,'fontsize',14)
hold on
%generate direction field
[xx,yy]=meshgrid(0:0.1:2, 0:0.1:2)
xdot=k1./(1+yy.^n) - k3*xx - k5*xx;
ydot=k2 - k4*yy + k5*xx;
%quiver(xx,yy,xdot,ydot,0.9);
L = sqrt(xdot.^2 + ydot.^2); % vector lengths
quiver(xx,yy,xdot./L,ydot./L,0.5, 'Color', 'black');
%s1 nullcline
ns12=linspace(0,2,100);
ns11=k1./((k3+k5)*(1.+ns12.^n));
ns21=linspace(0,2,100);
ns22=(k2+k5*ns21)/k4;
%plot nullclines
plot(ns11, ns12, 'k--', 'LineWidth',2)
plot(ns21, ns22, 'k--', 'LineWidth',2)
axis([0 2 0 2])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')
str1(1) = {'B'};
text(-.22,2,str1, 'Fontsize', 40)


%generate continuation diagram Figure 4.18
%Run N simulations to steady state for different values of k_1
N=30
k1_values=zeros(1,N);
S1_state=zeros(1,N);
for i=1:N
    k1=0+(40)*((i-1)/(N-1));
    [t,S]=ode45(ODEFUN, [0,100], [0,0], options);
    k1_values(i)=k1;
    S1_state(i)=S(length(t),1);
end
figure(8)
set(gca,'fontsize',14)
plot(k1_values,S1_state, 'k', 'linewidth', 3)
xlabel('k_1')
ylabel('Steady state S_1 concentration')

    
end

%declare right-hand-side
function dS = asyddt(t,S)

global k1 ;
global k2;
global k3;
global k4;
global k5;
global n;


dS=[k1/(1+S(2)^n) - k3*S(1) - k5*S(1); k2 - k4*S(2) + k5*S(1)];

    
end