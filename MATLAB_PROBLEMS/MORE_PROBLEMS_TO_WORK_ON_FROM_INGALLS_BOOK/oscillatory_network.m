%file oscillatory_network.m
%Model of oscillatory network from Figure 4.14. This code generates Figures
%4.15, 4.16, and 4.17

function oscillatory_network

%declare model parameters
global k0;
global k1;
global k2;
global n;

%declare dynamics and set simulation options
ODEFUN=@oscddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);

%Set simulation length
Tend=8;

%Two sets of parametrizations are considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parametrization for damped oscillations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0=8;
k1=1;
k2=5;
n=2;  

%generate simulation traces
[t1,S1]=ode45(ODEFUN, [0,Tend], [1.5,1], options);
[t2,S2]=ode45(ODEFUN, [0,Tend], [0,1], options);
[t3,S3]=ode45(ODEFUN, [0,Tend], [0,3], options);
[t5,S4]=ode45(ODEFUN, [0,Tend], [2,0], options);

%figure 4.15A
figure(1)
set(gca,'fontsize',14)
plot(t1, S1(:,1), 'k', t1, S1(:,2), 'k--', 'LineWidth',3)
axis([0 8 0 3.5])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')
str1(1) = {'A'};
text(-1,3.5,str1, 'Fontsize', 40)

%figure 4.15B
figure(2)
set(gca,'fontsize',14)
hold on
plot(S2(:,1),S2(:,2),'k','LineWidth',2);
plot(S3(:,1),S3(:,2),'k','LineWidth',2);
plot(S4(:,1),S4(:,2),'k','LineWidth',2);

%direction field
[xx,yy]=meshgrid(0:0.4:4, 0:0.4:4)
xdot=k0 - k1*xx.*(1+yy.^n)
ydot=k1*xx.*(1+yy.^n) - k2*yy
L = sqrt(xdot.^2 + ydot.^2); % vector lengths
quiver(xx,yy,xdot./L,ydot./L,0.5, 'Color', 'black');

%s1 nullcline
ns12=linspace(0,4,100);
ns11=k0./(k1*(1+ns12.^n));
%s2 nullcline
ns22=linspace(0,4,100);
ns21=k2*ns22./(k1*(1+ns22.^n));

plot(ns11, ns12, 'k--', 'LineWidth',3)
plot(ns21, ns22, 'k--', 'LineWidth',3)

axis([0 4 0 4])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')
str1(1) = {'B'};
text(-0.5,4,str1, 'Fontsize', 40)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parametrization for sustained (limit cycle) oscillations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
k0=8;
k1=1;
k2=5;
n=2.5;


Tend=40;

[t1,S1]=ode45(ODEFUN, [0,Tend], [0,1], options);
[t2,S2]=ode45(ODEFUN, [0,Tend], [0,3], options);
[t3,S3]=ode45(ODEFUN, [0,Tend], [2,0], options);


%figure 4.16A
figure(3)
set(gca,'fontsize',14)
plot(t1, S1(:,1), 'k', t1, S1(:,2), 'k--', 'LineWidth',3)
axis([0 8 0 4])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')
str1(1) = {'A'};
text(-1,4,str1, 'Fontsize', 40)

%figure4.16B
figure(4)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2),'k','LineWidth',2)
plot(S2(:,1),S2(:,2),'k','LineWidth',2)
plot(S3(:,1),S3(:,2),'k','LineWidth',2)

%direction field
[xx,yy]=meshgrid(0:0.4:4, 0:0.4:4)
xdot=k0 - k1*xx.*(1+yy.^n)
ydot=k1*xx.*(1+yy.^n) - k2*yy
L = sqrt(xdot.^2 + ydot.^2); % vector lengths
quiver(xx,yy,xdot./L,ydot./L,0.5, 'Color', 'black');

%s1 nullcline
ns12=linspace(0,4,100);
ns11=k0./(k1*(1+ns12.^n));
%s2 nullcline
ns22=linspace(0,4,100);
ns21=k2*ns22./(k1*(1+ns22.^n));

plot(ns11, ns12, 'k--', 'LineWidth',3)
plot(ns21, ns22, 'k--', 'LineWidth',3)

axis([0 4 0 4])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')
str1(1) = {'B'};
text(-0.5,4,str1, 'Fontsize', 40)

hold off

%figure 4.17
figure(5)
set(gca,'fontsize',14)
hold on

[t1,S1]=ode45(ODEFUN, [0 7.95], [1.8,1.6], options);
[t2,S2]=ode45(ODEFUN, [0 100], [1.8,1.6], options);

%plot spiral trajectory
plot(S1(:,1),S1(:,2),'k', 'LineWidth',2)

%plot limit cycle; reached by halfway point of long simulation
plot(S2(length(t2)/2:length(t2),1),S2(length(t2)/2:length(t2),2),'k','LineWidth',3)

%direction field
[xx,yy]=meshgrid(0:0.1:4, 0:0.1:3)
xdot=k0 - k1*xx.*(1+yy.^n) 
ydot=k1*xx.*(1+yy.^n) - k2*yy
L = sqrt(xdot.^2 + ydot.^2); % vector lengths
quiver(xx,yy,xdot./L,ydot./L,0.5, 'Color', 'black');

%s1 nullcline
ns12=linspace(0,3,100);
ns11=k0./(k1*(1+ns12.^n));
%s2 nullcline
ns22=linspace(0,3,100);
ns21=k2*ns22./(k1*(1+ns22.^n));

plot(ns11, ns12, 'k--', 'LineWidth',3)
plot(ns21, ns22, 'k--', 'LineWidth',3)

axis([1 2.75 1.1 2.3])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')

hold off


end

%dynamics for oscillatory network
function dS = oscddt(t,S)

global k0;
global k1;
global k2;
global n;


dS=[k0 - k1*S(1)*(1+S(2)^n);k1*S(1)*(1+S(2)^n) - k2*S(2)];

end
