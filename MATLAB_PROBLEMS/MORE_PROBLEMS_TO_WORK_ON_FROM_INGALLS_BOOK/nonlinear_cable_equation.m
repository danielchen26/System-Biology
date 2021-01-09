%file nonlinear_cable_equation.m
%Model of excitation over an excitable membrane
%Figure 8.11

function nonlinear_cable_equation

global ga;
global r;

ga= 30;  %3 Siemans/meter
r=0.05; %cm, =0.01mm: squid giant axon is 0.5mm = 0.05cm

%set up solution mesh
m = 0;
N=300;
x = linspace(-30,30,N);
t = linspace(0,200,200);

%solve the partial differential equations
sol = pdepe(m,@nonlinear_cable_pde_eqn,@nonlinear_cable_pde_initial,@nonlinear_cable_pde_bc,x,t);

V = sol(:,:,1);
w = sol(:,:,2);

%produce Figure8.12A
figure(1)
set(gca,'fontsize',14)
surf(x(floor((N+1)/2):floor(3*N/4)),t,V(:,floor((N+1)/2):floor(3*N/4)), 'EdgeColor','none')    
xlabel('Position (cm)')
ylabel('Time (msec)')
zlabel('Membrane Voltage (mV)')
colormap hsv
colormap(gray)
shading interp

%produce Figure 8.12B
figure(2)
set(gca,'fontsize',14)
hold on
plot(x, V(1,:),'k', 'linewidth', 2);
plot(x, V(50,:),'k--', 'linewidth', 2);
plot(x, V(150,:),'k:', 'linewidth', 2);
axis([0 15 -80 40])
xlabel('Position (cm)')
ylabel('Membrane voltage (mV)')
legend('t=0 msec', 't=50 msec', 't=150 msec')
str1(1) = {'B'};
text(-2.5,60,str1, 'Fontsize', 40)

end

%set equations
function [c,b,s] = nonlinear_cable_pde_eqn(x,t,u,DuDx)

global ga;
global r;

Cap=20;
V_Nernst_K=-84;
gbar_K=8;
V_Nernst_Ca=120;
gbar_Ca=4.4;
gbar_leak=0.5;
V_Nernst_leak=-60;
v1=-1.2;
v2=18;
v3=2;
v4=30;
phi=0.04;

m_inf = @(x) ( 0.5*(1+tanh((x-v1)/v2)) );
w_inf = @(x) ( 0.5*(1+tanh((x-v3)/v4)) );
tau_w=0.8;

c = [1;1];
b = [(ga*r/(2*Cap));0].*DuDx;
s = [(1/Cap)*(-gbar_Ca*m_inf(u(1))*(u(1) -V_Nernst_Ca) - gbar_K*u(2)*(u(1) - V_Nernst_K)- gbar_leak*(u(1) - V_Nernst_leak)); phi*(w_inf(u(1))-u(2))/tau_w];
end

%set initial condition
function value = nonlinear_cable_pde_initial(x)

value = [60*(1/(1+1*x^8))-61.85746; 0.004247517];
end

%set boundary conditions
function [pl,ql,pr,qr] = nonlinear_cable_pde_bc(xl,ul,xr,ur,t)

pl = [0; 0];
ql = [1; 1];
pr = [0;0];
qr = [1;1];

end