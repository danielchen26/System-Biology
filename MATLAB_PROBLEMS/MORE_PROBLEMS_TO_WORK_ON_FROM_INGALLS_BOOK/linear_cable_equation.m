%file linear_cable_equation.m
%Model of excitation over a passive membrane
%Figure 8.10

function linear_cable_equation

%declare model parameters
global E;
global Cap;
global gm;
global ga;
global r;

%assign parameter values
E=-93.6;%mV
Cap=0.98; %microF/cm^2
gm=0.0144; %mS/cm^2
ga= 30;  % mSiemans/cm
r=0.01; %cm

%set up solution mesh
m = 0;
x = linspace(-0.5,0.5,150);
t = linspace(0,1,100);

%solve partial differential equation
sol = pdepe(m,@linear_cable_pde_eqn,@linear_cable_pde_initial,@linear_cable_pde_bc,x,t);

%declare solution
u = sol(:,:,1);
    
%produce Figure8.11A
figure(1)
set(gca,'fontsize',14)
surf(10*x,t,u)    
xlabel('Position (mm)')
ylabel('Time (msec)')
zlabel('Membrane voltage (mV)')
colormap hsv
colormap(gray)
shading interp

%produce Figure 8.11B
figure(2)
hold on
set(gca,'fontsize',14)
plot(10*x, u(1,:), 'k', 'LineWidth', 2);
plot(10*x, u(5,:),'k--', 'LineWidth', 2);
plot(10*x, u(20, :),'k:', 'LineWidth', 2);
xlabel('Position (mm)')
ylabel('Membrane voltage (mV)')
legend('t = 0 msec', 't = 0.1 msec', 't=0.2 msec')
str1(1) = {'B'};
text(-6.4,-60,str1, 'Fontsize', 40)
end

%set equation
function [c,b,s] = linear_cable_pde_eqn(x,t,u,DuDx)
global E;
global Cap;
global gm;
global ga;
global r;

c = Cap;
b = (ga*r/2)*DuDx;
s = gm*(E-u);
end

%set initial condition
function value = linear_cable_pde_initial(x)
global E;
value = E+30/(1+10*(10*x)^2);
end

%set boundary conditions
function [pl,ql,pr,qr] = linear_cable_pde_bc(xl,ul,xr,ur,t)

global E;
pl = ul-E;
ql = 0;
pr = ur-E;
qr = 0;
end

