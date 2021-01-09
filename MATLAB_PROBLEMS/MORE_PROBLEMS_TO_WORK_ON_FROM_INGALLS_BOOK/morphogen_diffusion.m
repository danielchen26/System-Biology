%file morphogen_diffusion.m
%model of morphogen diffusion in a one-dimensional domain.
%problem 8.6.10

function morphogen_diffusion


m = 0;

%set up solution mesh
N=30;
L=10;
x = linspace(0,L,N);
t = linspace(0,20,20);

%solve PDE
sol = pdepe(m,@morphogen_eqn,@morphogen_initial,@morphogen_bc,x,t);
u = sol(:,:,1);

%plot
figure(1)
set(gca,'fontsize',14)
surf(x,t,u)    
xlabel('Position (cm)')
ylabel('Time (msec)')
zlabel('concentration')

figure(2)
set(gca,'fontsize',14)
hold on
plot(x, u(1,:),'k', 'linewidth', 2)
plot(x, u(,:),'g--', 'linewidth', 2)
plot(x, u(15,:),'r', 'linewidth', 2)
xlabel('Position (cm)')
    
end

%set equation
function [c,b,s] = morphogen_eqn(x,t,u,DuDx)

c = 1;
b = 10*DuDx;
s = -u;
end

%set initial condition
function value = morphogen_initial(x)

value = [0];
end

%set boundary conditions
function [pl,ql,pr,qr] = morphogen_bc(xl,ul,xr,ur,t)

a=1/10;

pl = [a];
ql = [1];

pr = [0];
qr = [1];
end
