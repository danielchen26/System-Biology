%file Turing_pattern.m
%model of pattern formation by morphogen diffusion 
%problem 8.6.11

function Turing_pattern
m = 0;

%set up solution mesh
N=80;
L=20;
x = linspace(-L,L,N);
t = linspace(0,50,50);

%solve equation for Problem 8.6.11(a)
sol = pdepe(m,@turing_eqn_a,@turing_initial,@turing_bc,x,t);
u = sol(:,:,1);
v = sol(:,:,2);

%plot solution
figure(1)
set(gca,'fontsize',14)
surf(x,t,u) 

figure(2)
set(gca,'fontsize',14)
surf(x,t,v)

figure(3)
set(gca,'fontsize',14)
hold on
plot(x, u(length(t),:),x, v(length(t),:), 'linewidth', 2);

%solve equation for Problem 8.6.11(b)
sol = pdepe(m,@turing_eqn_b,@turing_initial,@turing_bc,x,t);
u = sol(:,:,1);
v = sol(:,:,2);

%plot solution
figure(4)
set(gca,'fontsize',14)
surf(x,t,u) 

figure(5)
set(gca,'fontsize',14)
surf(x,t,v)

figure(6)
set(gca,'fontsize',14)
hold on
plot(x, u(length(t),:),x, v(length(t),:), 'linewidth', 2);

   
end

%set equation for problem 8.6.11(a)
function [c,b,s] = turing_eqn_a(x,t,u,DuDx)

%activator-inhibitor
eps=0.1;
Df=10;
alpha=1;
beta=1;
gamma=1;
delta=1.1;


c = [1;1];
b = [eps*Df;Df].*DuDx;
s = [alpha*u(1)^2/u(2)-beta*u(1); gamma*u(1)^2-delta*u(2)];

end

%set equation for problem 8.6.11(b)
function [c,b,s] = turing_eqn_b(x,t,u,DuDx)

%substrate-depletion
alpha=0.9;
gamma=alpha;
beta=0.9;
D=6;
eps=0.075;
r1=0.2;
r2=0.02;

u1=u(1);
u2=u(2);
c = [1;1];
b = [eps*D;D].*DuDx;
s = [alpha*u1+u2-r1*u1*u2-r2*u1*u2^2; -gamma*u1-beta*u2+r1*u1*u2+r2*u1*u2^2];
end

%set initial condition
function value = turing_initial(x)

value = [1; 1/(x^2+1)];
end

%set boundary condition
function [pl,ql,pr,qr] = turing_bc(xl,ul,xr,ur,t)

pl = [0; 0];
ql = [1; 1];
pr = [0;0];
qr = [1;1];
end

