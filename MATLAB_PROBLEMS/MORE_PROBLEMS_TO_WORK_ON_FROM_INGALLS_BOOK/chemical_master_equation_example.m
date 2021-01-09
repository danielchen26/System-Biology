%file chemical_master_equation_example
% solutions of the chemical master equation for the simple system
% A --> B, propensity k1*N_A
% B --> A, propensity k2*N_B

%k1=1, k2=3.  Maple gives the matrix exponential as:

%initial condition for the three possible states: 
%(N_A, N_B) = (2,0), (1,1), (0,2)
IC=[1/3 1/3 1/3]';

%declare vectors to store solution values
X=zeros(3,100);
A=zeros(3,3,100);

%declare kinetic parameters and system size
N=2
k1=3
k2=1

%generate the 'A' matrix for the linear ODE system
A=zeros(N+1,N+1);
for i=1:N+1
    A(i,i)=-(N-(i-1))*k1-(i-1)*k2;
end

for i=1:N
    A(i,i+1)=i*k2;
    A(i+1,i)=(N-(i-1))*k1;
end


%produce figure 7.40
%solutions are stated explicitly in terms of the matrix exponential of the
%'A' matrix
xd={'N=(2,0)','N=(1,1)','N=(0,2)'};
figure(1)
%panel A
subplot(1,3,1)
set(gca,'fontsize',14)
bar(expm(0*A)*IC, 'k')
set(gca,'ylim',[0,0.6],'xticklabel',xd);
axis([0 4 0 0.6])
ylabel('Probability')
%panel B
subplot(1,3,2)
set(gca,'fontsize',14)
bar(expm(0.1*A)*IC, 'k')
set(gca,'ylim',[0,0.6],'xticklabel',xd);
axis([0 4 0 0.6])
ylabel('Probability')
%panel C
subplot(1,3,3)
set(gca,'fontsize',14)
bar(expm(1*A)*IC, 'k')
set(gca,'ylim',[0,0.6],'xticklabel',xd);
axis([0 4 0 0.6])
ylabel('Probability')


%Figure 7.41
%repeat for different system sizes (i.e. values of N)
figure(2)
hold on

N=2

A1=zeros(N+1,N+1);

for i=1:N+1
    A1(i,i)=-(N-(i-1))*k1-(i-1)*k2;
end

for i=1:N
    A1(i,i+1)=i*k2;
    A1(i+1,i)=(N-(i-1))*k1;
end

IC=zeros(N+1,1);

IC(1)=1;

%produce Figure 7.41A
subplot(1,3,1)
xd={'0','1','2'};
set(gca,'fontsize',14)
bar(expm(10*A1)*IC, 'k')
set(gca,'ylim',[0,0.7],'xticklabel',xd);
axis([0 4 0 0.7])
ylabel('Probability');
xlabel('N_B');
str1(1) = {'A'};

N=20
A2=zeros(N+1,N+1);

for i=1:N+1
    A2(i,i)=-(N-(i-1))*k1-(i-1)*k2;
end

for i=1:N
    A2(i,i+1)=i*k2;
    A2(i+1,i)=(N-(i-1))*k1;
end

IC=zeros(N+1,1);
IC(1)=1;

%produce Figure 7.41B
subplot(1,3,2)
set(gca,'fontsize',14)
bar(expm(10*A2)*IC, 'k')
ylabel('Probability');
xlabel('N_B');
axis([0 21 0 0.25])
str1(1) = {'B'};


N=200
A3=zeros(N+1,N+1);

for i=1:N+1
    A3(i,i)=-(N-(i-1))*k1-(i-1)*k2;
end

for i=1:N
    A3(i,i+1)=i*k2;
    A3(i+1,i)=(N-(i-1))*k1;
end

IC=zeros(N+1,1);
IC(1)=1;

%produce Figure 7.41C
subplot(1,3,3)
set(gca,'fontsize',14)
bar(expm(10*A3)*IC, 'k')
axis([0 201 0 0.07])
ylabel('Probability');
xlabel('N_B');
str1(1) = {'C'};


