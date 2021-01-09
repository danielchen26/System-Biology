%file stochastic_brusselator.m
%stochastic implementation (Gillespie SSA) of the brusselator model
%Figure 7.47


clear all

%set parameter values
k1=5000;
k2=50;
k3=0.00005;
k4=5;


%set initial condition for molecule counts and for time t, 
X=[0 1000 2000]


%set number of reactions
%for Figure 7.47, Tend=500000 (takes a while to run)
Tend=30000;


j=1;
for j=1:Tend
%calculate propensities
a1=k1;
a2=k2*X(j,2);
a3=k3*X(j,3)*X(j,2)*(X(j,2)-1)/2;
a4=k4*X(j,2);

asum=a1+a2+a3+a4;

%update counter
j=j+1;

%update time
X(j,1)=X(j-1,1)+log(1/rand(1))/asum;

%state transition
mu=rand(1);
if 0 <= mu && mu < a1/asum
X(j,2)=X(j-1,2)+1;
X(j,3)=X(j-1,3);
    else if a1/asum <= mu && mu  < (a1+a2)/asum
        X(j,2)=max(X(j-1,2)-1,0);
        X(j,3)=X(j-1,3)+1;
        else if (a1+a2)/asum <= mu && mu < (a1+a2+a3)/asum
            X(j,2) = X(j-1,2)+1;
            X(j,3)=max(X(j-1,3)-1,0);
            else
            X(j,2)=max(X(j-1,2)-1,0);
            X(j,3) =X(j-1,3);
            end
        end
end
     %j=j+1;
end



%produce Figure 7.47
figure(1)
subplot(1,2,1)
c=[0 0 0;0.4 0.4 0.4;0.8 0.8 0.8]
%setup axes
set(0,'DefaultAxesColorOrder',c)
stairs(X(:,1), X(:,3), 'k', 'linewidth', 2)
hold on
stairs(X(:,1), X(:,2), 'color', [0.5 0.5 0.5],'linewidth', 2)
axis([0 4.5 0 7000])
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
legend('X', 'Y')
str1(1) = {'A'};
text(-1.2,7000,str1, 'Fontsize', 40)
hold off

subplot(1,2,2)
stairs(X(:,2), X(:,3), 'k', 'linewidth', 2)
hold on
axis([0 7000 0 7000])
set(gca,'fontsize',12)
ylabel('Number of X Molecules')
xlabel('Number of Y Molecules')
str1(1) = {'B'};
text(-1900,7000,str1, 'Fontsize', 40)
hold off


