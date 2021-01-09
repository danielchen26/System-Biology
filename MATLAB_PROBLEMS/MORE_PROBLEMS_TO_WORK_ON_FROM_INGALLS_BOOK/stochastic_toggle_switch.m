% file stochastic_toggle_switch.m
% stochastic simulation of genetic toggle switch
% problem 7.8.26

clear all

%set kinetic parameter values
d=1;
a=5;
Km=1;
beta=4;


%set initial condition for molecule counts and for time t, 

X=[0 0 a]

%set number of reaction.  
Tend=20000;

j=1;
for j=1:Tend

%calculate propensities
%expression
a1=a/(Km+X(j,3));
a2=a/(Km+X(j,2));
%decay:
a3=d*X(j,2);
a4=d*X(j,3);

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
X(j,3)=X(j-1,3)+1;
X(j,2)=X(j-1,2);
else if (a1+a2)/asum <= mu && mu < (a1+a2+a3)/asum
X(j,2)=X(j-1,2)-1;
X(j,3)=X(j-1,3);
else
X(j,3)=X(j-1,3)-1; 
X(j,2)=X(j-1,2);
    end
    end
    
end

end

figure(1)
subplot(1,2,1)
stairs(X(:,1), X(:,2), 'g', 'linewidth', 2)
hold on
stairs(X(:,1), X(:,3), 'b', 'linewidth', 2)
hold on
%axis([0 4.5 0 7000])
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
legend('X', 'Y')
hold off

subplot(1,2,2)
stairs(X(:,2), X(:,3), 'k', 'linewidth', 2)
hold on
%axis([0 7000 0 7000])
set(gca,'fontsize',12)
ylabel('Number of X Molecules')
xlabel('Number of Y Molecules')
hold off


