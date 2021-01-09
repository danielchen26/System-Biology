% file stochastic_vilar_oscillator.m
% stochastic simulation of genetic oscillator
% from Vilar et al. (2002) PNAS 99 5988-5992.
% problem 7.8.27

clear all

%kinetic parameter values

eps=0.1;
ga=250;
ba=5;
Ka=0.5; 
alpha0=0.1;
da=1; 
kc=200; 
gr=50; 
br=10;
Kr=1; 
dr=eps*da;

%set initial condition for molecule counts and for time t, 
%state: [time, A, R, C]
X=[0 0 0 0];

%set number of reactions  
Tend=20000;

j=1;
for j=1:Tend
%while X(j,1) < 0.2

A=X(j,2); R=X(j,3); C=X(j,4);

%calculate propensities
a1=(ga/ba)*(alpha0+(A/Ka))/(1+(A/Ka));
a2=da*A;
a3=kc*A*R;
a4=(gr/br)*((A/Kr))/(1+(A/Kr));
a5=dr*R;
a6=da*C;

asum=a1+a2+a3+a4+a5+a6;


%update counter
j=j+1;

%update time
X(j,1)=X(j-1,1)+log(1/rand(1))/asum;

%state transition
mu=rand(1);
if 0 <= mu && mu < a1/asum
X(j,2)=X(j-1,2)+ba;
X(j,3)=X(j-1,3);
X(j,4)=X(j-1,4);
else if a1/asum <= mu && mu  < (a1+a2)/asum
X(j,2)=X(j-1,2)-1;
X(j,3)=X(j-1,3);
X(j,4)=X(j-1,4);
else if (a1+a2)/asum <= mu && mu < (a1+a2+a3)/asum
X(j,2)=X(j-1,2)-1;
X(j,3)=X(j-1,3)-1;
X(j,4)=X(j-1,4)+1;
    else if (a1+a2+a3)/asum <= mu && mu < (a1+a2+a3+a4)/asum
X(j,2)=X(j-1,2);
X(j,3)=X(j-1,3)+br;
X(j,4)=X(j-1,4);
else if (a1+a2+a3+a4)/asum <= mu && mu < (a1+a2+a3+a5)/asum
X(j,2)=X(j-1,2);
X(j,3)=X(j-1,3)-1;
X(j,4)=X(j-1,4);
    else
X(j,2)=X(j-1,2);
X(j,3)=X(j-1,3)+1;
X(j,4)=X(j-1,4)-1;
    end
    end
    end
    end
end

end

figure(1)
stairs(X(:,1), X(:,2), 'g', 'linewidth', 2)
hold on
stairs(X(:,1), X(:,3), 'b', 'linewidth', 2)
hold on
%axis([0 4.5 0 7000])
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
legend('A', 'R')
hold off



