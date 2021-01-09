%file local_sensitivity.m
%generate figure 4.22

function local_sensitivity

%curve:
V=linspace(2.2,10,300);
s=3./(V-2)

%tangent line:
V2=linspace(2.7,5.3)
s2=1.5-(V2-4)*(0.75)


figure(1)
hold on
set(gca,'fontsize',14)
plot(V,s,'k', 'LineWidth',3);
plot(V2,s2,'k--', 'LineWidth',3);
axis([2 8 0 4])
xlabel('V_{max} (mM/min)')
ylabel('s (mM)')

end