%file surfacetangent.m
%plot the surface and tangent plane in Figure 4.13.  Plot can be rotated to
%match figure.

%set mesh
s=linspace(0.5,4,40)
[s1,s2]=meshgrid(s,s)

%generate plot
figure(1)
surf(s1,s2, s1.*s2./(1+s1+s2+ s1.*s2), s1.*s2./(1+s1+s2+ s1.*s2), 'EdgeColor','none')
hold on
ss=linspace(1,2.75,30)
[ss1,ss2]=meshgrid(ss,ss)
colormap hsv
colormap(gray)
shading flat

set(gca,'fontsize',14)
surf(ss1,ss2, 0.001+4/9 + (6/81).*(ss1-2) + (6/81).*(ss2-2), -(4/9 + (6/81).*(ss1-2) + (6/81).*(ss2-2)), 'EdgeColor','none')
xlabel('s_1')
ylabel('s_2')
zlabel('z')