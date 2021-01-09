%file potentialwells.m
%plot potential surfaces in Figure 4.11.  The plots can be rotated to
%match the figure.

%single-well potential

%produce grid
sx=linspace(-1,1,80)
sy=linspace(-2,2,80)
[s1,s2]=meshgrid(sx,sy)
%plot surface
figure(1)
hold on
surface(s1,s2, s1.^2 + (0.5*s2).^2, 'EdgeColor','none')
axis([-2 2 -2 2 0 2])
colormap hsv
colormap(gray)
shading flat

%plot (approximate) trajectories
x1=zeros(1,10);
y1=linspace(0,2,10)
plot3(x1,y1,(x1).^2 + (0.5*y1).^2, 'k', 'linewidth', 2)

x3=linspace(0,1,10);
y3=zeros(1,10)
plot3(x3,y3,(x3).^2 + (0.5*y3).^2, 'k', 'linewidth', 2)

x2=linspace(0.75, 0, 50)
y2=2.3*sqrt((x2))
plot3(x2,y2,(x2).^2 + (0.5*y2).^2, 'k', 'linewidth', 2)

set(gca,'Visible','off')

%double-well potential
sx=linspace(-2.75,2.75,80)
sy=linspace(-0.75,0.75,80)
[s1,s2]=meshgrid(sx,sy)
figure(2)
hold on
surface(s1,s2, (0.2*s1.^2-1).^2 + s2.^2, 'EdgeColor','none')
axis([-2.5 2.5 -2 2 0 2])
colormap hsv
colormap(gray)
shading flat

x1=zeros(10);
y1=linspace(0,0.75,10)
plot3(x1,y1,1+y1.^2, 'k', 'linewidth', 2)

x2=linspace(-0.14, -2.2, 50)
y2=0.01*(x2+2.2).^6
plot3(x2,y2,(0.2*x2.^2-1).^2 + y2.^2, 'k', 'linewidth', 2)

x23=linspace(-1.14, -2.2, 50)
y23=0.62*(x23+2.2).^3
plot3(x23,y23,(0.2*x23.^2-1).^2 + y23.^2, 'k', 'linewidth', 2)

x3=linspace(0.14, 2.2, 50)
y3=0.01*(x3-2.2).^6
plot3(x3,y3,(0.2*x3.^2-1).^2 + y3.^2, 'k', 'linewidth', 2)

x33=linspace(1.14, 2.2, 50)
y33=-0.62*(x33-2.2).^3
plot3(x33,y33,(0.2*x33.^2-1).^2 + y33.^2, 'k', 'linewidth', 2)

set(gca,'Visible','off')
