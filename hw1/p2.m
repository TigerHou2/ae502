%% AE 502 HW1 P2, SP21
% Tiger Hou

%% initialization
close all
clear;clc
latexify

G = 1;
% mu = 0.012;  % Earth-Moon
mu = 5.974e24 / 1.989e30;  % Sun-Earth

m2 = mu;
m1 = 1 - mu;
R = 1;

alpha = m2 / (m1+m2);
beta  = m1 / (m1+m2);

x1 = -alpha*R; % construct the system with barycenter at the origin
x2 =   beta*R;

%% create L1-L2 contours
% hw = 0.2;
hw = 0.02;
x = linspace(x2-hw,x2+hw,600);
y = linspace(-hw,hw,600);
[X,Y] = meshgrid(x,y);

R1 = sqrt((X-x1).^2+Y.^2);
R2 = sqrt((X-x2).^2+Y.^2);
Z = X.^2 + Y.^2 + 2*(1-mu)./R1 + 2*mu./R2;

%% plot results for L1-L2
figure(1)
new_bone = bone;
new_bone = new_bone(28:60,:);
colormap(new_bone)
hold on
% set(gca,'ColorScale','log')
% [M,c] = contourf(X,Y,Z,[1,3.12, 3.171, 3.187, 3.24, 3.3]);
[M,c] = contourf(X,Y,Z,[1,2.999, 3.0007, 3.000891, 3.0012]);
clabel(M,c)
xlabel('x')
ylabel('y')
scatter(x2,0,'b','filled')
hold off
axis equal
setgrid
latexify(16,12,12)

%% create L3-L5 contours
x = linspace(-1,1.5,1500);
y = linspace(-1.25,1.25,1500);
[X,Y] = meshgrid(x,y);

R1 = sqrt((X-x1).^2+Y.^2);
R2 = sqrt((X-x2).^2+Y.^2);
Z = X.^2 + Y.^2 + 2*(1-mu)./R1 + 2*mu./R2;

%% plot results for L3-L5
figure(2)
hold on
% set(gca,'ColorScale','log')
colormap(new_bone)
% [M,c] = contourf(X,Y,Z,[1, 3, 3.1, 3.18, 3.5, 5]);
[M,c] = contourf(X,Y,Z,[1,2.999, 3.00, 3.2, 3.5, 4, 5]);
clabel(M,c)
xlabel('x')
ylabel('y')
scatter(x1,0,'r','filled')
scatter(x2,0,'b','filled')
hold off
axis equal
setgrid
latexify(16,12,12)
