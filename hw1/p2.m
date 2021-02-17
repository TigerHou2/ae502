%% AE 502 HW1 P2, SP21
% Tiger Hou

%% initialization
close all
clear;clc

G = 1;
mu = 0.012;

m2 = mu;
m1 = 1 - mu;
R = 1;

alpha = m2 / (m1+m2);
beta  = m1 / (m1+m2);

x1 = -alpha*R; % construct the system with barycenter at the origin
x2 =   beta*R;

%% create contour
x = linspace(-1,1.5,2000);
y = linspace(-1.25,1.25,2000);
[X,Y] = meshgrid(x,y);

R1 = sqrt((X-x1).^2+Y.^2);
R2 = sqrt((X-x2).^2+Y.^2);
Z = X.^2 + Y.^2 + 2*(1-mu)./R1 + 2*mu./R2;

%% plot results
figure(1)
hold on
scatter(x1,0,'r','filled')
scatter(x2,0,'b','filled')
set(gca,'ColorScale','log')
[M,c] = contourf(X,Y,Z,[3,3.05,3.1,1.15,3.2,3.5,4,4.5,5,5.5]);
clabel(M,c)
% contourf(X,Y,Z,[3.05,3.1])
hold off
axis equal
