%% AE 502 HW1 P3, SP21
% Tiger Hou

%% initialization
close all
clear;clc
addpath('../../tools')
latexify

mu = 0.012150585609262;
x0 = [  1.118824382902157, 0.0, ...
        0.014654873101278, 0.0, ...
        0.180568501159703, 0.0]';
P = 1.706067405636607;

G = 1;
m2 = mu;
m1 = 1 - mu;
R = 1;

alpha = m2 / (m1+m2);
beta  = m1 / (m1+m2);

x1 = -alpha*R; % construct the system with barycenter at the origin
x2 =   beta*R;

w = sqrt(G*(m1+m2)/R^3);

%% state transition matrix

phi_t0 = eye(6);
Y0 = [reshape(phi_t0,36,1);x0];

dudx = @(x,y,z) x*w^2 ...
     - (m1*(2*x - 2*x1))/(2*((x - x1)^2 + y^2 + z^2)^(3/2)) ...
     - (m2*(2*x - 2*x2))/(2*((x - x2)^2 + y^2 + z^2)^(3/2));
 
dudy = @(x,y,z) w^2*y ...
     - (m1*y)/((x - x1)^2 + y^2 + z^2)^(3/2) ...
     - (m2*y)/((x - x2)^2 + y^2 + z^2)^(3/2);
 
dudz = @(x,y,z) ...
     - (m1*z)/((x - x1)^2 + y^2 + z^2)^(3/2) ...
     - (m2*z)/((x - x2)^2 + y^2 + z^2)^(3/2);

dudx_dx = @(x,y,z) w^2 ...
        - m2/((x - x2)^2 + y^2 + z^2)^(3/2) ...
        - m1/((x - x1)^2 + y^2 + z^2)^(3/2) ...
        + (3*m1*(2*x - 2*x1)^2)/(4*((x - x1)^2 + y^2 + z^2)^(5/2)) ...
        + (3*m2*(2*x - 2*x2)^2)/(4*((x - x2)^2 + y^2 + z^2)^(5/2));

dudy_dy = @(x,y,z) w^2 ...
        - m2/((x - x2)^2 + y^2 + z^2)^(3/2) ...
        - m1/((x - x1)^2 + y^2 + z^2)^(3/2) ...
        + (3*m1*y^2)/((x - x1)^2 + y^2 + z^2)^(5/2) ...
        + (3*m2*y^2)/((x - x2)^2 + y^2 + z^2)^(5/2);

dudz_dz = @(x,y,z) ...
          (3*m1*z^2)/((x - x1)^2 + y^2 + z^2)^(5/2) ...
        - m2/((x - x2)^2 + y^2 + z^2)^(3/2) ...
        - m1/((x - x1)^2 + y^2 + z^2)^(3/2) ...
        + (3*m2*z^2)/((x - x2)^2 + y^2 + z^2)^(5/2);

dudx_dy = @(x,y,z) ...
          (3*m1*y*(2*x - 2*x1))/(2*((x - x1)^2 + y^2 + z^2)^(5/2)) ...
        + (3*m2*y*(2*x - 2*x2))/(2*((x - x2)^2 + y^2 + z^2)^(5/2));

dudy_dz = @(x,y,z) ...
          (3*m1*y*z)/((x - x1)^2 + y^2 + z^2)^(5/2) ...
        + (3*m2*y*z)/((x - x2)^2 + y^2 + z^2)^(5/2);

dudz_dx = @(x,y,z) ...
          (3*m1*z*(2*x - 2*x1))/(2*((x - x1)^2 + y^2 + z^2)^(5/2)) ...
        + (3*m2*z*(2*x - 2*x2))/(2*((x - x2)^2 + y^2 + z^2)^(5/2));

dudy_dx = @(x,y,z) dudx_dy(x,y,z);
dudz_dy = @(x,y,z) dudy_dz(x,y,z);
dudx_dz = @(x,y,z) dudz_dx(x,y,z);
    
x_dot = @(x,y,z,vx,vy,vz) ...
        [   vx; vy; vz; ...
            2*vy + dudx(x,y,z); ...
           -2*vx + dudy(x,y,z); ...
            dudz(x,y,z)];

F = @(x,y,z) ...
    [0 0 0 1 0 0; ...
     0 0 0 0 1 0; ...
     0 0 0 0 0 1; ...
     dudx_dx(x,y,z) dudx_dy(x,y,z) dudx_dz(x,y,z)  0 2 0;...
     dudy_dx(x,y,z) dudy_dy(x,y,z) dudy_dz(x,y,z) -2 0 0;...
     dudz_dx(x,y,z) dudz_dy(x,y,z) dudz_dz(x,y,z)  0 0 0];
 
phi_dot = @(yy) F(yy(37),yy(38),yy(39)) * reshape(yy(1:36),6,6);

Y_dot = @(yy) ...
        [   reshape(phi_dot(yy),36,1); ...
            x_dot(yy(37),yy(38),yy(39),yy(40),yy(41),yy(42))];

tspan = 2*[0,P];

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,Y] = ode45(@(t,yy) Y_dot(yy), tspan, Y0, options);

phi_tp = reshape(Y(end,1:36)',6,6);

% sanity check - does the orbit complete 1 revolution at t = P?
max(Y(end,37:42)'-Y0(37:42))

%% compute eigenvalues of the STM at time P
%   to find stable and unstable manifolds

[V,D] = eig(phi_tp);

% we only care about the first two eigenvectors
wu = real(V(:,1)); % unstable manifold
ws = real(V(:,2)); %   stable manifold


%% plot manifolds
n = 10;
t_arr = linspace(0,2*P,n+1);
t_arr(1) = [];
du = 0.001; % disturbance to unstable manifold
ds = 0.001; % disturbance to stable manifold

% plot unstable manifold
for i = 1:n
    
    % propagate to a certain timestep
    tspan = [0,t_arr(i)];
    Y0 = [reshape(phi_t0,36,1);x0];
    [~,Y] = ode45(@(t,yy) Y_dot(yy), tspan, Y0, options);
    phi_tp = reshape(Y(end,1:36)',6,6);
    
    % compute direction of manifolds
    wdir_u = phi_tp * wu;
    wdir_s = phi_tp * ws;

    % normalize directions
    wdir_u = wdir_u / norm(wdir_u);
    wdir_s = wdir_s / norm(wdir_s);
    
    % get the starting position
    xx = Y(end,37:42)';
    
    % plot unstable manifold
    xw = xx + du*wdir_u;
    tspan = [0,2*P];
    Y0 = [reshape(phi_t0,36,1);xw];
    [~,Y] = ode45(@(t,yy) Y_dot(yy), tspan, Y0, options);
    figure(1)
    hold on
    plot3(Y(:,37),Y(:,38),Y(:,39),'r','LineWidth',1.2);
    hold off
    
    % plot stable manifold
    xw = xx + ds*wdir_s;
    tspan = [0,2*P];
    tspan = fliplr(tspan);
    Y0 = [reshape(phi_t0,36,1);xw];
    [~,Y] = ode45(@(t,yy) Y_dot(yy), tspan, Y0, options);
    figure(2)
    hold on
    plot3(Y(:,37),Y(:,38),Y(:,39),'b','LineWidth',1.2);
    hold off
    
end

figure(1)
view(3)
hold on
moon = scatter3(x2,0,0,'b','filled');
hold off
grid on
grid minor
title(['Unstable Manifolds, d = ' num2str(du)])
xlabel('x')
ylabel('y')
zlabel('z')
legend(moon,'The Moon','Location','Best')
latexify(16,12,16)

figure(2)
view(3)
hold on
moon = scatter3(x2,0,0,'r','filled');
hold off
grid on
grid minor
title(['Stable Manifolds, d = ' num2str(ds)])
xlabel('x')
ylabel('y')
zlabel('z')
legend(moon,'The Moon','Location','Best')
latexify(16,12,16)
