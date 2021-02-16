%% AE 502 HW1 P1, SP21
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

w = sqrt(G*(m1+m2)/R^3);
T = 40;
dt = 0.001;
steps = ceil(T/dt);

% Lagrange points
% https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec18.pdf
L1 = R*(1-(alpha/3)^(1/3));
L2 = R*(1+(alpha/3)^(1/3));
L3 = -R*(1+5*alpha/12);

% u0 = [ 0.5,  0.5,  0.0,  1.0,  0.5,  0.0]';
u0 = [ 0.5,  0.5*sqrt(3),  0.0,  0.0,  0.0,  0.0; ... % L4
       0.5, -0.5*sqrt(3),  0.0,  0.0,  0.0,  0.0; ... % L5
       (1+1e-6)*L3, 0.0,  0.0,  0.0,  -1.00002e-1,  0.0; ... % weird L3
       L3, 1e-3, 0, 0, +6.1e-5, 0 ]'; % L3 with perturbation + vel adj

N = size(u0,2);
uu = zeros(6,N,steps);

%% propagation

for jj = 1 : steps
    
    uu(:,:,jj) = u0;
    % ============= RK4 =============
    k1 = dt * ff( u0,N,m1,m2,x1,x2,w);
    k2 = dt * ff((u0 + k1/2),N,m1,m2,x1,x2,w);
    k3 = dt * ff((u0 + k2/2),N,m1,m2,x1,x2,w);
    k4 = dt * ff((u0 + k3)  ,N,m1,m2,x1,x2,w);
    u0 = u0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end


%% plotting - 3D

dat = permute(uu,[1,3,2]);

figure(1)
hold on
scatter3(x1,0,0,'r','filled')
scatter3(x2,0,0,'b','filled')
for i = 1:N
    plot3(dat(1,:,i),dat(2,:,i),dat(3,:,i))
end
hold off
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid minor


%% function definitions

function rv_dot = ff(rv,N,m1,m2,x1,x2,w)
% takes the 6xN position & velocity vector and computes the acceleration
%   where N is the number of particles to track

rv_dot = zeros(6,N);

for i = 1:N
    x = rv(1,i);
    y = rv(2,i);
    z = rv(3,i);

    dudx = x*w^2 ...
         - (m1*(2*x - 2*x1))/(2*((x - x1)^2 + y^2 + z^2)^(3/2)) ...
         - (m2*(2*x - 2*x2))/(2*((x - x2)^2 + y^2 + z^2)^(3/2));
    dudy = w^2*y ...
         - (m1*y)/((x - x1)^2 + y^2 + z^2)^(3/2) ...
         - (m2*y)/((x - x2)^2 + y^2 + z^2)^(3/2);
    dudz = - (m1*z)/((x - x1)^2 + y^2 + z^2)^(3/2) ...
           - (m2*z)/((x - x2)^2 + y^2 + z^2)^(3/2);

    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    rv_dot(4:6,i) = [dudx + 2*w*rv(5,i);...
                     dudy - 2*w*rv(4,i);...
                     dudz];
end

end