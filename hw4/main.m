close all
clear;clc
addpath('vallado')

%% parse orbit data

dat = importdata('project4_data.txt');

TE = dat(:,1); % seconds
RA = dat(:,2); % radians
DE = dat(:,3); % radians

%% constants and settings

mu = 398600.44; % km^3/s^2
R = 6378.1; % km
lat = pi/6;
lon = 0;
omega = rad2deg(7.2936e-5);
site = R * [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]';
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% problem 1

% we take the 1st measurement in each track
% indices = 1 + 5*(0:2);
indices = [1, 8, 15];
te_vec = TE(indices);
ra_vec = RA(indices);
de_vec = DE(indices);
site1 = rotz(te_vec(1)*omega) * site;
site2 = rotz(te_vec(2)*omega) * site;
site3 = rotz(te_vec(3)*omega) * site;

% use Gaussian IOD from Vallado
[r,v] = anglesg(de_vec(1),de_vec(2),de_vec(3),...
                ra_vec(1),ra_vec(2),ra_vec(3),...
                te_vec(1),te_vec(2),te_vec(3),...
                site1, site2, site3, R, mu, 1);
            
clc
            
[a,e,i,o,w,f] = Get_Orb_Params(r,v,mu);
disp( 'Problem 1:' )
disp(['    a = ' num2str(a) ' km'])
disp(['    e = ' num2str(norm(e))])
disp(['    i = ' num2str(rad2deg(i)) ' deg'])
disp(['    o = ' num2str(rad2deg(o)) ' deg'])
disp(['    w = ' num2str(rad2deg(w)) ' deg'])
disp(['    f = ' num2str(rad2deg(f)) ' deg'])
disp(' ')


%% problem 2

x0_hat = [r,v]';
x_star = x0_hat;
x0_bar = zeros(6,1);
STM = eye(6);
P0_bar = diag(1e10*ones(6,1));
Ri = diag([(1/3600*pi/180)^2,(1/3600*pi/180)^2]);
L = zeros(6,6);
N = zeros(6,1);

%% define EOMs for ode45

xdot = @(X) [X(4:6);-mu*X(1)/norm(X(1:3))^3; ...
                    -mu*X(2)/norm(X(1:3))^3; ...
                    -mu*X(3)/norm(X(1:3))^3; ];
Q = @(X) [      mu/norm(X(1:3))^5*(3*X(1)^2-norm(X(1:3))^2), ...
                3*mu*X(1)*X(2)/norm(X(1:3))^5, ...
                3*mu*X(1)*X(3)/norm(X(1:3))^5; ...
                ...
                3*mu*X(2)*X(1)/norm(X(1:3))^5, ...
                mu/norm(X(1:3))^5*(3*X(2)^2-norm(X(1:3))^2), ...
                3*mu*X(2)*X(3)/norm(X(1:3))^5; ...
                ...
                3*mu*X(3)*X(1)/norm(X(1:3))^5, ...
                3*mu*X(3)*X(2)/norm(X(1:3))^5, ...
                mu/norm(X(1:3))^5*(3*X(3)^2-norm(X(1:3))^2); ...
            ];
A = @(X) [   zeros(3), eye(3); ...
             Q(X), zeros(3)];
STM_dot = @(X,STM)  reshape( A(X) * reshape(STM,6,6), 36, 1);
combined_dot = @(combined) [xdot(combined(1:6));...
                            STM_dot(combined(1:6),combined(7:42))];
                        
%% back propagate to epoch

[~,out] = ode45(@(t,y) xdot(y), [te_vec(2),0], x_star, options);
combined_orig = out(end,:)';
combined_orig = [combined_orig(1:6); reshape(STM,36,1)];
combined = combined_orig;

%% perform batch processing

len = length(TE);

while (1)
    
    tp = 0;
    
    for i = 1:len
        t = TE(i);
        Y = [RA(i),DE(i)]';
        
        [timestep,combined] = ...
            ode45(@(t,y) combined_dot(y), [tp,t], combined, options);
        combined = combined(end,:)';
        STM = reshape(combined(7:42),6,6);
        x = combined(1);
        y = combined(2);
        z = combined(3);
        
        R = rotz(t*omega) * site;
        Px = x - R(1);
        Py = y - R(2);
        Pz = z - R(3);
        P = norm([Px,Py,Pz]);
            
        Gy = [atan2(Py,Px),asin(Pz/P)]';
        Ht = [  -sin(Gy(1))*cos(Gy(1)) / Px, ...
                cos(Gy(1))^2 / Px, ...
                0, 0, 0, 0; ...
                -tan(Gy(2))*Px/P^2, ...
                -tan(Gy(2))*Py/P^2, ...
                cos(Gy(2))/P, ...
                0, 0, 0];
        
        yi = Y - Gy;
        Hi = Ht * STM;
        
        L = L + Hi' * inv(Ri) * Hi;
        N = N + Hi' * inv(Ri) * yi;
    
        tp = t;
        
    end
    
    x0_hat = L \ N;
    metric = norm(x0_hat);
    P0 = inv(L);
    
    if metric < 1e-6
        break
    end
    
    L = inv(P0_bar);
    N = inv(P0_bar) * x0_bar;
    
    combined_orig = combined_orig + [x0_hat;zeros(36,1)];
    combined = combined_orig;
    x0_bar = x0_bar - x0_hat;
    
end

r = combined(1:3);
v = combined(4:6);

%% propagate and show results

[~,out] = ode45(@(t,y) xdot(y), [TE(1),0], [r;v], options);
u = out(end,:)';

r0 = u(1:3);
v0 = u(4:6);
[a,e,i,o,w,f] = Get_Orb_Params(r0,v0,mu);
disp( 'Problem 2:' )
disp(['    a = ' num2str(a) ' km'])
disp(['    e = ' num2str(norm(e))])
disp(['    i = ' num2str(rad2deg(i)) ' deg'])
disp(['    o = ' num2str(rad2deg(o)) ' deg'])
disp(['    w = ' num2str(rad2deg(w)) ' deg'])
disp(['    f = ' num2str(rad2deg(f)) ' deg'])

[~,orbit] = ode45(@(t,y) xdot(y), 0:20:50000, [r;v], options);
[~,data] = ode45(@(t,y) xdot(y), 0:20:TE(15), [r;v], options);

figure(1)
hold on
plot3(orbit(:,1), orbit(:,2),orbit(:,3),'r','LineWidth',1.2)
plot3(data(:,1), data(:,2),data(:,3),'b','LineWidth',1.2)
scatter3(0,0,0,10,'k','filled')
hold off
view(3)
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
axis equal
setgrid
latexify(16,14,14)