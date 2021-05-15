close all
clear;clc
addpath('vallado')

%% parse orbit data

dat = importdata('project4_data.txt');

TE = dat(:,1); % seconds
RA = dat(:,2); % radians
DE = dat(:,3); % radians

%% constants

mu = 398600.44; % km^3/s^2
R = 6378.1; % km
lat = pi/6;
lon = 0;

site = R * [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]';


%% problem 1

% we take the 1st measurement in each track
indices = 1 + 5*(0:2);
te_vec = TE(indices) - TE(1);
ra_vec = RA(indices);
de_vec = DE(indices);
site1 = rotz(te_vec(1)/3600/24*360) * site;
site2 = rotz(te_vec(2)/3600/24*360) * site;
site3 = rotz(te_vec(3)/3600/24*360) * site;

% use Gaussian IOD from Vallado
[r,v] = anglesg(de_vec(1),de_vec(2),de_vec(3),...
                ra_vec(1),ra_vec(2),ra_vec(3),...
                te_vec(1),te_vec(2),te_vec(3),...
                site1, site2, site3, R, mu, 1);
            
clc
            
[a,e,i,o,w,f] = Get_Orb_Params(r,v,mu);
disp( 'Problem 1:' )
disp(['a = ' num2str(a) ' km'])
disp(['e = ' num2str(norm(e))])
disp(['i = ' num2str(rad2deg(i)) ' deg'])
disp(['o = ' num2str(rad2deg(o)) ' deg'])
disp(['w = ' num2str(rad2deg(w)) ' deg'])


%% problem 2

TE = TE-TE(1);
len = length(TE);

x0_hat = [r,v]';
STM = eye(6);
P0_bar = diag(1e10*ones(6,1));
R = diag((1/3600*pi/180)^2);
L = inv(P0_bar);
N = inv(P0_bar)*x0_hat;

xp = x0_hat;
tp = 0;



while (1)
    
    for i = 1:len
        t = TE(i);
        Y = [RA(i),DE(i)]';
        
        r = norm(xp(1:3));
        x = xp(1);
        y = xp(2);
        z = xp(3);
        
        xdot = [xp(1:3) ; -mu*x/r^3; ...
                          -mu*y/r^3; ...
                          -mu*z/r^3 ];
        Q = [   mu/r^5*(3*x^2-r^2), ...
                3*mu*x*y/r^5, ...
                3*mu*x*z/r^5; ...
                ...
                3*mu*y*x/r^5, ...
                mu/r^5*(3*y^2-r^2), ...
                3*mu*y*z/r^5; ...
                ...
                3*mu*z*x/r^5, ...
                3*mu*z*y/r^5, ...
                mu/r^5*(3*z^2-r^2); ...
                ];
        A = [   zeros(3), eye(3); ...
                Q, zeros(3)];
        if i ~= 1
            STM = ode45(@(t,y) A*STM, [tp,t], STM);
        else
            STM = eye(6);
        end
        
        R = rotz(t/3600/24*360) * site;
        Ht = [  -sin(Y(1))*cos(Y(2)) / (x-R(1)), ...
                cos(Y(1))^2 / (x-R(1)), ...
                0, 0, 0, 0; ...
                -tan(Y(2))*(x-R(1))/norm(R)^2, ...
                -tan(Y(2))*(y-R(2))/norm(R)^2, ...
                cos(Y(2))/norm(R)^2, ...
                0, 0, 0];
        
    end
    
    break
end



