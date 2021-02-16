close all
clear;clc

syms w x y z x1 x2 m1 m2 real

r1 = sqrt((x-x1)^2+y^2+z^2);
r2 = sqrt((x-x2)^2+y^2+z^2);
U = w^2/2*(x^2+y^2) + m1/r1 + m2/r2;

dudx = simplify(diff(U,x))
dudy = simplify(diff(U,y))
dudz = simplify(diff(U,z))

% assuming w = 1
r1 = sqrt((x-x1)^2+y^2+z^2);
r2 = sqrt((x-x2)^2+y^2+z^2);
U = 1/2*(x^2+y^2) + m1/r1 + m2/r2;

dudx = simplify(diff(U,x))
dudy = simplify(diff(U,y))
dudz = simplify(diff(U,z))