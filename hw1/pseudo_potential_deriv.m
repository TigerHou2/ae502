close all
clear;clc

syms w x y z x1 x2 m1 m2 real

r1 = sqrt((x-x1)^2+y^2+z^2);
r2 = sqrt((x-x2)^2+y^2+z^2);
U = w^2/2*(x^2+y^2) + m1/r1 + m2/r2;

dudx = simplify(diff(U,x))
dudy = simplify(diff(U,y))
dudz = simplify(diff(U,z))

dudx_dx = simplify(diff(dudx,x))
dudx_dy = simplify(diff(dudx,y))
dudx_dz = simplify(diff(dudx,z))

dudy_dx = simplify(diff(dudy,x))
dudy_dy = simplify(diff(dudy,y))
dudy_dz = simplify(diff(dudy,z))

dudz_dx = simplify(diff(dudz,x))
dudz_dy = simplify(diff(dudz,y))
dudz_dz = simplify(diff(dudz,z))