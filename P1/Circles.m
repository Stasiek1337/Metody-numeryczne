function [k] = Circles(x, y, r)
% opis

i = rand; j = rand; k = rand;
t = 0:pi/150:2*pi;
x = r * cos(t) + x;
y = r * sin(t) + y;
k = fill(x, y,[i, j, k]);
