function [r,theta] = cart2polar(x)
% convert a 2d vector from cartesian coordinate to polar
r=norm(x(1:2));
theta=imag(log(x(1) + j*x(2)));
