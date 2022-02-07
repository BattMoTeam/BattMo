function [ sw ] = lineup(y1, y2, x1, x2, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dy = y2-y1;
dx = (x2-x1);

m = dy/dx;

x = x-x1;

sw = (x >= x1 & x<= x2) .* m.*x + y1 + ...
    (x>x2) .* y2 + ...
    (x<x1) .* y1;

end
