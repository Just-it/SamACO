function [ F ] = F5( x )
%F5 Summary of this function goes here
%   [-30,30]
x = x';
F = zeros(1, size(x, 2));
for i = 1 : size(x, 1) - 1
    F = F + 100 * (x(i + 1, :) - x(i, :) .^ 2) .^ 2 + (1 - x(i, :)) .^ 2;
end
end

