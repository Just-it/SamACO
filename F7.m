function [ F ] = F7( x )
%F7 Summary of this function goes here
%   [-1.28,1.28]
x = x';
F = zeros(1, size(x, 2));
for i = 1 : size(x, 1)
    F = F + i * x(i, :) .^ 4;
end
end

