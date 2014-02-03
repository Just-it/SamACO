function [ F ] = F8( x )
%F8 Summary of this function goes here
%   [-500,500]
x = x';
F = zeros(1, size(x, 2));
for i = 1 : size(x, 1)
    F = F + x(i, :) .* sin(sqrt(abs(x(i, :)))); 
end
F = -F;
end

