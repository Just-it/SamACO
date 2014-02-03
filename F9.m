function [ F ] = F9( x )
%F9 Summary of this function goes here
%   [-5.12,5.12]
x = x';
F = sum(x .^ 2 - 10 * cos(2 * pi * x) + 10 * ones(size(x)));
end

