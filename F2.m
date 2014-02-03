function [ F ] = F2( x )
%F3 Summary of this function goes here
%   [-10,10]
x = x';
F = sum(abs(x)) + prod(abs(x));
end

