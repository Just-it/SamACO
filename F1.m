function [ F ] = F1( x )
%F2 Summary of this function goes here
%   [-100, 100]
F = sum(x .^ 2, 2)';
end