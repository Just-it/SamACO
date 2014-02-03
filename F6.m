function [ F ] = F6( x )
%F6 Summary of this function goes here
%   [-100,100]
x = x';
F = sum(floor(x + 0.5 * ones(size(x))) .^ 2);

end

