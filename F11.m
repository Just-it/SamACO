function [ F ] = F11( x )
%F11 Summary of this function goes here
%   [-600,600]
x = x';
temp = repmat(sqrt(1 : size(x, 1))', 1, size(x, 2));
F = (sum(x .^ 2)) ./ 4000 - prod(cos(x ./ temp)) + 1 * ones(1, size(x, 2));
F = F';
end

