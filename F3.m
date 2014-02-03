function [ F ] = F3( x )
%F4 Summary of this function goes here
%   [-100,100]
x = x';
F = zeros(1, size(x, 2));
f1 = zeros(1, size(x, 2));
for i = 1 : size(x, 1)
    f1 = f1 + x(i, :);
    F = F + f1 .^ 2;
end
end