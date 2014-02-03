function [ F ] = F10( x )
%F10 Summary of this function goes here
%   [-32,32]
x = x';
F = zeros(1, size(x, 2));
F1 = zeros(size(x(1, :)));
F2 = zeros(size(x(1, :)));
D = size(x, 1);
for i = 1 : D
    F1 = F1 + cos(2 * pi * x(i, :)) / D;
    F2 = F2 + x(i, :) .^ 2 / D;
end
F = -20 * exp(-0.2 * sqrt(F2)) - exp(F1) + 20 + exp(1);
end

