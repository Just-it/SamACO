function [ F ] = F4( x )
%F4 Summary of this function goes here
%   [-100,100]
x = x';
F=max(abs(x));
end