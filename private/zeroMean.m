function [zm] = zeroMean(x)

[row, ~] = size(x);
zm = x-ones(row, 1) * mean(x);
   