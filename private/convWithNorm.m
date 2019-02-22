function [v] = convWithNorm(data, kernel)
% 	CONVWITHNORM   Convolves data with kernel using MATLAB's inbuilt conv(data, kernel, 'same')
% 	function and normalizes at both ends of the data vector.
% 		[V] = CONVWITHNORM(DATA, KERNEL)
%
%
% 	Created by Jonas B Zimmermann on 2014-10-06.
% 	Copyright (c) 2014 Donoghue Lab, Neuroscience Department, Brown University.
% 	All rights reserved.
% 	@author  $Author: Jonas B Zimmermann$

nst = floor(length(kernel) / 2);
% normalize kernel
kernel = kernel(:) ./ sum(kernel);

dflip = false;
if size(data, 1) == 1
	data = data';
	dflip = true;
end

% get weights for both ends of kernel, in case it's not symmetric
weight1 = sum(kernel) ./ cumsum(kernel);
weight2 = sum(kernel) ./ cumsum(flip(kernel));

v = zeros(size(data));

for i_c = 1:size(data, 2)
% convolve data with kernel
v(:, i_c) = conv(data(:, i_c), kernel, 'same');
% correct convolved data at either ends
v(1:(nst+1), i_c) = v(1:(nst+1), i_c) .* weight1((end-nst):end);
v(end-nst:end, i_c) = v(end-nst:end, i_c) .* weight2(end:-1:end-nst);
end

if dflip
	v = v';
end

end %  function
