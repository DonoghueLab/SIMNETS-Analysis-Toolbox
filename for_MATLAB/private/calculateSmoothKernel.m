function [smoothKernel, smoothSteps] = calculateSmoothKernel(smoothKernel, smoothSteps)
% 	CALCULATESMOOTHKERNEL   Given a kernel and a number of steps, returns a vector with kernel weights
% 		[KERN] = CALCULATESMOOTHKERNEL(SMOOTHSTEPS)
%
% 	Inputs:
%		smoothSteps:		scalar integer {5}, or length of smoothKernel if it is a vector
% 		smoothKernel:		{'gauss'}, 'box', or numeric vector. For gauss, sigma will be 1/10 of the window length
% 							(so either side covers 5 standard deviations).
% 							smoothKernel will be normalized to 1.
% 	Outputs:
% 		smoothKernel:		normalized kernel, calculated from input smoothKernel
% 		smoothSteps:		length(smoothKernel)
%
% 	Created by Jonas B Zimmermann on 2015-01-30.
% 	Copyright (c) 2015 Donoghue Lab, Neuroscience Department, Brown University, Providence, RI.
% 	All rights reserved.
% 	@author  $Author: Jonas B Zimmermann$
if nargin < 2
	if nargin < 1
		smoothKernel = 'gauss';
	end
	smoothSteps = 5;
end

if ischar(smoothKernel)
	if isempty(smoothSteps) || ~isnumeric(smoothSteps)
		smoothSteps = 5;
	end
end
if ischar(smoothKernel) && strcmpi(smoothKernel, 'box')
	smoothKernel = ones(smoothSteps, 1);
elseif ischar(smoothKernel)
	if strcmpi(smoothKernel, 'half_gauss')
		nst = floor(smoothSteps / 2);
		sigma = nst / 5;
		smoothKernel = exp(-.5*(((-nst:nst)/sigma).^2));
		smoothKernel(1:nst) = 0;
		smoothKernel = smoothKernel / sum(smoothKernel);
	else
		if ~strcmpi(smoothKernel, 'gauss')
			warning('TTPlot:calculateSmoothKernel:smoothKernel','Smoothing kernel not recognized (''%s''). Using default (''gauss'') instead.', smoothKernel);
		end
		nst = floor(smoothSteps / 2);
		sigma = nst / 5;
		smoothKernel = exp(-.5*(((-nst:nst)/sigma).^2));
	end
elseif isnumeric(smoothKernel) && ~isempty(smoothSteps);
	smoothSteps = length(smoothKernel);
else
	error('TTPlot:calculateSmoothKernel:smoothKernel', 'Argument ''smooth_kernel'' has to be ''box'', ''gauss'', or a numeric vector.');
end

% Finish by normalizing kernel:
smoothKernel = smoothKernel / sum(smoothKernel);

end %  function
