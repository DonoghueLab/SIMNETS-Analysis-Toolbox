function ov = parseBinArg(iv, default)
% 	PARSEBINARG   Parses input for 'yes'/'no'/true/false and returns true/false
% 		[OV] = PARSEBINARG(IV, DEFAULT)
%
% 	Parses argument iv for 'yes'/'no'/true/false/1/0 and returns true or false.
% 	Optional argument default determines return value if iv is empty (flase by default)
%
% 	@author Jonas B Zimmermann
% 	@version $Id:$
%
% 	Created by Jonas B Zimmermann on 2014-08-07.
% 	Copyright (c) 2014 Jonas B Zimmermann & Brown University.
% 	All rights reserved.

if nargin==1
	default = false;
end
if isempty(iv)
	ov = default;
	return;
end
if (ischar(iv) && strcmpi(iv, 'yes')) || (~ischar(iv) && iv)
	ov = true;
else
	ov = false;
end

end %  function
