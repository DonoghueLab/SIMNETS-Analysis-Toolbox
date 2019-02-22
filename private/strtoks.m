function y=strtoks(x,s,needle)
% y = STRTOKS(x) returns a cell array of strings consisting of the
% space delimited parts of X.
% y = STRTOKS(x,s) uses S instead of space.
% y = STRTOKS(...,needle) only returns those tokens that contain NEEDLE
% as a substring.

% Author: Daniel A. Wagenaar
% Released under GPL Version 2, June 1991

if nargin<2
  s=[];
end

if iscell(x)
  y=x;
  return
end

n=1;
y=cell(0,1);
while length(x)>0
  if ~isempty(s)
    [ z, x ] = strtok(x,sprintf(s));
  else
    [ z, x ] = strtok(x);
  end
  if length(z)
    y{n}=z;
    n=n+1;
  end
end

if nargin>=3
  z=cell(0,1);
  m=1;
  for n=1:length(y)
    if strfind(y{n},needle)
      z{m}=y{n};
      m=m+1;
    end
  end
  y=z;
end
