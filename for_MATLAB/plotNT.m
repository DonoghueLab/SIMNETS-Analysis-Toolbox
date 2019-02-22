function [varargout] = plotNT(NT,id,varargin)
% 	PLOTNT plots trajectory (or snapshots) of neural data
% 	[ALLH] = PLOTNT(NT,ID,VARARGIN)
% 	ALLH:		handles for plot objects (optional)
% 	NT:			K x M x N matrix of trajectories, where K is number of trials,
% 				M is dimension of space, N is number of time steps
% 	ID:			K vector of ids belonging to each trial in NT
% 	VARARGIN:	'colors': matrix of colors (deprecated)
% 				'symbol': plot symbol (deprecated)
% 				'lwidth': line width (deprecated)
% 				'msize': marker size (deprecated)
% 				'avg_trajectories': plot average trajectories
% 				'fhandle': figure handle. If not given we plot into whatever is frontmost.
% 				'propfunction': function that gets called to determine plot object properties.
% 						Should take an id/event code as argument and return a structure with
% 						plot properties.
% 				'use_patchline': true|{false} use patchline function instead of normal plot
% 						to allow transparent lines
% [ALLH] = PLOTNT(NT,ID,VARARGIN)
%
% @author Carlos Vargas-Irwin, Jonas Zimmermann
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.
if nargout > 1
	error('This function has only one output!');
end
if nargin < 2 || isempty(id)
	id = 1:size(NT, 1);
end
kw = getopt(['colors symbol=''o'' lwidth=2 msize=5 fhandle ahandle add_legend=1 '...
' propfunction average_only=0 avg_trajectories=0 smooth_trajectories smooth_kernel=''gauss'' ' ...
' smooth_steps=5 use_patchline curr_index plot_history set_user_data extra_data legend_location=''southeastoutside'''],varargin);
msize = kw.msize;
lwidth = kw.lwidth;
symbol = kw.symbol;

if ~isempty(kw.add_legend) && ((ischar(kw.add_legend) && strcmpi(kw.add_legend, 'yes')) || kw.add_legend)
	kw.add_legend = 1;
else
	kw.add_legend = 0;
end

kw.use_patchline = parseBinArg(kw.use_patchline, false);
kw.set_user_data = parseBinArg(kw.set_user_data, true);
if kw.use_patchline && ~exist('patchline', 'file')
	kw.use_patchline = false;
	warning('SSIMSToolBox:plotNT:patchlineNotInPath', 'You asked me to use ''patchline'' for plotting, however this function is not on the MATLAB path. Will use normal plot functions instead.');
end

if isempty(kw.colors)
	nuid = length(unique(id));
	if nuid <= 1
		nuid = 2;
	end
	cval = lines(nuid);
else
	cval = kw.colors;
end

if ~isempty(kw.fhandle) && (isempty(kw.ahandle) || ~ishandle(kw.ahandle))
	set(0, 'CurrentFigure', kw.fhandle);
	kw.ahandle = get(kw.fhandle, 'currentAxes');
	if isempty(kw.ahandle)
		kw.ahandle = axes('parent', kw.fhandle);
	end
end
if isempty(kw.ahandle) || ~ishandle(kw.ahandle)
	kw.ahandle = gca;
end
hold(kw.ahandle, 'on');

if ~isempty(kw.propfunction)
	if ischar(kw.propfunction)
		propf = str2func(kw.propfunction);
	elseif isa(kw.propfunction, 'function_handle')
		propf = kw.propfunction;
	else
		propf=[];
	end
else
	propf = [];
end

kw.smooth_trajectories = parseBinArg(kw.smooth_trajectories, false);
if kw.smooth_trajectories
	[kw.smooth_kernel, kw.smooth_steps] = calculateSmoothKernel(kw.smooth_kernel, kw.smooth_steps);
end

idSize = size(id);
if idSize(1) > 1
    id = id';
end

[~, ndim, tsteps] =  size(NT);
ugid = unique(id);
if tsteps > 1
    allh = gobjects(1,length(id));
else
    allh = gobjects(1,length(ugid));
end

if ~isempty(kw.plot_history) && ~isempty(kw.curr_index)
	if kw.curr_index > tsteps
		kw.curr_index = tsteps;
	end
	start_i = max(1, kw.curr_index - kw.plot_history );
	plot_times = start_i:kw.curr_index ;
else
	plot_times = 1:tsteps;
end


classcount = 0;
for gon=(ugid)

	classcount = classcount + 1;
	f = (id == gon);
	h = [];
    if (kw.average_only)
        dat=squeeze(NT(f,:,:));
        md=mean(dat);
        sd=std(dat);
        [x, y, z]  = ellipsoid(md(1),md(2),md(3),sd(1)/2,sd(2)/2,sd(3)/2,30);
        h = surfl(x, y, z, 'parent', kw.ahandle);
    else
	if ndim==2
		x = squeezeKeepFirst((NT(f, 1, :)));
		y = squeezeKeepFirst((NT(f, 2, :)));
		if kw.avg_trajectories
			if kw.smooth_trajectories
				x = convWithNorm(mean(x, 2), kw.smooth_kernel);
				y = convWithNorm(mean(y, 2), kw.smooth_kernel);
			else
				x = mean(x, 2);
				y = mean(y, 2);
			end
		else
			if kw.smooth_trajectories
				x = convWithNorm(x, kw.smooth_kernel);
				y = convWithNorm(y, kw.smooth_kernel);
			end
		end
		if (tsteps > 1) && kw.use_patchline
			x = x(plot_times, :);
			y = y(plot_times, :);
			for i_h = 1:size(x, 2)
				h(i_h) = patchline(x(:,i_h), y(:,i_h), 'parent', kw.ahandle, 'facecolor', 'none');
			end
		else
			h = plot(x(plot_times, :), y(plot_times, :), symbol, 'parent', kw.ahandle);
		end
	else
		x = squeezeKeepFirst((NT(f, 1, :)));
		y = squeezeKeepFirst((NT(f, 2, :)));
		z = squeezeKeepFirst((NT(f, 3, :)));


		if kw.avg_trajectories
			if kw.smooth_trajectories
				x = convWithNorm(mean(x, 2), kw.smooth_kernel);
				y = convWithNorm(mean(y, 2), kw.smooth_kernel);
				z = convWithNorm(mean(z, 2), kw.smooth_kernel);
			else
				x = mean(x, 2);
				y = mean(y, 2);
				z = mean(z, 2);
            end
               	
		else
			if kw.smooth_trajectories
				x = convWithNorm(x, kw.smooth_kernel);
				y = convWithNorm(y, kw.smooth_kernel);
				z = convWithNorm(z, kw.smooth_kernel);
			end
			%h = plot3(x, y, z, symbol, 'parent', kw.ahandle);
		end
		if (tsteps > 1) && kw.use_patchline
			x = x(plot_times, :);
			y = y(plot_times, :);
			z = z(plot_times, :);

			for i_h = 1:size(x, 2)
				h(i_h) = patchline(x(:,i_h), y(:,i_h), z(:,i_h), 'parent', kw.ahandle, 'facecolor', 'none');
			end
		else
			h = plot3(x(plot_times, :), y(plot_times, :), z(plot_times, :), symbol, 'parent', kw.ahandle);
        end
        
           
	end
	propvec = isprop(h, 'Color');
	if any(propvec)
		set(h(propvec), 'color', cval(mod(classcount-1, size(cval, 1))+1,:));
	end
	propvec = isprop(h, 'linewidth');
	if any(propvec)
		set(h(propvec), 'linewidth', lwidth);
	end
	propvec = isprop(h, 'markersize');
	if any(propvec)
		set(h(propvec), 'markersize', msize);
	end
	propvec = isprop(h, 'markerfacecolor');
	if any(propvec)
		set(h(propvec), 'markerfacecolor', cval(mod(classcount-1, size(cval, 1))+1,:));
	end
	propvec = isprop(h, 'markeredgecolor');
	if any(propvec)
		set(h(propvec), 'markeredgecolor',cval(mod(classcount-1, size(cval, 1))+1,:));
	end


    end
    if ~isempty(propf)
		% if there is an extra_data struct, also pass it to the property function
		if isempty(kw.extra_data)
			applyPropertiesToHandle(h, propf(gon));
		else
 			applyPropertiesToHandle(h, propf(gon, 'extra_data', kw.extra_data(f)));
		end
		if kw.set_user_data
			if length(h) > 1
				set(h(2:end), 'UserData', []);
			end
			if isstruct(get(h(1), 'UserData'))
				ud = get(h(1), 'UserData');
				ud.n_tr = sum(f);
				set(h(1), 'UserData', ud);
			end
		else
			set(h, 'UserData', []);
		end
		if kw.avg_trajectories && (tsteps > 1) && (isempty(kw.plot_history) || (kw.plot_history > 0) )
			set(h, 'marker', 'none');
		end
	end

    if tsteps > 1
    	allh(f) = h;
    else
        allh(classcount) = h;
    end
end
axis off
if length(allh)==1
    handle_with_userdata = ~isempty(get(allh, 'userdata'));
else
    handle_with_userdata = cellfun(@(x)~isempty(x),(get(allh, 'userdata')));
end
if kw.add_legend && any(handle_with_userdata)
	ud = get(allh(handle_with_userdata), 'userdata');
	uh = allh(handle_with_userdata);
	if length(handle_with_userdata) == 1
		uuh = 1;
		uud = ud;
		if isstruct(uud) && isfield(uud, 'label')
			uud = uud.label;
		end
	else
		uds = getLabelData(ud);
		[uud, uuh] = unique(uds);

	end
	if iscellstr(uud) || ischar(uud)
		legend(uh(uuh), uud, 'Location', kw.legend_location);
	end


end
if ndim>2
	view(2)
end
axis tight
if nargout > 0
	varargout{1} = allh;
end

end

function x = squeezeKeepFirst(x)
	x = reshape(x, size(x, 1), size(x, 3), [])';
end

function s = getLabelData(ud)
	s = cell(numel(ud), 1);
	for ii = 1:numel(ud)
		if ischar(ud{ii})
			s{ii} = ud{ii};
		elseif isstruct(ud{ii}) && isfield(ud{ii}, 'label')
			s{ii} = ud{ii}.label;
		else
			s{ii} = '';
		end
	end
end
