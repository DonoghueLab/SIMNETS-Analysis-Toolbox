function [h] = plotMeanNT(NT,id,symbol,lwidth,msize,varargin)
% plotMeanNT(NT,id,symbol,lwidth,msize)
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.


if(nargin == 5)
    jval = jet(length(unique(id)));
    %jval = jval(1:60,:);
    skip =  floor(length(jval)/ (length(unique(id))) );
    cval(1,:) = [0 0 0];
    for c = 1:length(unique(id))
        cval(c,:) = jval(skip*c,:);
    end
else
    cval = varargin{1};
end


[trials, dim, timesteps] =  size(NT);

classcount = 0;
for gon=(unique(id))
    
    classcount = classcount + 1;
    
    f = find([ id ]== gon);
    
    if timesteps == 1;
        h = plot3(mean(NT(f,1)),mean(NT(f,2)),mean(NT(f,3)),symbol);
    else
        h = plot3(mean(reshape(NT(f,1,:),length(f),timesteps)),mean(reshape(NT(f,2,:),length(f),timesteps)),mean(reshape(NT(f,3,:),length(f),timesteps)),symbol);
        
    end
    
    
    %set(h,'linewidth',1,'markersize',15,'color','b','markerfacecolor','b','markeredgecolor','y')
    set(h,'linewidth',lwidth,'markersize',msize,'color',cval(classcount,:),'markerfacecolor',cval(classcount,:),'markeredgecolor',cval(classcount,:));
    hold on
    
    
end
axis off

set(gca, 'fontsize', 15)
% legend('TC Pre.','TC Pow.','Disk key','Disk Pow.')
axis tight
