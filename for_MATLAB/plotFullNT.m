function [allh] = plotFullNT(NT,id,symbol,lwidth,msize,varargin)
% plotMeanNT(NT,trialid,symbol,lwidth,msize, *optionalcolorscheme*)


if(nargin == 5)
    jval = jet;
    jval = jval(1:60,:);
    skip =  floor(length(jval)/ (length(unique(id))-1) );
    cval(1,:) = [0 0 0]
    for c = 1:length(unique(id))-1
        cval(c+1,:) = jval(skip*c,:);
    end
else
    cval = varargin{1};
end


allh = [];


[trials  dim timesteps ] =  size(NT)

classcount = 0;
for gon=(unique(id))

    classcount = classcount + 1;
    
        f = find([ id ]== gon);   
        
        for n = 1:length(f)
            
            
%            h = plot3((NT(f(n),:,1)),(NT(f(n),:,2)),(NT(f(n),:,3)),symbol);
            
if timesteps == 1;
            h = plot3((NT(f(n),1)),(NT(f(n),2)),(NT(f(n),3)),symbol);
else 
            h = plot3(reshape(NT(f(n),1,:),1,timesteps),reshape(NT(f(n),2,:),1,timesteps),reshape(NT(f(n),3,:),1,timesteps),symbol);

end


            %set(h,'linewidth',1,'markersize',15,'color','b','markerfacecolor','b','markeredgecolor','y')
            set(h,'linewidth',lwidth,'markersize',msize,'color',cval(classcount,:),'markerfacecolor',cval(classcount,:),'markeredgecolor',cval(classcount,:));
            hold on
            
             allh = [allh h];
        end
             
end
axis off

set(gca, 'fontsize', 15)
% legend('TC Pre.','TC Pow.','Disk key','Disk Pow.')
axis tight
