function [ cindex centroids] = selectSignificanceClusterNumber( silh, shufCI, crange,NSPACEcluster)

% Overview: This function prompts the user to select a new cluster number. 
% Previous k-means cluster number was determined to be non-significant.
% User needs to select new cluster number in plot, or continue with the
% non-significant cluster number.
%
% INPUT: silhoutte values for NS map; the CI of the null Silhouette
% distribution; range of tested cluster numbers; NSmap
%
% OUTPUT: new cluster index values and their centroids
%
%
% Validated using MatLab R2018a
% Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin, 
% SIMNETS: a computationally efficient and scalable framework for identifying networks of functionally similar neurons
% 
% Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J
% "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
% 
% @author Jacqueline Hynes
% Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
% Questions? Contact Jacqueline Hynes@Brown.edu


%% 1) PLOT: silhouette values
    h1=figure; 
    ax(1) = jbfill(1:max(crange) ,shufCI(2,:) ,shufCI(1,:) ,1:max(crange));
    ax(2) = plot(1:max(crange), silh, 'k.-', 'linewidth', 2);
    plot(1:max(crange), silh, 'r.', 'markersize', 30);
    box off; xlabel('Cluster numbers');
    ylabel('Avg. Silhouette value');
    set(gca,'fontsize',14)
    set(gcf, 'position',[354   598   867   493]);

    title('Warning: non-significant cluster number. Use cursor to select new cluster value', 'color', 'r', 'fontweight', 'bold'); 
    warning('Spurious Optimal cluster number: confirm cluster number');
    figure(h1)
     
% Select a new cluster number before the timer runs out or else, continue
% without selecting a cluster number
    
    Input_timeout = 30; 
    t = timer('ExecutionMode','singleShot','StartDelay',5,'TimerFcn',@(~,~)close(h1)); start(t)

    [x, y, button] = ginput_timeout(1); % Wait for user to use cursor to select pixels
    close(1); 
    newK = round(x);
    if all(x) & all(y) & all(button) 
        
        newK = 1;
        disp('NO USER INPUT... no clusters selected')
        [cindex ] = kmeans(NSPACEcluster, newK,'Replicates',100, 'Display', 'off');
        
    elseif mod(x, round(x))>0
        
        fprintf('YOU MISSED... but we will round to %d\n', newK); 
        [cindex ] = kmeans(NSPACEcluster, newK,'Replicates',100, 'Display', 'off');
        delete(t)
    elseif mod(x, round(x))==0
        fprintf('GOOD AIM...new cluster number is %d\n', newK);
        [cindex ] = kmeans(NSPACEcluster, newK,'Replicates',100, 'Display', 'off');
        delete(t)
    end

      
end 
 
