
% Overview: runSIMNETS.m from the SSIMS and SIMNETS toolbox 
%
% See Matlab tutorial for demo and explanation of steps involved in
% SSIMNETS Algorithm. 
%
% 
%
% Validated using MatLab R2018a
%
% Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin, 
% SIMNETS: a computationally efficient and scalable framework for identifying networks of functionally similar neurons
% 
% Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J
% "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
% 
% @author Jacqueline Hynes
% Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
% Questions? Contact <carlos_vargas_irwin@brown.edu>; <Jacqueline
% Hynes@Brown.edu>

close all 
 
%% Load Data: 

  %---- load demo data '/demo' or your own data to run analysis--- %
        oad('SIMNETS_realData_centerOutM1.mat')                    % REAL M1 DATA
 
    spiketrains = {}
    events = []; 

%% PARAMETERS: 


    % Required Parameters
    start_offset = 0;   % Start of timewindow relative to events

    win_len =  ;       % Start of timewindow relative to events

    q =  ;             % temporal resolution of analysis (1/q = ms). 

    perp =   ;         % nearest neighbors in tSNE
    

    % Optional Parameters
    
    clusterdim = 10;           % dimensions for clustering

    displaydim =3;              % dimensions for visualization 

    iterations = 0;     % Turn ON/OFF statistical test by specifiying number of iterations (NB. This can be time-consuimg)

    corrType = 'pearson';      % 'Pearson' 'Kendall'  'Tau'

    drTechnique = 'tsne';      % Dimensionality Reduction Method: 

    crange = 10 ;               % cluster number paramater sweep range

    silhCI = 99;               % signifiance threshold for silhouette test (in percent)

  
%% Run SSIMNETS Function

tic 

[ SNETS ] = SSIMNETS(spiketrains, events+start_offset, win_len, q, perp, 'displaydim',displaydim ,'clusterdim', clusterdim, 'iterations',iterations, 'crange',crange, 'silhCI',silhCI, 'corrType',corrType, 'drTechnique', 'tsne' );

toc
 
    %--- SIMNETS PLOTS ---% 
    % 1a: Plot neurons with labeled with the k-means clusters indices 
            subplot(2,2,1)
            plotNT(SNETS.NSPACE,neuronID,'color', colsPD); 
            title('SIMNETS Neuron space - preferred directions')
            colorbar; colormap hsv; caxis([0 359 ]);
            hcb = colorbar;
            title(hcb,'pd')
         
    % 1b: Plot neurons with labeled with the k-means clusters indices 
            subplot(2,2,2) 
            plotNT(SNETS.NSPACE,SNETS.clusterindex)
            title('SIMNETS Neuron map - k-means clusters')
            legend({'cluster 1','cluster 2','cluster 3'},'location',  'bestoutside')
           

