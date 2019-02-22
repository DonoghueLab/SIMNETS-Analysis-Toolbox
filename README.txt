SSIMS and SIMNETS Toolbox V4.0
==================


Carlos Vargas-Irwin
Jacqueline Hynes
Jonas B. Zimmermann
David Brandman

Donoghue Lab, Brown University, 2012-2019.

1. SIMNETS Live Matlab Tutorial Overview
How to use this tutorial: This live matlab script is part of the Spike train Similarity Space (SSIMS) and Similarty Network (SIMNETS) neural analysis tool-box. It's purpose is to provide the user with a 'how-to-guide' in running the main SIMNETS function (SSIMNETS.m) and to the illustrate the various operations performed by its sub-functions. Certain features of this tutorial will work best with Matlab Version >= R2018.  When analyzing your own data, we reccommend using the 'SSIMNETS.m' function (rather than this Live Tutorial Script) for increased speed. This following sections of this tutorial will cover: 1) SIMNETS Inputs, 2) Running SIMNETS, 3) SIMNETS Outputs, and 4) looking under-the-hood of SIMNETS. 
Getting Started: This code requires the 'SSIMS & SIMNETS toolbox' + the helper functions. The toolbox is completely implemented in MATLAB. To get started, just add the sub-folders `for_MATLAB` and 'demo' and 'matlab helpers' to your MATLAB path. However, as certain basic MATLAB functions are rather slow, the core functionality has also been implemented in C/C++ and can be compiled as `mex` files.  The C/C++  optimized version of the tool box will improve performance dramatically (2 orders of magnitude faster). Setting up the compilier environment will take about 10-15 minutes, but this is time saved down-the-line.See doc\INSTALL.md for instructions on installing the c++ optimzed version of the toolbox for your specific OS. 
This live script will run the SIMNETS algorithm using one of two possible demo datasets that are provided with this tool-box (uncomment preferred demo dataset in Section 2). 
                                                                  --------------------------------------------------------------------------------------------------
What is the SIMNETS algorithm? The SIMNETS algorithm is designed to compute a Neuron (functional) Similarity Map for a population of simultaneously recorded neurons, calculate the optimal number of sub-nets (i.e., neuron clusters) in the map, and then graphically display a 3D visualization of the Neuron Similarity Map containing the sub-nets. Note: this algorithm does not require knowledge of the neuron tuning functions (or even that they have any defined tuning functions) or knowledge of the information encoded in on individual trials, i.e., the specific internal or external covariate(s) that the neurons are responding to.

INPUTS: Although there are several input paramters that can be changes in the SIMNETS algorithm, 5 are required: 
Neural Data (2 inputs): the user can choose betwen one of two possible demo datasets that are provided with this tool-box (uncomment preferred demo dataset in Section 2). If you want to analyze your own neural dataset, you will need: 
A list of spike timestamps (in seconds) for each of the simultaneously recorded neurons (N). Ideally, the spike data will be spike-sorted into single-unit. SIMNETS will technically also work with the multi-unit data, however, you should factor this into your interpretation of the 'Neuron' similarity maps. 
A list of event timestamps (in seconds) for each of the events of interest (just start times). The events of interest can coorrespond to stimulus onset times, movement onset times, decision correlates, or even the onset of an interesting neural phenomenon (e.g., seziures activity, lfp gamma oscillation,  a particular phase of theta oscillatory activity). Events could also represent different time epochs within a trial that is composed on multiple stages, e.g., movement planning and execution. The events timestamps could also correspond to arbritrary events, e.g., every 1 second of a 30 minute period of sleep data, for example.  
Parameters (4 inputs): The user should vary these required parameters and other optional parameters to get a intuition for how their settings can effect the structure of the SIMNETS Neuronal Similarity Map and how their manipulation can faciliate hypothesis testing. Some of the parameters that will have the biggest impact on the Neuron Similarity maps include:                          
Win_len (in seconds): Spike train window length (in seconds)
start_offset (in seconds): Start of spike trains relative to the event of interest (e.g., 0.1 seconds before the start of movement event). 
 q : Temporal sensitivity of spike train comparison. Part of a cost-based edit-distance function [4] used to estimate similarity between spike trains that determine   the cost of 'shifting' a spike in time relative to insterting/deleting spikes. Low q values (e.g., q = 0) only register firing rate differences between the two spike trains being compared, whereas higher q values register differences in spike positions across the two spike trains being compared with increasing precision (e.g.,  q:100 = 10ms).
Perplexity: number of neurons used to define a 'local neighbour' in the neuron space during dimensionality reduction with tSNE [5].  The goal of tSNE is to preserve the local neighbourhood structure found in the high-dimensional space when projecting to lower-dimesions.

                                                             --------------------------------------------------------------------------------------------------

Cite: please cite the DOI for the SIMNETS paper and the DOI for the Software Repository (doi: xxxx) when using this software and/or this analysis framework for analyzing your own data. 
[1] SSIMS and SIMNETS toolbox: DOI: https://doi.org/  
[2] Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin (2018). "SIMNETS: a computationally efficient and            scalable framework for identifying networks of functionally similar neurons" . DOI: https://doi.org/10.1101/463364)       
[3] Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J. (2015).  "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
[4 ] Victor, Jonathan D. “Spike train metrics” Current opinion in neurobiology vol. 15,5 (2005): 585-92.
[5] .J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.

@author Jacqueline Hynes. Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
Questions? Contact <Carlos Vargas_irwin@brown.edu>  or <Jacqueline Hynes@Brown.edu>. We are happy to help with any trouble shooting or provide guidance on how to best analyze your own data.
Version history
---------------
#   4.0.0  22 Feb 2019
*   3.0.10: 3 November 2016
  Add build instructions for macOS 12 and MATLAB 2016b
*   3.0.9:  15 September 2016
  Major overhaul of the toolbox structure.
  Removed legacy functions, improved function signatures
  Add example with real data
  Improved installation instructions
  This is a pre-release to test functionality before wider distribution
*   3.0:    11 May 2016
	Rewrite of most of the toolbox. We now use armadillo for linear algebra functions.
	There are also efficient functions to extract spike trains in windows, based
	on custom C++ classes efficiently handling spike trains.
	Build instructions for Windows greatly improved
*   2.2:    17 Novemeber 2014
    First public release.


References
----------

Van der Maaten, Laurens J P and Geoffrey E Hinton (Nov. 2008). “Visualizing High-Dimensional Data Using t-SNE”. In: Journal of Machine Learning Research 9, pp. 2579–2605.

Victor, J D and K P Purpura (1996). “Nature and precision of temporal coding in visual cortex: a metric-space analysis”. In: Journal of Neurophysiology 76.2, pp. 1310–26.
