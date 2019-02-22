SSIMS and SIMNETS Toolbox V4.0
==================


Carlos Vargas-Irwin
Jacqueline Hynes
Jonas B. Zimmermann
David Brandman

Donoghue Lab, Brown University, 2012-2019.

The purpose of this toolbox is to generate dimensionality-reduced spike train similarity (SSIM) maps and Neuron Similarity Maps from spike train data, facilitating visualization and further analysis. The toolbox integrates two algorithms to achieve this:
* spike train distance metric described by Victor and Purpura (1996)
* t-Distributed Stochastic Neighbor Embedding (tSNE) (Van der Maaten and Hinton 2008).


## GETTING STARTED ## 

### Quick-Guide ###

    1. Add the sub-folder `for_MATLAB` to your MATLAB path 
    2. (optional) installation of compiled code for increased performance (see 'install.md' below)
    3. Use 'runSIMNETS.m' to analyze your own data (See 'SIMNETS_Live_tutorial.mlx' for guide)
    3. (optional) see the SIMNETS_live_tutorial for explanation of algorithm. 
  
The toolbox is completely implemented in MATLAB. To get started, just add the sub-folder `for_MATLAB` to your MATLAB path. However, certain basic MATLAB functions are rather slow. 
Therefore, core functionality has ALSO been implemented in C/C++ and can be compiled as `mex` files. This will improve performance dramatically (~2 orders of magnitude). See compilation instructions for your platform below.
                                                        --------------------------------------------------------------------------------------------------

## Cite ##
please cite the DOI for the SIMNETS paper and the DOI for the Software Repository (doi: xxxx) when using this software and/or this analysis framework for analyzing your own data. 
[1] SSIMS and SIMNETS toolbox: DOI: https://doi.org/  
[2] Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin (2018). "SIMNETS: a computationally efficient and            scalable framework for identifying networks of functionally similar neurons" . DOI: https://doi.org/10.1101/463364)       
[3] Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J. (2015).  "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
[4 ] Victor, Jonathan D. “Spike train metrics” Current opinion in neurobiology vol. 15,5 (2005): 585-92.
[5] .J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.

@author Jacqueline Hynes. Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
Questions? Contact <Carlos Vargas_irwin@brown.edu>  or <Jacqueline Hynes@Brown.edu>. 
We are happy to help with any trouble shooting or provide guidance on how to best analyze your own data.
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
[1] SSIMS and SIMNETS Neural Analysis Toolbox software DOI
[2] Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin (2018). "SIMNETS: a computationally efficient and            scalable framework for identifying networks of functionally similar neurons" . DOI: https://doi.org/10.1101/463364)       
[3] Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J. (2015).  "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
[4] Victor, J D and K P Purpura (1996). “Nature and precision of temporal coding in visual cortex: a metric-space analysis”. In: Journal of Neurophysiology 76.2, pp. 1310–26.
[5] .J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.
Van der Maaten, Laurens J P and Geoffrey E Hinton (Nov. 2008). “Visualizing High-Dimensional Data Using t-SNE”. In: Journal of Machine Learning Research 9, pp. 2579–2605.

