
# SSIM and SIMNETS Analysis Toolbox

<p align="center"> <img src="images/SIMILARITY MAPS_schematics-03.png" alt="Fig 1. SIMNETS" class="inline" width="800" height="300"/>  </p>

### Quick-Guide ####
The purpose of this toolbox is to generate dimensionality-reduced **S**pike train **Sim**ilarity (SSIM) and Neuron **Sim**ilarity **Net**work (SIMNET) Maps from spike train data, facilitating visualization and further analysis. 

To get started: 
1. Download the tool box from Github: [Toolbox Repository](https://github.com/DonoghueLab/SSIM-and-SIMNETS-Analysis-Toolbox) 

2. Add toolbox to your MATLAB path.

3. **Optional, but recommended** Installation of complied Matlab code for significantly(!) increased performance see 'install.md'. 

4. Open **SIMNETS_Live_tutorial.mlx'** in Matlab: for detailed guidance on how to use SIMNETS (dim-reduced Neuron Similarity Map) and all of its sub-functions with two different demo datasets. 

5. Open **SSIMS_democenter_out.m** in Matlab: for guidance on using SSIMS (dim-reduced Ensemble Activity Spiketrain Simliarty Maps) with a single demo dataset.  

6. For more details on both methods, see pre-print [2], publication [3], or our webpage: [Donoghue Lab - Github Page - Analysis Toolbox](https://donoghuelab.github.io/SSIM-and-SIMNETS-Analysis-Toolbox/) 



#### INSTALLATION OF C/C++ OPTIMIZED MATLAB VERSION OF TOOLBOX (optional): ####
The toolbox is completely implemented in MATLAB. To get started, just add the sub-folder `for_MATLAB` to your MATLAB path. However, certain basic MATLAB functions are rather slow. Therefore, core functionality has ALSO been implemented in C/C++ and can be compiled as `mex` files. This will improve performance dramatically (~2 orders of magnitude). 

#### CITATION OF WORK: ####
Please cite the DOI for the SIMNETS pre-print and the DOI for the Software Repository (doi: 'assignment pending') when using this software and/or this analysis framework for analyzing your own data. 

[1] SSIM and SIMNETS Analysis Toolbox: DOI: ( 'assignment pending')  

[2] bioRxiv Pre-print: Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin (2018). **SIMNETS: a computationally efficient and scalable framework for identifying networks of functionally similar neurons** . DOI: https://doi.org/10.1101/463364) 

[3] Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J. (2015).  **Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014).** 

This toolbox integrates two algorithms to achieve the dimensionality-reduced Similarity Maps:

[4] Victor, J D and K P Purpura (1996). **Nature and precision of temporal coding in visual cortex: a metric-space analysis”. In: Journal of Neurophysiology** 76.2, pp. 1310–26.

[6]Van der Maaten, Laurens J P and Geoffrey E Hinton (Nov. 2008). **Visualizing High-Dimensional Data Using t-SNE**. In: Journal of Machine Learning Research 9, pp. 2579–2605.

#### Version history ####
------------------------

*   4.0.0:  22 Feb 2019 . @JBHynes
  Added SIMNETS (NEURON SIMILARITY)TOOLBOX, updated build/install instructions. Added toolbox to Github. First public release of SIMNETS. 
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



## QUESTIONS: 
@author Jacqueline Hynes. Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.

Questions? Contact <Carlos_Vargas_irwin@brown.edu> or <Jacqueline_Hynes@Brown.edu>
We are happy to help with any trouble shooting and provide guidance on how to best analyze/interpret your own data.
