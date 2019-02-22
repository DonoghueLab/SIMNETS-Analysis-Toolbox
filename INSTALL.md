# Installation of the SSIMS and SIMNETS Toolbox #

## GETTING STARTED ## 
The toolbox is completely implemented in MATLAB. To get started, just add the
sub-folder `for_MATLAB` to your MATLAB path. However, certain basic MATLAB functions are rather slow. 
Therefore, core functionality has ALSO been implemented in C/C++ and can
be compiled as `mex` files. This will improve performance dramatically
(~2 orders of magnitude). See compilation instructions for your platform below.

    1. Add the sub-folder `for_MATLAB` to your MATLAB path 
    2. (optional) installation of compiled code for increased performance (see below)
    3. Use 'runSIMNETS.m' to analyze your own data (See 'SIMNETS_Live_tutorial.mlx' for guide)

## DEMO DATA AND TUTORIAL ##
The 'demo' folder contains two sample datasets and a script that can be run 
out-of-the-box without installation of the complied version. However, we *strongly* reccomend using
the complied version of the software. The live tutorial matlab script <SIMNETS_LiveTutorial.mlx> 
is set up to run with a demo dataset of simulated neurons or real M1 neurons.

                 'SIMNETS_LiveTutorial.mlx'

## FUNCTION DOCUMENTATION FOR MEX AND M FILES ##
MATLAB functions (and their `.mex` equivalents) are documented in the `.m` files.
See the 'demo' folder for a self-contained example with actual data.
The C/C++ implementation is rudimentarily documented. HTML and LaTeX representations can be 
created with [`doxygen'] (http://www.stack.nl/~dimitri/doxygen/), e.g. with 'make docs'


## HOW TO COMPILE MEX FILES INCLUDED IN SSIMS AND SIMNETS TOOLBOX ##
Several heavily used and computationally expensive functions have been rewritten
in C/C++ and are available as MEX files. MEX files need to be recompiled for the 
computer they are to be used on. MATLAB provides instructions and documentation 
on how to setup the compilation environment for different machines and operating
systems. You will need to install 2-3 programs to set-up the compilation environment.The specific
programs you will need to install will depend on your machine type (Mac or windows)and your 
Matlab version. 

All machine types will rely on the Armadillo [linear algebra library] (http://arma.sourceforge.net), 
which has to be installed first. 



## On Mac/linux: ##

## Mac OS X ##
This procedure has been tested to work with Mac OS X 10.9, 10.10, and 10.11, 10.12
with Xcode 6.2, 7.3, 7.3, and 9.2 respectively. MATLAB 2016a, 2017a, 2017a is confirmed to work
on all these systems, MATLAB 2015b has only been tested on 10.11.An example configuration file for 
MATLAB 2015b on Mac OS X 10.11 using `clang-omp` installed with [homebrew](http://brew.sh) has been included
with this toolbox (`mex_C++_maci64.xml`).For best results, use [openmp](http://openmp.org) for parallelized code
execution.

I - SETTING UP ENVIRONMENT AND COMPILING

1. INSTALL [Xcode](https://developer.apple.com/xcode/), from the [App Store](https://itunes.apple.com/us/app/xcode/id497799835?mt=12) or by download. Run it at least once to accept license agreement. Make sure that Xcode knows about its command line tools: open Xcode, open Xcode->Preferences, go to the 'Locations' pane, ensure that 'Xcode 7.3.1' (or whatever your version of Xcode is installed) is selected under 'Command Line Tools'.

    TROUBLE SHOOTING: Depending on your version of MATLAB, only some versions of Xcode are compatible (e.g. MATLAB 2014b only works with Xcode 4.6+ or 5.0+, but not 6.0+ or later. MATLAB 2015b supports Xcode 5.1+ or later, but [requires additional setup for 7.0+](http:/www.mathworks.com/support/sysreq/files/SystemRequirements Release2015b_SupportedCompilers.pdf)). MATLAB 2016a is confirmed to work with Xcode 6.0+ or 7.0+; MATLAB 2017a is confirmed to work with xcode 9.2 (and some earlier versions). If later on during compilation you receive an error message stating '`cstring.h` not found' or similar, try to run `xcode-select --install` from Terminal.

2. INSTALL [`homebrew`](http://brew.sh), the missing package manager for Mac OS.
    (Open `Terminal` and follow the instructions at the above website.)
    
3. INSTALL [Armadillo](http://arma.sourceforge.net), a fast linear algebra package for C++.
    Still in `Terminal`, enter

                brew install homebrew/science/armadillo

    This step may take a long time..
    
4. CHECK PATH IN TERMINAL: In `Terminal`, navigate to your home folder (`cd ~`). Enter `ls -al` to check if the file [`.matlab7rc.sh`] exists (look in right column; note that it is a hidden file, so normally will not appear in Finder). The name of the file is important:

        MATLAB *only* looks for `~/.matlab7rc.sh`!)
       
    ** If the file does not exist, copy the file `$toolbox/doc/matlab7rc.sh` to `~/.matlab7rc.sh`:

              cp $toolbox/doc/install_mac/matlab7rc.sh ~/.matlab7rc.sh

    where `$toolbox` corresponds to the path of this toolbox in your environment. Check for file again. In terminal, navigate to your home folder (`cd ~`) and enter `ls -al` to check if the file exists.

    ** If it does exist, open it in a text editor:

              open -a xCode .matlab7rc.sh

and change the line defining `LDPATH_PREFIX` in the `mac|maci|maci64)` section (around line 195) to read
    
            LDPATH_PREFIX='/usr/local/lib/gcc/6'
    
    (The path to `gcc` may change, as of writing, this is where it will be installed by `homebrew`).

5. RESTART: Open or restart MATLAB
6. IN MATLAB PROMPT: Change working directory to the toolbox path:

            cd('whatever your toolbox path is')
                    
            
7. IN MATLAB PROMPT: On the MATLAB prompt, enter

            copyfile localdefs_sample.mk localdefs.mk
        
            edit localdefs.mk
        

8. In MATLAB EDITOR: Make changes to the opened file your system. In particular, adjust `MATLABPATH` to point to your installation of MATLAB (~ line:15).

9. In `Terminal`, navigate to the toolbox folder and compile:

             cd 'whever your toolbox is'
             
             make
     

This should compile all necessary files. Compiled `*.mex*` files will be automatically placed in `for_MATLAB` in the toolbox package. Enjoy the speed!



II - SETTING UP PARALLEL PROCESSING
To use parallel processing, your compiler needs to support parallelized code. The `clang/LLVM` compiler bundled with Xcode <= 7.3.1 does not support OpenMP. However, it is easy to install a different version of `LLVM` and tell MATLAB to use it.

1. INSTALL (IN TERMINAL):  the current version of `llvm` from `brew` (from `Terminal` call):

            brew install llvm

2. COPY (IN MATLAB): copy `doc/install_mac/clang++_openmp_maci64.xml` in this package to your MATLAB preferred directory ( type `prefdir()` in matlab prompt to determine directory). Navigate to 'doc/install_mac/' and in matlab prompt, enter command:

            copyfile clang++_openmp_maci64.xml  '~/.matlab/'
         
  , where '~/.matlab/' is the preferred matlab path.

3. IN MATLAB: Tell MATLAB to use this configuration:

            mex -setup:'~/.matlab/R2017a/clang++_openmp_maci64.xml C++'
        
   where, '...' is the matlab path to the file moved in the previous step.

4. If you copied the `.matlab7rc.sh` startup script in step 4 (part I) above, continue with step 6. If you made changes to a
previous file, you will have to add a folder to `LDPATH_PREFIX`; this line (around line:195) should now read:

        LDPATH_PREFIX='/usr/local/lib/gcc/6:/usr/local/opt/llvm/lib'
        
5. RESTART MATLAB
6. IN MATLAB (EDITOR): Edit your `localdefs.mk` file to uncomment the options after the `## Using parallelization techniques.`  around line:30 to line:41. Line: 32 and Line:39 are commments (do not uncomment).
7. TERMINAL: In `Terminal`, navigate to the `toolbox` folder and enter

            make clean
            clean
8. RUN: In MATLAB, run `SSIMSToolboxTestParallel`. If it was compiled correctly it will display the number of parallel threads.


Note that an upgrade to macOS 12 may break things: linking against `libomp` included with `llvm` will crash MATLAB. Instead, use the build configuration in `clang++_openmp_macOS12_maci64.xml` and make changes to `localdefs.mk`, as shown in the sample file. Furthermore, `llvm` may have to be reinstalled (`brew reinstall llvm`). With questions regarding the compilation step, <dr.jonas@brown.edu> may be able to help.

## Linux ##
For now, you're on your own. Ensure the dependencies (`armadillo`, `clang` or another
parallelizable, OpenMP-compatible compiler) are satisfied. Copy `localdefs_sample.mk` to
`localdefs.mk` and make appropriate changes. Run

make all




## ON WINDOWS ##
I - SETTING UP ENVIRONMENT AND COMPILING

1. DOWNLOAD the [Armadillo library](http://arma.sourceforge.net/download.html)
      and unpack the archive. Rename the unpacked folder to `armadillo` and move
      it into the `lib` folder of this toolbox.

2.CHECK LIBRARYS: Make sure armadillo is correctly configured to use MATLAB's LAPACK
      and BLAS libraries. Armadillo's configuration file is at
      `armadillo\include\armadillo_bits\config.hpp`.

        TROUBLE SHOOTING TIPS: With our installation (MATLAB 2015b 64bit),
        line 69 for version 7.400.2 had to be uncommented to read

            #define ARMA_BLAS_LONG_LONG

        Also, line 62 had to be commented out:

            //#define ARMA_BLAS_UNDERSCORE

3.CHECK INSTALLATION: For compilation, we have used MinGW64 as it is available through MATLAB's
      AddOn system. Make sure it is installed correctly and mex uses it (to change,
      execute `mex -setup C++` and `mex -setup C` from the MATLAB prompt).
      Within MATLAB, navigate to the SSIMS toolbox folder. Adjust the

            build_files_win.m

      script to match your system's circumstances; however, if you use `armadillo`
      and 'MinGW64' as described, you should not have to change anything.
      Run `build_files_win.m` to compile. Build products will be automatically
      copied into the `for_MATLAB` folder.


Currently, parallelization using the [OpenMP API](http://openmp.org/wp/) has
not been tested under Windows. If you have compiled this toolbox successfully
with parallelization, please let one of the following people know:carlos Vargas-irwin, <carlos_vargas_irwin@brown.edu>;Jacqueline Hynes, <jacqueline_hynes@brown.edu>; or Jonas Zimmermann, <dr.jonas@brown.edu>. 



