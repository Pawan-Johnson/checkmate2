# CheckMATE
CheckMATE (Check Models At Terascale Energies) is a program package which accepts simulated event files in many formats for any given model. The program then determines whether the model is excluded or not at 95% C.L. by comparing to many recent experimental analyses. Furthermore the program can calculate confidence limits and provide detailed information about signal regions of interest. It is simple to use and the program structure allows for easy extensions to upcoming LHC results in the future.

More details can be found at: https://checkmate.hepforge.org/

<!-- TODO: Add Citation -->

## Requirements
1. Python: The current version is compatible with both Python2 (2.7.4 or later) and Python3 (possibly >= 3.4).<br>
   The  'future' library is required (`pip install future`). You can control which Python version is used at configuring with `--with-python=/path/to/python`. Other than that your default `python` will be used. 
   
   Note: that the features like multi-bin signal region and likelihoods require Python3. 
2. root: Anything later than version 6.20 should work fine. <br>
   Instructions to install or the precompiled root binaries can be obtained from: https://root.cern/install

3. Delphes: Version 3.5 or later https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/Releases

    Note: Should you need version 3.4 you need several changes in the code which you will figure out from compiler errors.
    
4. HEPMC2: https://gitlab.cern.ch:8443/hepmc/HepMC/-/tags
  
5. Pythia (optional): We currently only support 8.2 series. Installation instructions can be obtained from: https://pythia.org/
  
6. MadGraph5_aMC@NLO (optional): Supports Version >=2.7 should work just fine. With 3 series there are some hiccups but you should get the result anyway. Follow the instructions at: http://madgraph.phys.ucl.ac.be/ 
     
     
## Installation

This is a brief note on how to install CheckMATE from the master branch on https://github.com/CheckMATE2/checkmate2

To build CheckMate from source:

1. Obtain the Source:
   
    To get the source, download the release .tar.gz or .zip package in the release page:
    ```
    https://github.com/protocolbuffers/protobuf/releases/latest
    ```
    Alternatively you can get the source by cloning the repository using the `git clone` command
    ```git
    git clone https://github.com/CheckMATE2/checkmate2
    ```
    Now `cd` into the directory containing the source code.
    ```bash
    cd checkmate2
    ```
2. Update the configuration files by using `autoreconf`:
   ```bash
   autoreconf
   ```
3. Configure the code for compilation:
   A sample configuration call would be
   ```bash
   ./configure --with-delphes=/path/to/delphes --with-hepmc=/path/to/hepmc --with-madgraph=/path/to/madgraph --with-pythia=/path/to/pythia
   ```  
4. Compile the code using `make`:
   ```bash
   make
   ```
    Optional: Use `make -j4` if you would like to compile on 4 cores.
  
5. To test the Installation (if you linked Pythia):
   ```bash
    cd bin
    ./CheckMATE 13tev_test.in
    ```
     
## New features

1. Multibin fits available in atlas_2010_14293 and atlas_1908_03122 (more to come):<br>
     For command line use "-mb simple/cls/full"; in input file "Multibin: simple/cls/full"
     simple: simple fitting of shape in several bins as defined in a search, reasonably fast;
     cls: calculates CLs using full likelihood provided by the collaboration, slower but still manageable
     full: calculates CLs and upper limits using full likelihood but quite slow (extremely slow for atlas_2010_14293 - running time in days). 
     The results are stored in the multibin_limits folder.
     
2. Neural nets implemented with ONNX Runtime, available in atlas_2106_09609 and atlas_2211_08028:<br>
     During configure use "--with-onnx=/path/to/onnxruntime"; see: https://onnxruntime.ai/; tested with versions 1.12.1 and 1.13.1 (linux-x64)
     
3. Boosted decision trees in atlas_2010_14293 brought by ROOT TMVA.
     

## Current

2023-01-11   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ multibin profile likelihood

2022-11-28   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ Python3 compatibility

2022-11-17   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ added atlas_2010_14293
    ~ statistical combinations
    ~ TMVA support

2022-11-04   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ providing master path for external analysis resources (like in LLP)
    ~ Delphes 3.5 required

2022-09-30   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ implementation of NN method in atlas_2106_09609

2021-10-06   Jong Soo Kim <jsk@th.physik.uni-bonn.de>
        
    ~ added added atlas_1502_05686, atlas_2103_11684, atlas_2004_10894, atlas_2106_09609
	atlas_1911_06660

2021-08-23   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_phys_2014_007

2021-08-20   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_1807_07447

2021-07-16   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_1911_12606

2021-03-14   Arran Freegard <acf1g14@soton.ac.uk>
        
    ~ added atlas_1403_5294 and atlas_higg_2013_03

2021-02-12   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_1908_03122

2020-11-03   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_conf_2018_041

2021-02-10   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_2101_01629
    
020-04-08   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ added atlas_conf_2020_048

2020-09-27   Manimala Chakraborti, Ipsita Saha

    ~ added atlas_1803_02762

2020-09-27   Manimala Chakraborti, Ipsita Saha

    ~ added atlas_conf_2019_020

2020-09-29   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_2004_14060

2020-04-08   Krzysztof Rolbiecki <krolb@fuw.edu.pl>

    ~ removed python Root check due to a change in naming conventions (python->pyroot) around v6.20 

2019-11-28   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_conf_2019_040

2019-11-22   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_1909_08457

2019-10-03   Krzysztof Rolbiecki <krolb@fuw.edu.pl>
        
    ~ added atlas_1706_03731
