# **CONTENTS**

1. [Installation Instructions](https://github.com/pkoukos/GoodTuringMD#1-install)
2. [Name](https://github.com/pkoukos/GoodTuringMD#2-name)
3. [Description](https://github.com/pkoukos/GoodTuringMD#3-description)
  1. [Extremely quick start guide](https://github.com/pkoukos/GoodTuringMD#i-extremely-quick-start-guide)
  2. [Obtain and install R](https://github.com/pkoukos/GoodTuringMD#ii-obtain-and-install-r)
  3. [Launch R](https://github.com/pkoukos/GoodTuringMD#iii-launch-r)
  4. [Launch the script](https://github.com/pkoukos/GoodTuringMD#iv-launch-the-script)
  5. [First time run](https://github.com/pkoukos/GoodTuringMD#v-first-time-run)
  6. [The analysis](https://github.com/pkoukos/GoodTuringMD#vi-the-analysis)
  7. [Troubleshooting](https://github.com/pkoukos/GoodTuringMD#vii-troubleshooting)
4. [Support](https://github.com/pkoukos/GoodTuringMD#4-support)
5. [Author](https://github.com/pkoukos/GoodTuringMD#author)


# **1. INSTALL**

Download [this](https://raw.githubusercontent.com/pkoukos/GoodTuringMD/master/Good_Turing.R) file and
follow the instructions [below](https://github.com/pkoukos/GoodTuringMD#0-extremely-quick-start-guide).

# **2. NAME**

Good_Turing.R - Application of Good-Turing statistics to quantify convergence of biomolecular simulations.

# **3. DESCRIPTION**

**A paper describing the method and the type of results obtained from
its application is available from [JCIM](http://dx.doi.org/10.1021/ci4005817).**

Good_Turing.R is an [R](http://www.r-project.org/) script which applies Good-Turing statistics on
RMSD matrices of macromolecular dynamics simulation trajectories, with the purpose of quantifying
convergence and/or sufficient sampling. It is available for all platforms R runs on, including 
Windows, Mac and *nix. For more info on R availability see [here]
(http://cran.r-project.org/doc/FAQ/R-FAQ.html#What-machines-does-R-run-on_003f).

The following steps are required in order to run the script :

### i. Extremely quick start guide 

For a more detailed overview skip this section and begin your reading at ii. If you feel confident read on.

If you already have the R package installed and you work with unix (GNU/Linux-MacOSX) then
just place the Good_Turing.R program to the directory containing your RMSD matrix and then
from the unix shell do the following :

1. Go to the directory containing the matrix plus the Good_Turing.R files.
2. Run R from the console by typing 'R'
3. Type : source('Good_Turing.R')
4. For the very first run and if any needed packages are missing you will be given the opportunity
   to install them on-the-fly. Do install them.
5. Type the name of the file containing the RMSD matrix when asked.
6. Hold your breath.

The expected format of the RMSD matrix is described in section vii below.

### ii. Obtain and install R

If you are on a windows machine follow the link provided above to obtain the R installer for your platform.

If you are on a Linux/GNU machine you have the option of compiling R from the source, or using pre-compiled
binaries for your OS. For Debian like systems something like : `sudo apt-get install r-base` should do the trick.

### iii. Launch R

Windows users should be able to launch R through Start > All Programs > R > R executable or a  relevant link.

Linux/GNU users should be able to launch R through their terminal emulator of choice simply by typing `R`.

### iv. Launch the script

Once in the R environment users should switch to the directory where the script is located using the R command
`setwd()`, enclosing the path in the parentheses with single or double quotes and use the command 
`source('Good_Turing.R')`. Alternatively the `source()` command can be used without switcing directories but the
relative path to the script must be provided in the parentheses enclosed by single or double quotes.

### v. First time run

Every time the script is executed it performs a check on the system library in order to determine whether its
dependencies are met. The two dependencies(v0.1) of the program are the packages :
* fastcluster
* minpack.lm

Strictly speaking the first is not a dependency since all it does is replace R's stock `hclust()` function for
hierarchical clustering with a much faster version. The second is necessary due to the need for a non linear
regression to be performed reliably. The function `nlsLM()` of the package is used to that end.

The user has the choice of installing the packages on his/her own or launching the script and let it determine
if the packages are installed. If they are not found then the user will be prompted to install them and after
their successful installation the execution of the program will proceed normally.

Users wishing to install the packages on their own can do so using the following command :
    
`install.packages(c('fastcluster', 'minpack.lm'))`
    
The previous command does not specify the location of the library nor which repository to use for the
installation so the user will be asked to specify these after issuing the command. Note of importance
for GNU/Linux users : If you would like to install the packages to a system wide location such as /usr/local/
should launch R with administrative privileges(ie `sudo R`), otherwise a personal library will be used.

### vi. The analysis

After the program has been launched and provided that the dependencies are met the program proceeds to evaluate
the RMSD matrix for compatibility(see the following section for details). After it has been determined that the
matrix is compatible with the program, the actual analysis begins. The progress of the analysis is reported in
STDOUT in the form of text messages. Successful completion of the program will result in the creation of two
files as well as a tar.

The two files are named :
* `good_turing.<filename>.prob.of.unobserved_vs_RMSD.dat` and
* `good_turing.<filename>.prob.of.unobserved_vs_RMSD.eps` where `<filename>` is the name of the provided file

The second file is an encapsulated postscript plot of the data of the first file. These files are the main
result of the analysis and they detail the probability of unseen/unobserved species/molecular conformations
for a range of RMSD values.

The above files are not always produced. Whether they are produced or not depends on the data contained in
the files archived in the tar. These are :
* `good_turing.<filename>.max_rmsds.dat`
* `good_turing.<filename>.max_rmsds.eps`
* `good_turing.<filename>.max_of_mins.dat`
* `good_turing.<filename>.max_of_mins.eps`

The data in these files is analysed in order to determine whether the length of the trajectory suffices for
quantifying convergence.

### vii. Troubleshooting
The program has been extensively tested with matrices produced with trajectories of the PDB/PSF-DCD world,
however any matrix that meets the following criteria is acceptable :
    
    *****************************************************************
    **                                                             **
    ** This program expects as input a plain ASCII file containing **
    ** a (NxN) RMSD matrix with its origin at the upper left-hand  **
    ** corner :                                                    **
    **                                                             **
    **                   a11 a12 a13 ... a1N                       **
    **                   a21 a22 a23 ... a2N                       **
    **                   a31 a32 a33 ... a3N                       **
    **                   ...................                       **
    **                   aN1 aN2 aN3 ... aNN                       **
    **                                                             **
    ** For example, the following is a portion of a valid input :  **
    **                                                             **
    **                0.000 0.803 0.826 ... 2.138                  **
    **                0.803 0.000 0.689 ... 2.074                  **
    **                0.826 0.689 0.000 ... 2.065                  **
    **                ...........................                  **
    **                2.138 2.074 2.065 ... 0.000                  **
    **                                                             **
    *****************************************************************

Runtime errors associated with bad input are accompanied by a message explaining exactly what is wrong with the
data and the above message.

# **4. Support**

You can always submit your questions/bug reports etc. and receive mailing list support at the [google group](http://groups.google.com/group/good-turing-md "GoodTuringMD Mailing List") of the program.

# **5. AUTHOR**

Good_Turing.R has been developed by Panagiotis Koukos, under the supervision of 
[Prof. Nicholas M. Glykos](http://utopia.duth.gr/~glykos/) at the 
[Department of Molecular Biology and Genetics](http://mbg.duth.gr/index.en.shtml)
of [Democritus University of Thrace](http://www.duth.gr/index.en.sxhtml).
