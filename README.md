# **NAME**

Good_Turing.R - Calculate probability of unobserved species for molecular dynamics RMSD matrices.

# **DESCRIPTION**

Good_Turing.R is an [R](http://www.r-project.org/) script which applies Good-Turing statistics on
RMSD matrices of macromolecular dynamics simulation trajectories, with the purpose of quantifying
convergence and/or sufficient sampling. It is available for all platforms R runs on, including 
Windows, Mac and *nix. For more info on availability see [here]
(http://cran.r-project.org/doc/FAQ/R-FAQ.html#What-machines-does-R-run-on_003f).

The following steps are required in order to run the script :

### 1. Obtain and install R

If you are on a windows machine follow the link provided above to obtain the R installer for your platform.

If you are on a Linux/GNU machine you have the option of compiling R from the source, or using pre-compiled
binaries for your OS. For Debian like systems something like : `sudo apt-get install r-base` should do the trick.

### 2. Launch R

Windows users should be able to launch R through Start > All Programs > R > R executable or a  relevant link.

Linux/GNU users should be able to launch R through their terminal emulator of choice simply by typing `R`.

### 3. Launch the script

Once in the R environment users should switch to the directory where the script is located using the R command
`setwd()`, enclosing the path in the parentheses with single or double quotes and use the command 
`source('Good_Turing.R')`. Alternatively the `source()` command can be used without switcing directories but the
relative path to the script must be provided in the parentheses enclosed by single or double quotes.

### 4. First time run

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

### 5. The analysis

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

### 6. Troubleshooting
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


# **AUTHOR**

Good_Turing.R has been developed by Panagiotis Koukos, under the supervision of 
[Prof. Nicholas M. Glykos](http://utopia.duth.gr/~glykos/) at the 
[Department of Molecular Biology and Genetics](http://mbg.duth.gr/index.en.shtml)
of [Democritus University of Thrace](http://www.duth.gr/index.en.sxhtml).
