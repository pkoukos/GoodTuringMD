################################################################################
##                                                                            ##
## Copyright (c) 2013, Panagiotis I. Koukos                                   ##
## All rights reserved.                                                       ##
##                                                                            ##
## Redistribution and use in source and binary forms, with or without         ##
## modification, are permitted provided that the following conditions         ##
## are met:                                                                   ##
##                                                                            ##
## 1. Redistributions of source code must retain the above copyright notice,  ##
##    this list of conditions and the following disclaimer.                   ##
##                                                                            ##
## 2. Redistributions in binary form must reproduce the above copyright       ##
##    notice, this list of conditions and the following disclaimer in the     ##
##    documentation and/or other materials provided with the distribution.    ##
##                                                                            ##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"##
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  ##
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ##
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE  ##
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR        ##
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF       ##
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS   ##
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN    ##
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    ##
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE ##
## POSSIBILITY OF SUCH DAMAGE.                                                ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
## Name        : Good_Turing.R                                                ##
## Version     : 0.1, September 2013                                          ##
## Input       : This program expects an RMSD matrix as input. The matrix must##
##               be square, ie it has to have the same number of lines and    ##
##               columns, and symmetrical. Columns have to be separated by    ##
##               whitespace (spaces, tabs, etc. )                             ##
## Output      : Successful completion of the program results in two files    ##
##               being produced :                                             ##
##               good_turing.<file.name>.prob.of.unobserved_vs_RMSD.dat and   ##
##               good_turing.<file.name>.prob.of.unobserved_vs_RMSD.eps       ##
##               where <file.name> is the full name of the provided file.     ##
##               as well as a tar containing 4 additional files. These are :  ##
##               good_turing.<file.name>.max_rmsds.dat                        ##
##               good_turing.<file.name>.max_rmsds.eps                        ##
##               good_turing.<file.name>.max_of_mins.dat                      ##
##               good_turing.<file.name>.max_of_mins.eps                      ##
##               For additional info about these files read the comments below##
## Description : The program loads the matrix and performs several checks on  ##
##               the data in order to determine whether it is suitable for the##
##               main analysis, which results in the creation of the two files##
##               described above. More detailed info in the comments below.   ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
## Dependencies.                                                              ##
## The program checks whether the required packages can be found in the local ##
## system, and if they are not then it installs them. Installation in the     ##
## default directories requires root priviliges for *nix systems.             ##
## 1. fastcluster : This package is not required but will significantly speed ##
##                  up the program, and as such is strongly recommended.      ##
## 2. minpack.lm  : This is the only true dependency of the program and is    ##
##                  necessary in order for a non linear curve fit to be       ##
##                  performed reliably.                                       ##
##                                                                            ##
################################################################################
if ('fastcluster' %in% rownames(installed.packages()) == FALSE &&
    'minpack.lm'  %in% rownames(installed.packages()) == FALSE) {
  install.packages(c('fastcluster', 'minpack.lm'))
} else if ('fastcluster' %in% rownames(installed.packages()) == FALSE) {
  install.packages('fastcluster')
} else if ('minpack.lm'  %in% rownames(installed.packages()) == FALSE) {
  install.packages('minpack.lm')
}

require(fastcluster)
require(minpack.lm)

################################################################################
##                                                                            ##
## Custom functions to be used during the execution of the program.           ##
##                                                                            ##
################################################################################

ComputeClusters <- function(my.matrix, rmsd.cutoff, rmsd.step) {
# Calculate probability of unobserved species for various RMSDs.
# Args :
#   my.matrix   : An RMSD matrix that will be used for the creation of a
#                 distance matrix, which will in turn be used for hierarchical
#                 clustering of the RMSD distances.
#   rmsd.cutoff : The 'height' at which the dendrogram created by hc are 'cut',
#                 thus determining how many values are grouped together.
#   rmsd.step   : The constant by which rmsd.cutoff is increased after each
#                 iteration of the while loop.
# Returns :
#   Two vectors containing the RMSDs and the probability of unobserved species
#   for each RMSD.
  hc <- hclust(as.dist(my.matrix), method="complete")
  results.x <- vector()
  results.y <- vector()
  j <- 1
  prob.unobserved <- 1
  while (prob.unobserved > 0) {
    clusters <- as.vector(cutree(hc, h=rmsd.cutoff))
    # Determine which clusters have frequency=1 and divide them with the sum of
    # the frequencies in order to obtain the probability of unobserved species.
    cluster.frequencies <- cbind(Frequency=sort(table(clusters)))
    occurrences.of.min  <- length(which(cluster.frequencies==1))
    prob.unobserved     <- occurrences.of.min / sum(cluster.frequencies)
    # For as long as the probability is greater than zero repeat the loop, but
    # stop as soon as it reaches zero.
    if (prob.unobserved > 0) {
      results.x[j] <- rmsd.cutoff
      results.y[j] <- prob.unobserved
    } else {
      break
    }

    rmsd.cutoff <- rmsd.cutoff + rmsd.step
    j <- j + 1
  }
  return(cbind(results.x, results.y))
}

DetermineRmsdStep <- function(my.matrix, init.rmsd, nofpoints) {
# Calculate the rmsd.step variable required for the above function.
# Args :
#   my.matrix : Same as for the above function.
#   init.rmsd : Same as for the rmsd.cutoff variable of the above function.
# Returns :
#   A numeric variable which is the rmsd.step to be used in the above function.
  hc <- hclust(as.dist(my.matrix), method="complete")

  occurrences.of.min <- 1
  lower.rmsd         <- 0
  upper.rmsd         <- 10
  # The first step for the calculation of the rmsd.step variable is determining
  # the lower RMSD value for which the probability of unobserved species equals
  # zero. This is achieved by calculating the probability of unobserved species
  # for an upper and a lower RMSD value, that approach each other until they
  # converge within 0.01.
  while (TRUE) {
    clusters            <- as.vector(cutree(hc, h=upper.rmsd))
    cluster.frequencies <- cbind(Frequency=sort(table(clusters)))
    occurrences.of.min  <- length(which(cluster.frequencies == 1))
    if (occurrences.of.min == 0) {
      previous.upper <- upper.rmsd
      if (lower.rmsd == 0) {
        upper.rmsd <- upper.rmsd / 2
      } else {
        upper.rmsd <- (lower.rmsd + upper.rmsd) / 2
      }
    } else {
      lower.rmsd <- upper.rmsd
      upper.rmsd <- previous.upper
    }

    if (abs(lower.rmsd - upper.rmsd) < 0.01) {
      break
    }
  }
  # After this value has been determined the minimum RMSD of the original matrix
  # is substracted from it, and it is divided by nofpoints. This is the number
  # of iterations the above function will perform.
  rmsd.step <- (upper.rmsd - init.rmsd) / nofpoints
  return(rmsd.step)
}

DetermineMinLength <- function(lengths){
# Determine which of the vectors in the list lengths has the smallest length.
# Args :
#   lengths : A list of vectors of variable lengths.
# Returns :
#   The vector with the smallest length.
  temp <- which(sapply(lengths, length)==min(sapply(lengths, length)))
  return(unlist(lengths[[temp[1]]]))
}

DetermineIfBinary <- function(filepath, max=1000) {
# Determine if a file is binary. Shamelessly copied from Spacedmans comment @
# http://stackoverflow.com/questions/16350164/
# Args :
#   filepath : The absolute or relative path of the provided file.
#   max      : Number of characters to read.
# Returns :
#   TRUE if a file is binary, FALSE otherwise.
  f <- file(filepath, "rb", raw=TRUE)
  b <- readBin(f, "int", max, size=1, signed=FALSE)
  close(f)
  return(max(b)>128)
}

DetermineFirstDiagonalMaxRmsd <- function(my.matrix, mat.length) {
# Determine the maximum RMSD of the first parallel to the diagonal of my.matrix.
# Args :
#   my.matrix  : An RMSD matrix.
#   mat.length : The number of lines of the matrix.
# Returns :
#   The max RMSD of the first parallel to the diagonal of the matrix.
  i        <- 1
  max.rmsd <- 0
  while (i < mat.length) {
    if (my.matrix[i, i + 1] > max.rmsd) {
        max.rmsd <- my.matrix[i, i + 1]
    }

    i <- i + 1
  }

  return(max.rmsd)
}

CreateSubmatrix <- function(my.matrix, sampling, max.rmsd.calc,
                            max.of.mins.calc) {
# Create a submatrix of the original and create the probability of unobserved
# species vs RMSD files.
# Args :
#   my.matrix        : An RMSD matrix.
#   sampling         : The factor by which my.matrix will be subsampled.
#   max.rmsd.calc    : Logical. If TRUE will calculate the max RMSD of the first
#                      diagonal of the subsampled matrices and return them
#                      without any prob. of unobs. files being created.
#   max.of.mins.calc : Logical. If TRUE will calculate the minimum RMSD of each
#                      line of the subsampled matrix and return the max among
#                      these values, without any prob. of unobs. files being
#                      created. Both max.rmsd.calc and max.of.mins.calc be FALSE
#                      in order for the function to be allowed to create the
#                      prob. of unobs. files.
  if (sampling == 1 && max.rmsd.calc == F) {
    results         <- ComputeClusters(my.matrix, min.rmsd, rmsd.step)
    results.x.final <- results[, 1]
    results.y.final <- results[, 2]
  } else {
    results.x         <- list()
    results.y         <- list()

    temp.max.rmsds    <- vector()
    temp.maxs.of.mins <- vector()
    temp.max.of.mins  <- vector()
    for (i in 1:sampling) {
      new.matrix <- my.matrix[seq(i, nrow(my.matrix), by=sampling),
                              seq(i, ncol(my.matrix), by=sampling) ]

      if (max.rmsd.calc == F && max.of.mins.calc == F) {
        results        <- ComputeClusters(new.matrix, min.rmsd, rmsd.step)
        results.x[[i]] <- results[, 1]
        results.y[[i]] <- results[, 2]
      }

      if (max.rmsd.calc == T) {
        temp.max.rmsds[i] <-
          DetermineFirstDiagonalMaxRmsd(new.matrix, sqrt(length(new.matrix)))
      }

      if (max.of.mins.calc == T) {
        for (j in 1:sqrt(length(new.matrix))) {
          sorted               <- sort(new.matrix[j, ])
          temp.maxs.of.mins[j] <- sorted[2]
        }
      temp.max.of.mins[i] <- max(temp.maxs.of.mins)
      }
    }

    if (max.rmsd.calc == T && max.of.mins.calc == T) {
      return(c(mean(temp.max.rmsds)  , sd(temp.max.rmsds),
               mean(temp.max.of.mins), sd(temp.max.of.mins)))
    } else if (max.rmsd.calc == T) {
      return(c(mean(temp.max.rmsds)  , sd(temp.max.rmsds)))
    } else if (max.of.mins.calc == T) {
      return(c(mean(temp.max.of.mins), sd(temp.max.of.mins)))
    }

    results.x.final <- DetermineMinLength(results.x)
    results.y.final <- vector()
    results.z.final <- vector()
    for (i in 1:length(results.x.final)) {
      tempvector <- vector()
      for (j in 1:sampling){
        tempvector[j] <- unlist(results.y[[j]][i])
      }

      results.y.final[i] <- mean(tempvector)
      results.z.final[i] <- sd(tempvector)
    }
  }

  postscript(file=paste('good_turing.', file.name, '.P_unobserved_vs_RMSD.eps',
                        sep=''))
  par(mar=c(7, 4, 4, 2))

  plot(results.x.final, results.y.final,
       main=paste('Probability of unobserved species vs RMSD\n', file.name),
       xlab='RMSD', ylab='P_unobserved', type='o', col='black')
  mtext(paste('Expected maximal RMSD upon doubling of the simulation time is',
        round(max.of.mins[sampling], 3), 'Angstrom.'), side=1, line=5)
  legend('topright', '(x,y)',
         c(paste(sampling, 'X', sep='')),
         lwd=c(2.5),
         lty=c(1))

  if (sampling > 1) {

    segments(results.x.final, results.y.final - results.z.final,
             results.x.final, results.y.final + results.z.final)
    segments(results.x.final - 0.01, results.y.final - results.z.final,
             results.x.final + 0.01, results.y.final - results.z.final)
    segments(results.x.final - 0.01, results.y.final + results.z.final,
             results.x.final + 0.01, results.y.final + results.z.final)

    write.table(cbind(round(results.x.final, 4),
                      round(results.y.final, 4),
                      round(results.z.final, 4)),
                file=paste('good_turing.',
                           file.name,
                           '.P_unobserved_vs_RMSD.dat',
                           sep=''),
                col.names=c('RMSD', 'P_un', 'Deviation'),
                row.names=F,
                sep='\t')
  } else {
    write.table(cbind(round(results.x.final, 4),
                      round(results.y.final, 4)),
                file=paste('good_turing.',
                           file.name,
                           '.P_unobserved_vs_RMSD.dat',
                           sep=''),
                col.names=c('RMSD', 'P_un'),
                row.names=F,
                sep='\t')
  }

  dev.off()
}

################################################################################
##                                                                            ##
##                         Main part of the program                           ##
##                                                                            ##
################################################################################

# This is the error message that is printed every time a bad matrix is read.
matrix.error.msg <- '
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
'

################################################################################
##                                                                            ##
##                    Declaration of the global variables                     ##
##                                                                            ##
################################################################################

# These are the samplings that will be used for the max RMSDs and max of mins
# calculations.
samplings <- c(1, 2, 3, 4, 6, 8, seq(10, 50, 4), seq(55, 100, 5),
               seq(110, 200, 10))
nofsamplings <- length(samplings)
# This is the number of iterations the ComputeClusters function will go
# through, in essence determining the number of points in the final plot.
rmsd.step.iterations <- 20
# This is the sampling factor after which the max rmsd and max of mins
# are deemed unrelialbe.
sampling.cutoff <- 100
# This is the number of iterations the non linear fitting function performs
# in order to determine the constants. Increase only if the function does
# not converge with the default value.
nls.iter <- 150
# This is the number of standard deviations the max.rmsd value should differ
# from A0.
sigma.factor <- 0.5
# Boolean variable which determines if weighted non linear fitting is performed.
# Default weights are the inverse of the variances.
weighted.fitting <- FALSE

# Find out the platform the program runs on.
os <- .Platform$OS.type

# Interactive file selection.
input.file <- file.choose()
cat('\n\n==================================================\n')
cat('1. Opening file for reading ...')

if (DetermineIfBinary(input.file) == T) {
  cat('              ERROR\n',
      '\n**************************************************\n',
      'It appears you  provided a binary  file. Aborting.\n',
      '**************************************************\n\n', sep='')
  cat(matrix.error.msg)
  stop()
} else {
  header.test <- try(scan(input.file, nmax=10, quiet=T), silent=T)
  if (length(grep('error', header.test, ignore.case=T))) {
    cat('              ERROR\n',
        '\n**************************************************\n',
        'It appears the matrix contains  a  header.  Please\n',
        'remove the lines corresponding  to  the header and\n',
        'resubmit the file.\n',
        '**************************************************\n\n', sep='')
    cat(matrix.error.msg)
    stop()
  }
}

# Get the relative/absolute path of the file.
if (os == 'windows') {
  temp.file.var <- unlist(strsplit(input.file, '[\\]'))
} else {
  temp.file.var <- unlist(strsplit(input.file, '/'))
}

# Remove the path to the file and keep only the file name.
if (length(temp.file.var) > 1) {
  file.name <- temp.file.var[length(temp.file.var)]
} else {
  file.name <- input.file
}

# Remove the file extension from the file name (if any) as well.
temp.file.var <- unlist(strsplit(file.name, '[.]'))
if (length(temp.file.var) > 1) {
  file.name <- ''
  for (i in 1:(length(temp.file.var) - 1)) {
    if (i == length(temp.file.var) - 1) {
      file.name <- paste(file.name, temp.file.var[i], sep='')
    } else {
      file.name <- paste(file.name, temp.file.var[i], '.', sep='')
    }
  }
}

cat('                 OK\n')
cat('2. Checking dimensions ...')
if (os == 'windows') {
  nofrows    <- length(scan(input.file, what='raw', sep='\n', quiet=T))
  blank.rows <- length(readLines(input.file)) - nofrows
} else {
  nofrows    <- as.integer(system(paste('wc -l < ', input.file),
                           intern=T))
  blank.rows <- nofrows - as.integer(system(paste('grep -cve \'^\\s*$\'',
                                                  input.file), intern=T))
  nofrows    <- nofrows - blank.rows
}

# Determine if the values are space or tab delimited and after that if the first
# character is the delimiter or not.
tab.dlm <- try(scan(input.file, nlines=1 + blank.rows, sep='\t', quiet=T),
               silent=T)
if (length(grep('error', tab.dlm, ignore.case=T))) {
  nofcols <- length(scan(input.file, nlines=1 + blank.rows, sep=' ', quiet=T))
  tab.dlm <- 0
} else {
  nofcols <- length(scan(input.file, nlines=1 + blank.rows, sep='\t', quiet=T))
  tab.dlm <- 1
}

if (tab.dlm == 0) {
  if (is.na(scan(input.file, nlines=1+blank.rows, sep=' ', quiet=T)[1]) == T) {
    nofcols <- nofcols - 1
  }
} else {
  if (is.na(scan(input.file, nlines=1+blank.rows, sep='\t', quiet=T)[1]) == T) {
    nofcols <- nofcols - 1
  }
}

# Dimension checks. Abort / Warn on failure.
if (nofrows != nofcols) {
  cat('                   ERROR\n\n')
  cat('Number of columns :', nofcols, '\n')
  cat('Number of rows    :', nofrows, '\n')
  cat('\n**************************************************',
      '\n      Non square matrix provided. Aborting.',
      '\n**************************************************\n\n', sep='')
  cat(matrix.error.msg)
  stop()
} else if (nofrows < 400) {
  cat('                 WARNING\n\n')
  cat('Number of columns :', nofcols, '\n')
  cat('Number of rows    :', nofrows, '\n')
  cat('\n**************************************************\n',
      'This program  is tuned  for large  scale  problems\n',
      '(with matrices  of  the  order of  thousands).  It\n',
      'can  not  presently  deal with  matrices less than\n',
      '400x400.\n',
      '**************************************************\n\n', sep='')
  stop()
}

if (nofrows >= 400) {
  cat('                      OK\n')
}
cat('3. Loading the matrix ...')
rmsd.matrix <- matrix(scan(input.file, n=nofrows*nofcols, quiet=T),
                      nofrows, nofcols, byrow=T)
cat(sprintf("            %5d x %5d\n", nofrows, nofcols))

cat('4. Sanity checks :\n')
# Perform the matrix checks.
# Is the matrix full of 0s ?
matrix.length <- nofrows * nofcols
cat('   Non-zero matrix ? ')
for (i in 2:matrix.length) {
  if (rmsd.matrix[i] != 0) {
    cat('                          yes\n')
    break
  } else if (i == matrix.length) {
    cat('                           NO\n',
        '\n**************************************************\n',
        'The matrix appears to  be full of zeros. Aborting.\n',
        '**************************************************\n\n', sep='')
    cat(matrix.error.msg)
    stop()
  }
}

# Calculate the minimun and the maximum RMSD of the matrix.
min.rmsd <- min(rmsd.matrix[upper.tri(rmsd.matrix)])
max.rmsd <- max(rmsd.matrix[upper.tri(rmsd.matrix)])

# Does the matrix contain negative numbers ?
cat('   Negative values present ?')
if (min.rmsd < 0) {
  cat('                   YES\n',
      '\n**************************************************\n',
      'Negative values detected  in the matrix. Aborting.\n',
      '**************************************************\n\n', sep='')
  cat(matrix.error.msg)
  stop()
} else {
  cat('                  none\n')
}

# Does the matrix contain abnormally large numbers ?
cat('   Sane RMSDs ?')
if (max.rmsd > 50) {
  cat('                                 NO\n',
      '\n**************************************************\n',
      '             RMSDs higher than 50A ?\n',
      '      Are you sure this is an RMSD matrix ?\n',
      '      Continuing and hoping for the best ...\n',
      '**************************************************\n\n', sep='')
} else {
  cat('                                yes\n')
}

# Does the diagonal of the matrix contain numbers other than zeros  ?
cat('   Origin of matrix at expected position ?')
for (i in 1:nofrows) {
  if (rmsd.matrix[i, i] > 0.01) {
    cat('      NO\n',
        '\n**************************************************\n',
        'Values greater than 0.01 detected in the diagonal.\n',
        'Aborting.\n',
        '**************************************************\n\n', sep='')
    cat(matrix.error.msg)
    stop()
  }
}
cat('     yes\n')

# Is the matrix symmetrical ?
cat('   Average deviation from symmetry ')
dev        <- 0
count.dots <- 1
for (i in 1:nofrows) {
  if (i == round( (nofrows / 10) * count.dots)) {
    cat('.')
    count.dots <- count.dots + 1
  }
  dev <- dev + sum(abs(rmsd.matrix[i, (i:nofrows)] -
                       rmsd.matrix[(i:nofrows), i]))
}

dev <- dev / (nofrows * (nofrows - 1) / 2 + nofrows)

if (dev >= 0.01) {
  cat(' FAIL\n',
      '\n**************************************************\n',
      '  Are you sure this is a symmetric RMSD matrix ?\n',
      '      Continuing and hoping for the best ...\n',
      '**************************************************\n\n', sep='')
} else {
  cat(' pass\n')
}

# Are there outlying values among the symmetrical pairs ?
cat('   Symmetry outliers ? ')
count.dots        <- 1
count.warnings    <- 0
errors            <- vector()
symmetry.errors.i <- vector()
symmetry.errors.j <- vector()
for (i in 1:nofrows) {
  if (i == round( (nofrows / 22) * count.dots)) {
    cat('.')
    count.dots <- count.dots + 1
  }
  errors <- which(rmsd.matrix[i, (i:nofrows)] -
                  rmsd.matrix[(i:nofrows), i] != 0)
  if (length(errors) > 0) {
    for (j in 1:length(errors)) {
      dif <- abs(rmsd.matrix[i, errors[j] + i - 1] -
                 rmsd.matrix[errors[j] + i - 1, i]) /
             (rmsd.matrix[i, errors[j] + i - 1] +
              rmsd.matrix[errors[j] + i - 1, i])
      if (!is.na(dif) && dif >= 0.2) {
        count.warnings <- count.warnings + 1
        symmetry.errors.i[count.warnings] <- i
        symmetry.errors.j[count.warnings] <- errors[j] + i - 1
      }
    }
  }
}

if (count.warnings == 0) {
  cat(' none\n')
} else {
  cat(sprintf(" %4d\n", count.warnings),
      '\n**************************************************\n',
      '    Symmetry - related  values  at  positions\n\n', sep='')
  for (i in 1:count.warnings) {
    cat('             ',
        sprintf("%5d", symmetry.errors.i[i]), '-',
        sprintf("%5d", symmetry.errors.j[i]), ', ',
        sprintf("%5d", symmetry.errors.j[i]), '-',
        sprintf("%5d", symmetry.errors.i[i]), '\n', sep='')
  }
  cat('\n             differ by more than 20%.\n',
      '**************************************************\n\n', sep='')
  rm(dif)
}

# All checks successful. Proceed normally.
rm(max.rmsd, matrix.length, dev, count.warnings,
   symmetry.errors.i, symmetry.errors.j)

# Determine the max of mins for the original matrix, ie for a subsampling factor
# of 1. Review the comments of CreateSubmatrix for an explanation of how this
# calculation is performed. The reasoning behind this analysis is that the max
# of mins that corresponds to each sampling should present a fairly accurate
# prediction of what differences to the already seen RMSDs one should expect to
# observe should the simulation time be doubled.
maxs.of.mins     <- vector()
max.of.mins      <- vector()
max.of.mins.devs <- vector()

for (i in 1:nofrows) {
  sorted          <- sort(rmsd.matrix[i,])
  maxs.of.mins[i] <- sorted[2]
}

max.of.mins[1]      <- max(maxs.of.mins)
max.of.mins.devs[1] <- NA

# Determine the RMSD step variable. For more info review the comments of the
# function that is called.
rmsd.step <- DetermineRmsdStep(rmsd.matrix, min.rmsd, rmsd.step.iterations)

max.rmsds     <- vector()
max.rmsd.devs <- vector()

cat('5. Sampling determination ')
count.dots <- 1
for (i in 1:nofsamplings) {
  results <- vector()
  results <- CreateSubmatrix(rmsd.matrix, samplings[i], T, T)

  max.rmsds[i]     <- results[1]
  max.rmsd.devs[i] <- results[2]
  if (i > 1) {
    max.of.mins[i]      <- results[3]
    max.of.mins.devs[i] <- results[4]
  }

  if (i == round( (nofsamplings / 19) * count.dots)) {
    cat('.')
    count.dots <- count.dots + 1
  }
}
cat('.  OK\n')
cat('6. Calculating probability curve ...')

max.rmsd.variances    <- max.rmsd.devs^2
max.of.mins.variances <- max.of.mins.devs^2

max.rmsd.min.var <- min(max.rmsd.variances[which(max.rmsd.variances > 0)])
max.mins.min.var <- min(max.of.mins.variances[which(max.of.mins.variances > 0)])

for (i in 1:nofsamplings) {
  if (is.na(max.rmsd.variances[i]) == T || max.rmsd.variances[i] == 0) {
    max.rmsd.variances[i]    <- max.rmsd.min.var
  }
  if (is.na(max.of.mins.variances[i]) == T || max.of.mins.variances[i] == 0) {
    max.of.mins.variances[i] <- max.mins.min.var
  }
}

# The non linear curve fitting that follows is performed on the max RMSDs and
# max of mins data in order to determine at which sampling the probability of
# unobserved vs RMSD analysis should be performed, or even if it should be
# performed at all. The equation itself looks like :
# f(x) = (x + A2) * ( (1 + |(x + A2) / A0| ^ A1) ^ (-1 / A1))
# and has applications in electronic circuit regulation, where it is used to
# limit the upper voltage that can pass through a diode, in order to prevent the
# circuit from being damaged.
if (weighted.fitting == T) {
  custom.regressionA <- nlsLM(max.rmsds ~
                              I(samplings + a2) *
                              I((1 + abs((samplings + a2) / a0) ^ a1)^(-1.0/a1)),
                              start=list(a0=1, a1=1, a2=1),
                              control=nls.lm.control(maxiter=nls.iter),
                              weights=I(1 / max.rmsd.variances))

  custom.regressionB <- nlsLM(max.of.mins ~
                              I(samplings + a2) *
                              I((1 + abs((samplings + a2) / a0) ^ a1)^(-1.0/a1)),
                              start=list(a0=1, a1=1, a2=1),
                              control=nls.lm.control(maxiter=nls.iter),
                              weights=I(1 / max.of.mins.variances))
} else {
  custom.regressionA <- nlsLM(max.rmsds ~
                              I(samplings + a2) *
                              I((1 + abs((samplings + a2) / a0) ^ a1)^(-1.0/a1)),
                              start=list(a0=1, a1=1, a2=1),
                              control=nls.lm.control(maxiter=nls.iter))

  custom.regressionB <- nlsLM(max.of.mins ~
                              I(samplings + a2) *
                              I((1 + abs((samplings + a2) / a0) ^ a1)^(-1.0/a1)),
                              start=list(a0=1, a1=1, a2=1),
                              control=nls.lm.control(maxiter=nls.iter))
}

# The a0 variable is the RMSD at which the max RMSDs reach a plateau.
a0 <- summary(custom.regressionA)$coefficients[1]
a1 <- summary(custom.regressionA)$coefficients[2]
a2 <- summary(custom.regressionA)$coefficients[3]

# Same as above but for the max of mins data.
b0 <- summary(custom.regressionB)$coefficients[1]
b1 <- summary(custom.regressionB)$coefficients[2]
b2 <- summary(custom.regressionB)$coefficients[3]

# This loop is used to determine at which point, if any, the max RMSD is greater
# than a0, ie at which point the max RMSDs level off. This point is then used to
# determine which is the smallest sampling that can be used whose max RMSD is
# greater than a0.
max.rmsd.devs[1] <- 0  # This is so that the else if check below is possible
outcome <- 0
for (i in 1:nofsamplings) {
  # If the sampling exceeds sampling.cutoff, then the data are deemed as
  # unsuitable for prob. of unobserved species vs RMSD analysis due to the 
  # very small size of the resulting matrices [for example, if the value of
  # the sampling.cutoff is 100, the matrices would be (1/100)th of the original].
  if (samplings[i] >= sampling.cutoff) {
    outcome <- 1
    break
  } else if ( (max.rmsds[i] + (sigma.factor * max.rmsd.devs[i])) >= a0) {
    top.sampling = samplings[i]
    for (j in 1:top.sampling) {
      if ( ! j %in% samplings) {
        results <- CreateSubmatrix(rmsd.matrix, j, T, T)

        temp.max.rmsds        <- results[1]
        temp.max.rmsd.devs    <- results[2]
        temp.max.of.mins      <- results[3]
        temp.max.of.mins.devs <- results[4]

        old.samplings <- samplings
        samplings     <- sort(c(samplings, j))

        sampling.diffs <- samplings %in% old.samplings
        insert.point   <- max(which(sampling.diffs == FALSE)) - 1

        max.rmsds     <- append(max.rmsds, temp.max.rmsds, after=insert.point)
        max.rmsd.devs <- append(max.rmsd.devs, temp.max.rmsd.devs,
                                after=insert.point)
        max.of.mins      <- append(max.of.mins, temp.max.of.mins,
                                   after=insert.point)
        max.of.mins.devs <- append(max.of.mins.devs, temp.max.of.mins.devs,
                                   after=insert.point)

        if ( (temp.max.rmsds + (sigma.factor * temp.max.rmsd.devs)) >= a0) {
          i <- which(samplings == j)
          break
        }
      } else if (j == top.sampling) {
        i <- which(samplings == j)
      }
    }

    if (length(samplings) > nofsamplings) {
      rm(old.samplings, sampling.diffs, insert.point, temp.max.rmsds, j,
         temp.max.rmsd.devs, temp.max.of.mins, temp.max.of.mins.devs,
         top.sampling)
    }

    CreateSubmatrix(rmsd.matrix, samplings[i], F, F)
    break
  }
}
max.rmsd.devs[1] <- NA  # Reversing the assignment of zero.

cat('            OK\n')
cat('7. Writing files ...')

dir.create(paste('good_turing.', file.name, sep=''))
setwd(paste('good_turing.', file.name, sep=''))

write.table(round(cbind(samplings, max.rmsds, max.rmsd.devs), 4),
            file=paste('good_turing.', file.name, '.max_rmsds.dat', sep=''),
            row.names=F,
            col.names=F)

write.table(round(cbind(samplings, max.of.mins, max.of.mins.devs), 4),
            file=paste('good_turing.', file.name, '.max_of_mins.dat', sep=''),
            row.names=F,
            col.names=F)

postscript(file=paste('good_turing.', file.name, '.max_rmsds.eps', sep=''))
plot(samplings,
     max.rmsds,
     ylim=c(min(max.rmsds) - max(max.rmsd.devs[2:length(max.rmsd.devs)]),
            max(max.rmsds) + max(max.rmsd.devs[2:length(max.rmsd.devs)])),
     main='Max RMSDs vs sampling',
     xlab='Sampling',
     ylab='RMSD')
segments(samplings, max.rmsds - max.rmsd.devs,
         samplings, max.rmsds + max.rmsd.devs)
segments(samplings - 0.5, max.rmsds - max.rmsd.devs,
         samplings + 0.5, max.rmsds - max.rmsd.devs)
segments(samplings - 0.5, max.rmsds + max.rmsd.devs,
         samplings + 0.5, max.rmsds + max.rmsd.devs)
lines(c(1:max(samplings)),
      predict(custom.regressionA, list(samplings=c(1:max(samplings)))),
      lty=1, col='blue')
dev.off()

postscript(file=paste('good_turing.', file.name, '.max_of_mins.eps', sep=''))
plot(samplings,
     max.of.mins,
     ylim=c(min(max.of.mins)-max(max.of.mins.devs[2:length(max.of.mins.devs)]),
            max(max.of.mins)+max(max.of.mins.devs[2:length(max.of.mins.devs)])),
     main='Max of minimum RMSDs vs sampling',
     xlab='Sampling',
     ylab='RMSD')
segments(samplings, max.of.mins - max.of.mins.devs,
         samplings, max.of.mins + max.of.mins.devs)
segments(samplings - 0.5, max.of.mins - max.of.mins.devs,
         samplings + 0.5, max.of.mins - max.of.mins.devs)
segments(samplings - 0.5, max.of.mins + max.of.mins.devs,
         samplings + 0.5, max.of.mins + max.of.mins.devs)
lines(c(1:max(samplings)),
      predict(custom.regressionB, list(samplings=c(1:max(samplings)))),
      lty=1, col='blue')
dev.off()

setwd('..')
tar(paste('good_turing.', file.name, '.tar', sep=''),
    paste('good_turing.', file.name, sep=''),
    tar='internal')
Sys.sleep(2)
unlink(paste('good_turing.', file.name, sep=''), recursive=T)

cat('                            OK\n')

if (outcome == 1) {
  cat('\nSummary :\n\n',
      'The maximal RMSDs  between the observed trajectory\n',
      'structures  have  not converged. This implies that\n',
      'the  length  of  the  given  trajectory  does  not\n',
      'suffice for meaningfully quantifying  convergence.\n',
      'The only comment  that  can safely be made is that\n',
      'upon  doubling  the  simulation  time  you  should\n',
      'expect  to  observe  structures  that  differ from\n',
      'those already observed by more than  approximately\n',
      paste(sprintf("%.1f", max(max.of.mins)), 'Angstrom.'), sep='')

  if ( (max(max.rmsds) * 1.1) <= a0) {
    cat('\n\nAdditional note : based  solely on the  given RMSD\n',
        'matrix, it  would  appear  that  the corresponding\n',
        'trajectory is  actively  evolving  and  is nowhere\n',
        'near convergence.', sep='')
  }
} else {
  cat('\nThe graph and corresponding data from this run are\n',
      'contained in the files with the extensions\n\n',
      '    good_Turing.*.P_unobserved_vs_RMSD.eps   and\n',
      '    good_Turing.*.P_unobserved_vs_RMSD.dat\n\n',
      '==================================================\n',
      'Summary :\n\n',
      'The maximal RMSDs of the trajectory converged with\n',
      'a sub-sampling factor of ', sprintf("%2d. ", samplings[i]),
      'The analysis suggests\n',
      'that  the  most  different  structure  you  should\n',
      'expect to observe if  you  double  the  simulation\n',
      'time    will    differ    by    no     more   than\n', sep='')

  if (i > 1 && max.of.mins.devs[i] >= 0.05 * max.of.mins[i]) {
    dev.to.char <- unlist(strsplit(as.character(max.of.mins.devs[i]), split=''))
    for (j in 3:length(dev.to.char)) {
      if (as.numeric(dev.to.char[j]) != 0) {
        break
      }
    }

    cat('approximately ', sprintf(paste("%.1f +- %.", j-2, "f", sep=''),
        max.of.mins[i], max.of.mins.devs[i]),
        ' Angstrom  (RMSD)  from\n',
        'those already observed.', sep='')
    rm(j)
  } else {
    cat('approximately ', sprintf("%.1f", max.of.mins[i]),
        ' Angstrom   (RMSD)   from   those\n',
        'already observed.', sep='')
  }
}

cat('\n==================================================\n\n\n')

rm(ComputeClusters, DetermineMinLength, DetermineRmsdStep, CreateSubmatrix,
   DetermineFirstDiagonalMaxRmsd, DetermineIfBinary, i, input.file, os, nofrows,
   results, min.rmsd, nofcols, rmsd.matrix, file.name, rmsd.step, maxs.of.mins,
   max.of.mins, sorted, temp.file.var, a0, a1, a2, b0, b1, b2, blank.rows,
   custom.regressionA, custom.regressionB, errors, header.test, sigma.factor,
   matrix.error.msg, max.mins.min.var, max.of.mins.devs, max.of.mins.variances,
   max.rmsd.devs, max.rmsd.min.var, max.rmsd.variances, max.rmsds, outcome,
   tab.dlm, samplings, nls.iter, rmsd.step.iterations, sampling.cutoff, 
   nofsamplings, count.dots, weighted.fitting)
