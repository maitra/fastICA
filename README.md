# fastICA

 C implementation of the basic fastICA algorithm

This is an implementation of the basic fastICA algorithm, originally written several years ago for some reason that I do not recall. (It may be that I did not find C code to implement this, which may be further validation of my poor online search abilities, but I am not sure).

Compile using `make run_fastICA` at the prompt.

Run using `./run_fastICA <parameters>` at the prompt.

Get help using `./run_fastICA -h`

````
NAME
        run_fastICA - Run the basic fastICA algorithm.
SYNOPSIS
        run_fastICA -n <int> -p <int> -# <int> -D <dir> -i <file> -v -h
OPTIONS
        -n <int>   number of observations
        -p <int>   number of dimensions
        -# <int>   desired number of independent components
        -D <dir>   working directory (default: OUTPUT)
        -i <file>  file containing the dataset
        -v verbose output, default is no verbosity
````
