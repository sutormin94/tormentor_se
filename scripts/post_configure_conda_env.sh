#!/usr/bin/env bash

$CONDA_EXE init $(basename $SHELL)
reset
eval "$($CONDA_EXE shell.bash hook)"
conda activate tormentor
echo $PATH
conda env config vars set PATH=$PATH:$PWD/bin/
conda activate tormentor

ln -s bin/circuclust vnom/dependencies/circuclust
ln -s bin/mars vnom/dependencies/mars
ln -s bin/usearch vnom/dependencies/usearch
