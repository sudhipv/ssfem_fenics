#!/usr/bin/env bash

## Make sure you have "foo.msh" file in the "/data": if not follow this
## have foo.geo file in "/data" folder and uncomment following two lines
#NP=4
#gmsh -2 foo.geo -part $NP -o foo.msh

## Using Matlab code to decompose patitioned mesh into points & elements
##cd preprocessors
matlabTerminal -r Kuf_IntrusiveND
