#!/bin/bash

rm .DS_Store
rm -rf ..Rcheck
rm -rf .Rproj.user
rm R/RcppExports.R
rm .Rhistory

cd src
rm *.o
rm *.o.tmp
rm *.so
rm RcppExports.cpp
rm .fuse*
rm .nfs*

