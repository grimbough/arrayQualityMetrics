#!/bin/bash
set -e
set -o errexit

SVNSRC=$HOME/madman/Rpacks
TEMP=$HOME/tmp

svn export --force $SVNSRC/arrayQualityMetrics $TEMP/arrayQualityMetrics
R CMD INSTALL $TEMP/arrayQualityMetrics
echo "Stangle('$TEMP/arrayQualityMetrics/inst/doc/arrayQualityMetrics.Rnw'); source('arrayQualityMetrics.R')" | R --no-save

