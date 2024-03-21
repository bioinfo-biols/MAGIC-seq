#!/bin/bash
sample=$1

sampath=`dirname ${sample}`
ifile1=`basename ${sample}`
array2=(${ifile1//_R1/ })
sample=${array2[0]}
