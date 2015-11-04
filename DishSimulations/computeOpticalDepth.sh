#!/bin/bash
ls delta_T* | grep -o "_nf0.[0-9]\{6\}" | grep -o "0.[0-9]\{6\}" > nfList.txt
ls delta_T* | grep -o "_z[0-9]\{3\}.[0-9]\{2\}" | grep -o "[0-9]\{3\}.[0-9]\{2\}" > zList.txt
python computeOpticalDepth.py