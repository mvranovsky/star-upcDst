#!/bin/bash

if [ "$#" -eq 0 ]; then
	du -hs /gpfs01/star/pwg/mvranovsk/Run17_P20ic/* | sort -hr
    exit 1
else
	for (( i=1; i<=$#; i++ )); do
		rm -r "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/${!i}"
	done
fi

