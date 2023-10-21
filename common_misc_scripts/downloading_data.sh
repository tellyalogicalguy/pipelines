#!/bin/bash
set -e
set -u
set -o pipefail

: << 'README'
Download data template:
Manages data storage and symlinking.
Downloads data to OUT_DIR and creates symlinks for the data in SYMLINK_DIR.
OUT_SUBDIRS is the folder structure after PROJECTPARENT_DIR.
So, whatever folder structure is defined for OUT_DIR_REL, will be created in ~/binf_analyses and ~/scratch/binf_analyses_data/ as the SYMLINK_DIR and OUT_DIR respectively.
README

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)"

#output directory
OUT_DIR_REL="${CUR_DIR}" #make the desired directory structure here
if [ ! -d $OUT_DIR_REL ]; then
	mkdir -p "$OUT_DIR_REL"
fi
OUT_DIR_ABS="$(cd "$OUT_DIR_REL" 2>&1 >/dev/null && pwd)/"
OUT_SUBDIRS="$(sed "s:$PROJECTPARENT_DIR::" <(echo $OUT_DIR_ABS))"
OUT_DIR="${SCRATCHPARENT_DIR}${OUT_SUBDIRS}"
SYMLINK_DIR="${PROJECTPARENT_DIR}${OUT_SUBDIRS}"
if [ ! -d $OUT_DIR ]; then
	mkdir -p ${OUT_DIR}
fi

#echo "OUT_DIR= "$OUT_DIR
#echo "SYMLINK_DIR= "$SYMLINK_DIR 

LINK=""
wget -P ${OUT_DIR} ${LINK} 2>> wget.log


#remove all the lines in the wget.log that are download progress indicators
awk '$1 !~ /^[/s]*[0-9]*K/{print $0}' wget.log > tmp_wget.log
cat tmp_wget.log > wget.log
rm tmp_wget.log

#creating symlinks in the home folder
if [ ! -d $SYMLINK_DIR ]; then
    mkdir -p ${SYMLINK_DIR}
fi
ln -s ${OUT_DIR}* "$SYMLINK_DIR"


echo "End of file"
