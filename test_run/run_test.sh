#!/bin/bash
set -e

rm -rf run Results_MothurMock_* taxonomy
rm -f bashMe.sh
find ../pipelineScripts -type d -exec chmod +rwx {} \;
find ../pipelineScripts -type f -exec chmod +rwx {} \;
mkdir -p taxonomy/curated
gunzip <../taxonomy/curated/16SMicrobial_curated.gzip >./taxonomy/curated/16SMicrobial_curated.txt
cp ../pipelineScripts/shell/bashMe.sh ./
ln -nsf ../pipelineScripts
ln -nsf ../db

bash bashMe.sh

echo -e "ANCHOR is over. Well done!"
rm -f ../taxonomy/curated/16SMicrobial_curated.txt
