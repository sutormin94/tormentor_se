#!/usr/bin/env bash

mkdir -p bin/
cd bin/
wget https://github.com/rcedgar/circuclust/releases/download/v1.0/circuclust_linux64

mv circuclust_linux64 circuclust
chmod +x circuclust

wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz
mv usearch11.0.667_i86linux32 usearch

chmod +x usearch

git clone https://github.com/lorrainea/MARS
cd MARS/

./pre-install.sh
make -f Makefile

cp mars ../
cd ..
rm -rf MARS