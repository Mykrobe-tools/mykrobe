#!/usr/bin/env bash
set -vex
wget -O mykrobe-data.tar.gz https://ndownloader.figshare.com/files/20996829
tar -xvf mykrobe-data.tar.gz
rm mykrobe-data.tar.gz
rm -fr src/mykrobe/data
mv mykrobe-data src/mykrobe/data
