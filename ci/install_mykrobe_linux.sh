#!/usr/bin/env bash
set -vexu

MYKROBE_ROOT_DIR=$1


apt update
apt install -y software-properties-common
apt-add-repository universe
apt update

apt install -y \
  build-essential \
  git \
  mongodb \
  python3-pip \
  python3-setuptools \
  python-is-python3 \
  wget


python -m pip install -U pip

cd $MYKROBE_ROOT_DIR
rm -rf mccortex
git clone --recursive -b geno_kmer_count https://github.com/Mykrobe-tools/mccortex mccortex
cd mccortex
make
cd ..
mkdir -p /data/db
python -m pip install -r requirements.txt
python -m pip install .
mongod --logpath mongo.log --quiet &>/dev/null &
sleep 3s
pytest --cov-report term-missing --cov=mykrobe
mongod --shutdown
rm mongo.log
# For whatever reason, mccortex is not getting put in the install location.
# Do it manually.
myk_dir=$(pip3 show mykrobe | awk '/^Location/ {print $NF}')
echo $myk_dir
cp mccortex/bin/mccortex31 $myk_dir/mykrobe/cortex/mccortex31
cd ..
mykrobe panels update_metadata --debug
mykrobe panels update_species --debug all

