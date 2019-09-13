## Installing ETE using the Anaconda/Miniconda environment
### Install Minconda  (you can ignore this step if you already have Anaconda/Miniconda)
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;
### Install ETE
conda install -c etetoolkit ete3 ete_toolchain
### Check installation
ete3 build check

## fishbase (optionnal)
pip3 install -U fishbase
