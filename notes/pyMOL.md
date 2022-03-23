## Installing pyMOL open source on Ubuntu

At some point I will write this up a little nicer.

```
git clone git@github.com:schrodinger/pymol-open-source.git

git clone git@github.com:rcsb/mmtf-cpp.git
cd mmtf-cpp/
mkdir build
cd build/
sudo apt-get install ninja-build
cmake -G Ninja ..
sudo ninja install
cd ../../

cd pymol-open-source/
sudo apt-get install python3 
sudo apt-get install pip3
python3 -m pip install --upgrade pip
pip install biopython Pillow numpy

sudo apt-get install libpng-dev
wget  https://anaconda.org/schrodinger/collada2gltf/2.1.4/download/linux-64/collada2gltf-2.1.4-h6bb024c_0.tar.bz2
tar xfv collada2gltf-2.1.4-h6bb024c_0.tar.bz2
cp bin/collada2gltf ~/bin/

pip3 install pmw
sudo apt-get install python3-tk

python3 setup.py --glut install --prefix=~/bin
~/bin/bin/pymol
```
