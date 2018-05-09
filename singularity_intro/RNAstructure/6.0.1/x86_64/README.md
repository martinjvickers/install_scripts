Example courtesy of Hugh Wollfenden (JIC) slightly modified.

## Build the image from recipe

```
sudo singularity build rnastructure-6.0.1.img ../src/rnastructure-6.0.1.def
```

The key interest is that the DATAPATH needs to be set. This can be two ways;

* By preappending `SINGULARITYENV_*` before `singularity exec`

e.g. 

`SINGULARITYENV_DATAPATH=/opt/software/RNAstructure/data_tables`

* Using the `%environment` section in the recipe

e.g.

```
%environment
  DATAPATH=/opt/software/RNAstructure/data_tables
```

## Sym-links

The file `rnastructure_exes` contains a list of all the RNAstructure
executables that are needed to run the software which are located in
`/usr/local/bin` directory in the container. 

The following command symlinks the `singularity.exec` file to the 
different exes within the container.

```
awk '{print "ln -s singularity.exec "$1}' rnastructure_exes
```

To test, this;

```
[mvickers@NBI-HPC x86_64]$ cd bin/
[mvickers@NBI-HPC bin]$ ./Fold
Incorrect number of required parameters given. (Found 0 but expected 2.)
USAGE: Fold <seq file> <ct file> [options]
Use any of the following options to get a help message: -h --help
```
