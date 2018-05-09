## Singularity scripts

These are software installs on the NBI HPC.

## Singularity quick notes

### Create a disk image

sudo singularity create -s 4096 rdiff-20171027.img

### Bootstrap the image

sudo singularity bootstrap rdiff-20171027.img ../src/rdiff-20171027.def

### Get a shell on the image

sudo singularity shell -w rdiff-20171027.img


