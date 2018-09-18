## Intro to Singularity

## Singularity quick notes

To run this on the NBI HPC, log into the HPC (`ssh slurm`).

Log into the software node;

```
ssh software7
```

Then source git and clone this repository in your home directory;

```
source git-1.8.1.2
git clone https://github.com/martinjvickers/install_scripts
cd install_scripts/singularity_intro
```



### Create a sandbox image

```
sudo singularity build --sandbox centos_sandbox centos_sandbox.def
```

### Build squashfs images

From a recipe;

```
sudo singularity build centos.img centos.def
```

From a sandbox;

```
sudo singularity build centos_sandbox.img centos_sandbox
```

### Get a shell on the image

A writeable one;

```
sudo singularity shell -w centos.img
```

or on a sandbox

```
sudo singularity shell -w centos_sandbox
```

