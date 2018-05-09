## Intro to Singularity

## Singularity quick notes

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

