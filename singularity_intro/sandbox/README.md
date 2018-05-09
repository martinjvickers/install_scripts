## Create a CentOS sandbox image

```
sudo singularity build --sandbox centos_sandbox centos_sandbox.def
```

## Make changes to the sandbox

```
sudo singularity shell -w centos_sandbox
```

NOTE: the need for `-w`, without it changes will not persist.

## Build the sandbox into a production image

This would create a read-only production squashfs compressed
image like the samtools examples. Simply add `-w` to make
writeable with (~20% space) overhead.

```
sudo singularity build production.img centos_sandbox
```
