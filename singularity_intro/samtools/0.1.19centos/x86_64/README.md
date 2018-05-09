Making a basic samtools image using the `yum` package manager in CentOS.

Includes correct directory layout for use in `/nbi/software/testing`.

```
sudo singularity build samtools-0.1.19centos.img ../src/samtools-0.1.19-centos.def
[mvickers@NBI-HPC x86_64]$ cd bin/
[mvickers@NBI-HPC bin]$ ln -s singularity.exec samtools
[mvickers@NBI-HPC bin]$ cd ../
[mvickers@NBI-HPC x86_64]$ bin/samtools 
WARNING: Could not chdir to home: /usr/users/JIC_c1/mvickers

Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-44428cd

Usage:   samtools <command> [options]

Command: view        SAM<->BAM conversion
         sort        sort alignment file
         mpileup     multi-way pileup
         depth       compute the depth
         faidx       index/extract FASTA
         tview       text alignment viewer
         index       index alignment
         idxstats    BAM index stats (r595 or later)
         fixmate     fix mate information
         flagstat    simple stats
         calmd       recalculate MD/NM tags and '=' bases
         merge       merge sorted alignments
         rmdup       remove PCR duplicates
         reheader    replace BAM header
         cat         concatenate BAMs
         bedcov      read depth per BED region
         targetcut   cut fosmid regions (for fosmid pool only)
         phase       phase heterozygotes
         bamshuf     shuffle and group alignments by name
```
