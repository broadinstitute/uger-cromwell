# uger-cromwell
Cromwell template for running on UGER with Singularity

This repo is meant to get you started with running cromwell on-prem on the UGER cluster as well as running in the cloud via PAPIv2 using the same WDL file.
It is meant as a supporting document to the WDL tutorial here: https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-tutorials
It is not meant to explain how WDL/Cromwell/UGER/Google work.

## Getting Started
Clone this repo to NFS shared storage on a Broad host (gold/silver/plat/login) such as your home directory or /broad/hptmp.
Assuming that Broad host is a UGER sumbit host please run:

```
$ use UGER
$ ish -l h_vmem=4G
$ cd /your/nfs/path/uger-cromwell

ex.

cd ~/uger-cromwell #if you cloned to your home directory
```

This will give you an instance on a UGER compute node where you can run cromwell and move you into the directory for this tutorial. To get started with cromwell you will first need to download the latest version from:
https://github.com/broadinstitute/cromwell/releases (the current tested version with this documentation is cromwell-52).

This can be accomplished with a command such as:

```
wget https://github.com/broadinstitute/cromwell/releases/download/52/cromwell-52.jar
```

You will also need a copy of womtool for the same version of cromwell.

```
wget https://github.com/broadinstitute/cromwell/releases/download/52/womtool-52.jar
```

You can now setup your environment to work with cromwell, uger, and singularity. Please run the following from within the uger-cromwell repo:

```
source 00-setup/setup.sh <your-broad-username>
```

This will create a temp NFS storage space for you and direct your singularity images to be stored there to not fill your home directory.

It is also advised that you prebuild any container which you will work with as it can take a bit of time.
In a new terminal:
```
$ ssh login
$ use UGER
$ ish -l h_vmem=4G -pe smp 4 -binding linear:4

$ source ~/uger-cromwell/00-setup/setup.sh <your-broad-username>
$ singularity shell docker://broadinstitute/gatk@sha256:8051adab0ff725e7e9c2af5997680346f3c3799b2df3785dd51d4abdd3da747b
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
Copying blob 423ae2b273f4 skipped: already exists
Copying blob de83a2304fa1 skipped: already exists
Copying blob f9a83bce3af0 skipped: already exists
Copying blob b6b53be908de skipped: already exists
Copying blob a69d35efa09a skipped: already exists
Copying blob 95d935bf5bc5 skipped: already exists
Copying config 1fe6830994 done
Writing manifest to image destination
Storing signatures
2020/08/07 20:36:36  info unpack layer: sha256:423ae2b273f4c17ceee9e8482fa8d071d90c7d052ae208e1fe4963fceb3d6954
2020/08/07 20:36:37  info unpack layer: sha256:de83a2304fa1f7c4a13708a0d15b9704f5945c2be5cbb2b3ed9b2ccb718d0b3d
2020/08/07 20:36:37  info unpack layer: sha256:f9a83bce3af0648efaa60b9bb28225b09136d2d35d0bed25ac764297076dec1b
2020/08/07 20:36:37  info unpack layer: sha256:b6b53be908de2c0c78070fff0a9f04835211b3156c4e73785747af365e71a0d7
2020/08/07 20:36:37  info unpack layer: sha256:a69d35efa09a62480558d949be2ca9e249dc35704ed4749f5da647e34e52c791
2020/08/07 20:37:05  info unpack layer: sha256:95d935bf5bc50aed34edabc910e0bd5e4dc0ed1e54537d9c7d9cf6ebd2ff51b3
INFO:    Creating SIF file...
Singularity> exit
```

You are now ready to run your first job back in your original terminal window.


## Examples

1. [Write your first WDL script running GATK HaplotypeCaller](01-helloHaplotypeCaller/README.md)
1. [Write a simple multi-step workflow](02-simpleVariantSelection/README.md)
1. [Run a sample variant discovery mini-pipeline](03-simpleVariantDiscovery/README.md)
1. [Use scatter-gather to joint call genotypes](04-jointCallingGenotypes/README.md)
