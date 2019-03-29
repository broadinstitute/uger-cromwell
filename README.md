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
$ ish -l h_vmem=4G =l os=RedHat7
$ cd /your/nfs/path/uger-cromwell

ex.

cd ~/uger-cromwell #if you cloned to your home directory
```

This will give you an instance on a UGER compute node where you can run cromwell and move you into the directory for this tutorial. To get started with cromwell you will first need to download the latest version from:
https://github.com/broadinstitute/cromwell/releases (the current tested version with this documentation is cromwell-36).

This can be accomplished with a command such as:

```
wget https://github.com/broadinstitute/cromwell/releases/download/36/cromwell-36.jar
```

You will also need a copy of womtool for the same version of cromwell.

```
wget https://github.com/broadinstitute/cromwell/releases/download/36/womtool-36.jar
```

You can now setup your environment to work with cromwell, uger, and singularity. Please run the following from within the uger-cromwell repo:

```
source 00-setup/setup.sh <your-broad-username>
```

This will create a temp NFS storage space for you and direct your singularity images to be stored there to not fill your home directory.


## Examples

1. [Write your first WDL script running GATK HaplotypeCaller](01-helloHaplotypeCaller/README.md)
1. [Write a simple multi-step workflow](02-simpleVariantSelection/README.md)
1. [Run a sample variant discovery mini-pipeline](03-simpleVariantDiscovery/README.md)
1. [Use scatter-gather to joint call genotypes](04-jointCallingGenotypes/README.md)
