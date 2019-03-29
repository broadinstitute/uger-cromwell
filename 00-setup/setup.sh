#!/bin/bash
reuse -q UGER
reuse -q Java-1.8
SINGULARITY_CACHEDIR=/broad/hptmp/$1
mkdir -p "${SINGULARITY_CACHEDIR}"
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR}"
