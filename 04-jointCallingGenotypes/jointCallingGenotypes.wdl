workflow jointCallingGenotypes {

  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict  
  Int memory_gb

  scatter (sample in inputSamples) {
    call HaplotypeCallerERC { 
      input: bamFile=sample[1], 
        bamIndex=sample[2],
        sampleName=sample[0], 
        RefFasta=refFasta, 
        RefIndex=refIndex, 
        RefDict=refDict,
        Memory_gb=memory_gb 
    }
  }
  call GenotypeGVCFs { 
    input: GVCFs=HaplotypeCallerERC.GVCF, 
        sampleName="CEUtrio", 
        RefFasta=refFasta, 
        RefIndex=refIndex, 
        RefDict=refDict,
        Memory_gb=memory_gb 
  }
}

task HaplotypeCallerERC {

  File RefFasta
  File RefIndex
  File RefDict
  Int Memory_gb
  String sampleName
  File bamFile
  File bamIndex

  command {
    gatk HaplotypeCaller \
        -ERC GVCF \
        -R ${RefFasta} \
        -I ${bamFile} \
        -O ${sampleName}_rawLikelihoods.g.vcf 
  }
  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
  }
  runtime {
    docker: "broadinstitute/gatk"
    memory_gb: Memory_gb
    memory: Memory_gb + "GB"
  }
}

task GenotypeGVCFs {
  #GenotypeGVCFs different in GATK4
  File RefFasta
  File RefIndex
  File RefDict
  Int Memory_gb
  String sampleName
  Array[File] GVCFs

  command {
    java -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R ${RefFasta} \
        -V ${sep=" -V " GVCFs} \
        -o ${sampleName}_rawVariants.vcf
  }
  output {
    File rawVCF = "${sampleName}_rawVariants.vcf"
  }
  runtime {
    docker: "broadinstitute/gatk3:3.8-1"
    memory_gb: "8"
    memory: "8GB"
  }
}
