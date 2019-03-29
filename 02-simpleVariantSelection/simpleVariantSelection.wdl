workflow SimpleVariantSelection {
    File refFasta
    File refIndex
    File refDict    
    String name
    Int memory_gb

    call haplotypeCaller {
          input: sampleName=name,
                      RefFasta=refFasta,
                      RefIndex=refIndex,
                      RefDict=refDict,
                      Memory_gb=memory_gb
    }
    call select as selectSNPs {
          input: sampleName=name,
                      RefFasta=refFasta,
                      RefIndex=refIndex,
                      RefDict=refDict,
                      type="SNP",
                      rawVCF=haplotypeCaller.rawVCF,
                      Memory_gb=memory_gb
    }
    call select as selectIndels {
          input: sampleName=name,
                      RefFasta=refFasta,
                      RefIndex=refIndex,
                      RefDict=refDict,
                      type="INDEL",
                      rawVCF=haplotypeCaller.rawVCF,
                      Memory_gb=memory_gb
    }
}

task haplotypeCaller {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File inputBAM
  File bamIndex
  Int Memory_gb
  command {
    gatk HaplotypeCaller \
        -R ${RefFasta} \
        -I ${inputBAM} \
        -O ${sampleName}.raw.indels.snps.vcf
  }
  output {
    File rawVCF = "${sampleName}.raw.indels.snps.vcf"
  }
  runtime {
        docker: "broadinstitute/gatk"
        memory_gb: Memory_gb
        memory: Memory_gb + "GB"
  }     
}

task select {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  String type
  File rawVCF
  Int Memory_gb

  command {
    gatk SelectVariants \
      -R ${RefFasta} \
      -V ${rawVCF} \
      -select-type ${type} \
      -O ${sampleName}_raw.${type}.vcf
  }
  output {
    File rawSubset = "${sampleName}_raw.${type}.vcf"
  }
  runtime {
        docker: "broadinstitute/gatk"
        memory_gb: Memory_gb
        memory: Memory_gb + "GB"
  }     
}
