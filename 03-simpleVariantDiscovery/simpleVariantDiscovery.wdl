workflow SimpleVariantDiscovery {
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
    call hardFilterSNP {
          input: sampleName=name, 
                      RefFasta=refFasta, 
                      RefIndex=refIndex, 
                      RefDict=refDict, 
                      rawSNPs=selectSNPs.rawSubset,
                      Memory_gb=memory_gb
    }
    call hardFilterIndel {
          input: sampleName=name, 
                      RefFasta=refFasta, 
                      RefIndex=refIndex, 
                      RefDict=refDict, 
                      rawIndels=selectIndels.rawSubset,
                      Memory_gb=memory_gb
    }
    call combine {
          input: sampleName=name, 
                      RefFasta=refFasta, 
                      RefIndex=refIndex, 
                      RefDict=refDict, 
                      filteredSNPs=hardFilterSNP.filteredSNPs, 
                      filteredIndels=hardFilterIndel.filteredIndels,
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

task hardFilterSNP {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File rawSNPs
  Int Memory_gb

  command {
    gatk VariantFiltration \
      -R ${RefFasta} \
      -V ${rawSNPs} \
      --filter-expression "FS > 60.0" \
      --filter-name "snp_filter" \
      -O ${sampleName}.filtered.snps.vcf
  }
  output {
    File filteredSNPs = "${sampleName}.filtered.snps.vcf"
  }
  runtime {
        docker: "broadinstitute/gatk"
        memory_gb: Memory_gb
        memory: Memory_gb + "GB"
  }
}

task hardFilterIndel {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File rawIndels
  Int Memory_gb

  command {
    gatk VariantFiltration \
      -R ${RefFasta} \
      -V ${rawIndels} \
      --filter-expression "FS > 200.0" \
      --filter-name "indel_filter" \
      -O ${sampleName}.filtered.indels.vcf
  }
  output {
    File filteredIndels = "${sampleName}.filtered.indels.vcf"
  }
  runtime {
        docker: "broadinstitute/gatk"
        memory_gb: Memory_gb
        memory: Memory_gb + "GB"
  } 
}

task combine {
  #needs to run in gatk3 because CombineVariants is not in gatk4
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File filteredSNPs
  File filteredIndels
  Int Memory_gb

  command {
     java -jar /usr/GenomeAnalysisTK.jar \
      -T CombineVariants \
      -R ${RefFasta} \
      -V ${filteredSNPs} \
      -V ${filteredIndels} \
      --genotypemergeoption UNSORTED \
      -o ${sampleName}.filtered.snps.indels.vcf
  }
  output {
    File filteredVCF = "${sampleName}.filtered.snps.indels.vcf"
  }
  runtime {
        docker: "broadinstitute/gatk3:3.8-1"
        memory_gb: Memory_gb
        memory: Memory_gb + "GB"
  } 
}
