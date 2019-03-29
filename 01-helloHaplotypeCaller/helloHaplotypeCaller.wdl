workflow helloHaplotypeCaller {
  call haplotypeCaller
}

task haplotypeCaller {
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File inputBAM
  File bamIndex
  Int memory_gb
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
        memory_gb: memory_gb
        memory: memory_gb + "GB"
  }

}
