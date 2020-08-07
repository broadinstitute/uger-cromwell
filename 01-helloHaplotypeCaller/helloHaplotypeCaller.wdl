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
  String memory
  Int cpu
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
        docker: "broadinstitute/gatk@sha256:8051adab0ff725e7e9c2af5997680346f3c3799b2df3785dd51d4abdd3da747b"
        memory: memory
        cpu: cpu
  }

}
