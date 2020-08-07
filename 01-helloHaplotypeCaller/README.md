# Basic WDL/Cromwell on UGER
We will be recreating the first WDL tutorial here: https://software.broadinstitute.org/wdl/documentation/article?id=7158 Please make sure to download and extract the zip data bundle to somewhere on shared NFS storage. To save you the trouble of downloading the zip from Google Drive and scp'ing it a Broad host a copy has been put in `/broad/data/wdl/tutorial/helloHaplotypeCaller.zip`. You can copy this to your home directory or hptmp and `unzip`.


Following the WDL tutorial we will generate an inputs.json file to tell cromwell where to look for our files, and we will validate our WDL for any syntax errors. Complete documentation for this process can be found in the previously linked WDL tutorial and the womtool documentation here: https://cromwell.readthedocs.io/en/develop/WOMtool/

Assuming all of the steps before were followed, we are currently sitting in this cloned repo somwhere on NFS storage on a UGER compute node. We then run the following to generate an inputs.json file for the first tutorial example:

```
java -jar womtool-52.jar inputs 01-helloHaplotypeCaller/helloHaplotypeCaller.wdl > 01-helloHaplotypeCaller/helloHaplotypeCaller_inputs.json
```

Editing 01-helloHaplotypeCaller/helloHaplotypeCaller_inputs.json with our favorite editor we can now specify our input file locations as well as how much memory UGER and PAPIvW will use when running this job.

An example using hptmp could look like the following:

```
{
  "helloHaplotypeCaller.haplotypeCaller.RefFasta": "/broad/hptmp/<your_username>/helloHaplotypeCaller/ref/human_g1k_b37_20.fasta",
  "helloHaplotypeCaller.haplotypeCaller.bamIndex": "/broad/hptmp/<your_username>/helloHaplotypeCaller/inputs/NA12878_wgs_20.bai",
  "helloHaplotypeCaller.haplotypeCaller.RefIndex": "/broad/hptmp/<your_username>/helloHaplotypeCaller/ref/human_g1k_b37_20.fasta.fai",
  "helloHaplotypeCaller.haplotypeCaller.RefDict": "/broad/hptmp/<your_username>/helloHaplotypeCaller/ref/human_g1k_b37_20.dict",
  "helloHaplotypeCaller.haplotypeCaller.inputBAM": "/broad/hptmp/<your_username>/helloHaplotypeCaller/inputs/NA12878_wgs_20.bam",
  "helloHaplotypeCaller.haplotypeCaller.sampleName": "NA12878"
  "helloHaplotypeCaller.haplotypeCaller.memory": "4 GB",
  "helloHaplotypeCaller.haplotypeCaller.cpu": "1"
}
```

We are now ready to submit our first WDL using cromwell to UGER.

```
$ java -Dconfig.file=uger.conf -jar cromwell-52.jar run 01-helloHaplotypeCaller/helloHaplotypeCaller.wdl -i 01-helloHaplotypeCaller/helloHaplotypeCaller_inputs.json
```

You should see a bunch of cromwell logs log to stdout and eventually lines similar to the following if everything worked correctly:
```
[2019-03-29 03:28:23,26] [info] DispatchedConfigAsyncJobExecutionActor [b9928bedhelloHaplotypeCaller.haplotypeCaller:NA:1]: job id: 9390797
[2019-03-29 03:28:23,27] [info] DispatchedConfigAsyncJobExecutionActor [b9928bedhelloHaplotypeCaller.haplotypeCaller:NA:1]: Status change from - to Running
[2019-03-29 03:44:27,03] [info] DispatchedConfigAsyncJobExecutionActor [b9928bedhelloHaplotypeCaller.haplotypeCaller:NA:1]: Status change from Running to Done
[2019-03-29 03:44:28,91] [info] WorkflowExecutionActor-b9928bed-2e7c-4f3c-b852-44b9f79e3404 [b9928bed]: Workflow helloHaplotypeCaller complete. Final Outputs:
{
  "helloHaplotypeCaller.haplotypeCaller.rawVCF": "uger-cromwell/cromwell-executions/helloHaplotypeCaller/b9928bed-2e7c-4f3c-b852-44b9f79e3404/call-haplotypeCaller/execution/NA12878.raw.indels.snps.vcf"
}
[2019-03-29 03:44:28,95] [info] WorkflowManagerActor WorkflowActor-b9928bed-2e7c-4f3c-b852-44b9f79e3404 is in a terminal state: WorkflowSucceededState
[2019-03-29 03:44:37,92] [info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "helloHaplotypeCaller.haplotypeCaller.rawVCF": "uger-cromwell/cromwell-executions/helloHaplotypeCaller/b9928bed-2e7c-4f3c-b852-44b9f79e3404/call-haplotypeCaller/execution/NA12878.raw.indels.snps.vcf"
  },
  "id": "b9928bed-2e7c-4f3c-b852-44b9f79e3404"
}
```
