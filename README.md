# Description

This repository contains a Nextflow implementation of the [CgpWGS_hg38 pipeline](https://github.com/panprostate/CgpWGS_hg38/tree/master).

## **Requirements**
All tools used in this pipeline are contained in the [CgpWGS_hg38 container](https://github.com/cancerit/dockstore-cgpwgs). It can be excuted via Singularity or Docker.


- Nextflow >= 23.04 (tested version - This script makes use of DSL 2 synthax, which became the default at 22.03.0-edge)
- [Reference files](https://github.com/panprostate/CgpWGS_hg38/blob/master/docs/sanger_pipeline.md) (As of now, the script will not check for these files and download them.)


**Please do not install Singularity through conda**.

Singularity must be owned by `root` in order to run this workflow. If you do not have access to a singularity installation owned by `root` [try requesting to your system administrator](https://sylabs.io/guides/3.7/user-guide/quick_start.html#installation-request).


## **How to run**
1. Setup a working directory and copy the `cgpWGS.nf` and `nextflow.config` files inside it.
2. Edit the "params{}" section in `nextflow.config` file to point to the appropriate directories (similar to the `nextflow.config` of the original pipeline). NOTE: Paths must either be absolute, or include the variable `$projectDir` if using paths relative to the script directory.
3. Edit the "process{}" section in `nextflow.config` file to adjust executor settings.
4. Create an input CSV table. The table must have 2 columns only: tumor file and normal file paths e.g. [./mytumor.bam,./mynormal.bam].
5. Execute `nextflow run cgpWGS.nf -with-singularity <path-to-container>` to execute the pipeline.

## **Tips**

- You can execute nextflow with the flag `-with-report` to create an HTML resource usage report. This is useful to set resource configuration.
- If your pipeline execution is interrupted/fails, you can execute nextflow with the flag `-resume` to resume where you left off.
