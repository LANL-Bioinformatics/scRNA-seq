## Automated Nextflow Workflow

To run the workflow: 

```
nextflow run sc_pipeline.nf -params-file /path/to/WorkflowSettings.json
```

Your system will need Nextflow installed, as well as a way to run containers.

All necessary scripts, third-party software, and dependencies are included in a Docker image: https://hub.docker.com/repository/docker/apwat/10x_sc/. 

By default, the workflow runs all processes in a Singularity container pulled and converted from this image. You may wish to set the cache for Singularity images by setting NXF_SINGULARITY_CACHEDIR.