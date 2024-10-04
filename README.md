## scRNA-seq

## Introduction

## Workflow

This workflow is automated using a Nextflow workflow described in `sc_pipeline.nf`.

### Requirements

Running `sc_pipeline.nf` requires an environment with Nextflow installed, in addition to container management software such as Singularity or Charliecloud (see [Containers](#containers)). This workflow was developed and tested using Nextflow v24.04.4.


### Containers

All necessary scripts, third-party software, and dependencies are included in a Docker image hosted at https://hub.docker.com/repository/docker/apwat/10x_sc/. The Dockerfile and Conda environment file used to create the image are archived in this repository. Nextflow will attempt to pull and convert this image to one of the runtimes below, and then run its processes in containers 

#### Singularity

The default container software used to run processes in `sc_pipeline.nf`. 
Setting the environment variable `NXF_SINGULARITY_CACHEDIR` will control where images are downloaded.

```
Tested with SingularityCE 3.11 -- older versions of Singularity may be unable to properly set up containers.
```

#### Charliecloud

To instead use Charliecloud as a container runtime, comment out the `singularity` scope in `nextflow.config`, and uncomment the `charliecloud` scope.

As of Nextflow v24.04.4, workflows with multiple processes running in the same Charliecloud container will attempt multiple pulls of the same image, which fail. To resolve this, download the image ahead of time:

```
export CH_IMAGE_STORAGE=/path/to/image_storage/
export NXF_CHARLIECLOUD_CACHEDIR=/path/to/image_storage/
ch-image pull apwat/10x_sc:0.5
```
-and then run the workflow as usual.

```
Tested with Charliecloud 0.37.
```

#### Other Container Runtimes

See [the Nextflow documentation for containers](https://www.nextflow.io/docs/latest/container.html).


### Configuration Files

The workflow expects two configuration files in JSON format: 
1. A file listing samples and associated metadata. For an example, see [Example_Input.json](/reference_files/Example_Input.json). Locations for the data files must be given as absolute paths.
2. A file setting workflow parameters. For an example, see [Example_Input_Settings.json](/reference_files/Example_Input_Settings.json). All parameters are listed in the `params` scope of [nextflow.config](/workflow/nextflow.config). Most parameters are optional except for `sampleInfo`, which must be the path to the sample description file described in `1`. The use of absolute paths is also recommended here, as relative paths will be interpreted relative to the launch directory.

   
### Running the Workflow

To run the workflow: 

```
nextflow run /path/to/sc_pipeline.nf -params-file /path/to/Input_Settings.json
```
