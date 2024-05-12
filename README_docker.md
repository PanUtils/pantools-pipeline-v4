# Build docker image
In folder where Dockerfile is located:
```bash
docker build -t pantools-pipeline-v4 .
```
    
# Run interactive container
Make sure you mount (-v "{folder-local-pc}:{folder-docker}") your configuration file, in- and output folders and your RAM-disk (/dev/shm/). Make sure the paths in your configuration file are pointing to the docker folder locations (/resources and /results in example)
```bash
docker run -v "/path/to/pantools_pipeline_v4_config.yaml:/pantools-pipeline-v4/config/config.yaml" -v "/path/to/pantools-pipeline-v4/resources:/resources/" -v "/path/to/results:/results/" -v "/dev/shm:/dev/shm" -it pantools-pipeline-v4 bash
```
This creates an interactive container where you can run all snakemake commands as described in the README for the pipeline. Type ```exit``` when you want to stop the container.