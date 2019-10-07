# Build and Run OPS Simulations in Docker

### Pull the image from Dockerhub
```docker
docker image pull amit112amit/alpine-ops:1.0`
```

### Sample Run Command
1. Create environment variable PYTHONPATH
2. Mount local directory inside the container. The local directory contains the driver and the input files.
3. Set the working directory inside the container
4. Run OPS simulation

```docker
docker run -e PYTHONPATH=/usr/lib/python3.7/site-packages:/simulation \
	-v /home/amit/WorkSpace/UCLA/simulations/Trial:/simulation/run \
	-w /simulation/run \
	amit112amit/alpine-ops:1.0 \
	python phasediagramsimulation.py -i 1
```
