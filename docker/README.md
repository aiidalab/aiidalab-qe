# Develop and maintain the docker images for the project

To build the docker images, run the following command (if you are using arm64 architecture, replace `linux/amd64` with `linux/arm64`):

```bash
docker buildx bake -f docker-bake.hcl --set qe.platform=linux/amd64 --load
```
