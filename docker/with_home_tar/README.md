## Docker image with home directory tar

When the user launches the container, a lot of time is spent on preparing the environment, setting up codes, pseudopotentials, etc. In order to reduce the launch time of the container, we tar already prepared home directory, and then untar it when the container is launched. This way, the user can start working with the container immediately.

```Dockerfile

To build the image, you can use the following command:

```bash
docker build -t aiidalab/qe-tar-home .
```

To run the container, you can use the following command:

```bash
docker run --rm -it -p 8888:8888 aiidalab/qe-tar-home
```

To compare with the image without the home directory, you can use the following command:

```bash
docker run --rm -it aiidalab/qe:amd64-latest /bin/bash
```




## Image size

| Type         | Data Size | Start Time  |
|--------------|-----------|-------------|
| Compress     | 5.21 GB   | 12 s        |
| No compress  | 6.4 GB    | 6 s         |


