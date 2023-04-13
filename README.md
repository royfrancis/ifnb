# ifnb

Build a docker image and push to dockerhub.

```
docker build -t shiny-ifnb .
docker tag shiny-ifnb royfrancis/shiny-ifnb:v1.0
docker tag shiny-ifnb royfrancis/shiny-ifnb:latest
docker push royfrancis/shiny-ifnb:v1.0
docker push royfrancis/shiny-ifnb:latest
```

DockerHub: https://hub.docker.com/r/royfrancis/shiny-ifnb
Serve: https://ifnb.serve.scilifelab.se/

