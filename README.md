# ifnb

Build a docker image and push to dockerhub.

```
# build
docker build -t shiny-ifnb .

# set tags
docker tag shiny-ifnb royfrancis/shiny-ifnb:v1.0
docker tag shiny-ifnb royfrancis/shiny-ifnb:latest

# run app
docker run --rm -p 3838:3838 royfrancis/shiny-ifnb

# push to dockerhub
docker push royfrancis/shiny-ifnb:v1.0
docker push royfrancis/shiny-ifnb:latest
```

Locally, the app is accessible in the browser at [http://127.0.0.1:3838](http://127.0.0.1:3838).

DockerHub: https://hub.docker.com/r/royfrancis/shiny-ifnb
Serve: https://ifnb.serve.scilifelab.se/

