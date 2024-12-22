build: build-base build-final

build-base:
  docker build --tag linken_base:latest --file Dockerfiles/Dockerfile_base .
  
build-final:
  docker build \
    --tag linken:latest \
    --build-arg GIT_HASH=$(git rev-parse --short HEAD) \
    --build-arg BUILD_DATE=$(date +%Y-%b-%d) \
    --file Dockerfiles/Dockerfile_build .

clean-build:
  docker rmi linken_base:latest linken:latest
  docker builder prune
