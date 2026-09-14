1) 

docker build -t ngs-seq-base -f .\Dockerfile .

2) 
docker run --rm -it -v "${PWD}:/work" ngs-seq-base

3) 
 docker exec -it jolly_allen bash