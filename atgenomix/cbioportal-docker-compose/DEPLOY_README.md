### prerequisite
- azure subscription with preconfigured azure storage account and fileshares serving as azure container instance 
  volumes. 
- [azcli](https://learn.microsoft.com/en-us/cli/azure/install-azure-cli)
- Docker Compose V2, preferably installed with [Docker Desktop](https://www.docker.com/blog/announcing-compose-v2-general-availability/)   

### docker compose  
Follow the following command to create docker context on azure, with reference to this [azure docs](https://learn.microsoft.com/en-us/azure/container-instances/tutorial-docker-compose)
```
az login
docker login azure
docker context create aci acienv
docker context use acienv
az acr login --name acienv
docker compose up
```
