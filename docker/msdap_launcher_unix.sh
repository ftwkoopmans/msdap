#!/bin/sh
#
# MS-DAP launch script
# https://github.com/ftwkoopmans/msdap

VERSION="1.0.7"


### OS
ISMACOS=$(uname -s | grep -i "darwin")

## if macOS
if [ ! -z "$ISMACOS" ]; then
  echo "macOS"
  command -v docker >/dev/null 2>&1 || { echo >&2 "ERROR: Docker is not installed"; exit 1; }

  ## Open Docker if, and only if, it is not running; check exit code from `docker info`
  DOCKER_STATE=$(docker info >/dev/null 2>&1)
  if [ $? != 0 ]; then
    open /Applications/Docker.app
    echo "Waiting for Docker to launch."
    # loop until docker info gives a non-error exit code
    until docker info >/dev/null 2>&1
    do
      echo ".."
      sleep 1
    done
  fi
else
## if Linux
  echo "linux"
  # https://stackoverflow.com/a/677212
  # docker install guide: https://docs.docker.com/engine/install/#server  (in the table of supported platforms, click on your OS. eg; Ubuntu)
  command -v docker > /dev/null 2>&1 || { echo >&2 "ERROR: Docker is not installed"; exit 1; }

  if [ ! docker info > /dev/null 2>&1 ]; then
    set -x
    docker info groups
    set +x
    echo >&2 "ERROR: Docker daemon must be up-and-running prior to running this script"
    echo >&2 "    And the user must be a member of the group docker"
    exit 1
  fi
fi


### test if image is available locally. if not, docker pull
MSDAP_IMAGE_LOCAL = $(docker images -q ftwkoopmans/msdap:"$VERSION") >/dev/null 2>&1
if [ -z "$MSDAP_IMAGE_LOCAL" ]; then
  docker pull ftwkoopmans/msdap:"$VERSION"
  # exit script on docker pull error
  [ $? -eq 0 ] || exit 1
fi


### stop and remove all previous instances by image name
CONTAINER_ID_REMOVE=$(docker ps -a | awk '{ print $1,$2 }' | grep msdap | awk '{ print $1 }')
if [ ! -z "$CONTAINER_ID_REMOVE" ]; then
  # https://stackoverflow.com/a/32074098
  docker rm $(docker stop "$CONTAINER_ID_REMOVE")
fi


## if the lsof tool is available, run it (as sudo) and check whether the 3838 port has been claimed by third-party software. If so, throw an error
## disabled for now to remove reliance on external tools that may not be available, docker should show an appropriate error message anyway ## command -v lsof >/dev/null 2>&1 && (sudo lsof -i :3838 | grep -q -E 'LISTEN|ESTABLISHED' && { echo 'ERROR: port 3838 is already in use by some other software, cannot start MS-DAP (which would open on the same port)'; exit 1; })


### open browser after an N second delay
URL_RSTUDIO="http://localhost:3839"
echo -e $'\e[32m' "\nRStudio & MS-DAP will be available at: ${URL_RSTUDIO}\n" $'\e[0m'

if [ ! -z "$ISMACOS" ]; then
  sleep 3 && open "${URL_RSTUDIO}" &
else
  echo $'please open your web browser and navigate to the aforementioned URL to access MS-DAP (automatic launching of browser not supported on Linux atm.)\n'
fi

echo $'use control+C (or close this terminal window) to stop the container\n'


### docker run (mounting current directory on host as /data within container)
TIMEZONE_UTC_OFFSET=$(date +%z)
docker run -p 3839:8787 -e DISABLE_AUTH=true -v $(pwd):/data -it ftwkoopmans/msdap:"$VERSION" /bin/bash -c "echo HOST_TIMEZONE_UTC_OFFSET=$TIMEZONE_UTC_OFFSET >> /usr/local/lib/R/etc/Renviron && /init"
