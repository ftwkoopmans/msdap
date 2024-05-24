# MS-DAP launch script
# https://github.com/ftwkoopmans/msdap

$VERSION = "1.0.8"


Write-Host "$((Get-Date).ToString("HH:mm:ss")) - Starting MS-DAP launcher script for version: $VERSION"


### test if docker is installed
$dd_path = "C:\Program Files\Docker\Docker\Docker Desktop.exe"
$dd_is_installed = Test-Path $dd_path -PathType Leaf
if(!$dd_is_installed) {
  Write-Host "ERROR: Docker Desktop is not installed; please download from https://www.docker.com . After installation with default settings, this is where we expect Docker Desktop to be found: ${dd_path}"
  pause
  exit 1
}


### start docker daemon if needed
$dd_is_running = Get-Process 'com.docker.proxy' -ErrorAction SilentlyContinue
if($dd_is_running -eq $null) {
  Write-Host "$((Get-Date).ToString("HH:mm:ss")) - Starting Docker Desktop"
  & $dd_path
  Write-Host -NoNewline "$((Get-Date).ToString("HH:mm:ss")) - Waiting until docker daemon is ready.."
  Start-Sleep -Seconds 5
} else {
  Write-Host -NoNewline "$((Get-Date).ToString("HH:mm:ss")) - Waiting until docker daemon is ready.."
}


### wait until docker is ready
$counter = 0
while($true) {
  (docker info) 2>&1 | Out-Null
  if($LASTEXITCODE -eq 0) {
    Write-Host ' done'
    break
  } else {
    # print a . every N seconds to signal the script is not hanging
    if($counter%5 -eq 0) {
      Write-Host -NoNewline '.'
    }
    Start-Sleep -Seconds 1
  }
  $counter++
}


### test if image is available locally. if not, docker pull
$msdap_image_local = (docker images -q ftwkoopmans/msdap:$VERSION)
if($msdap_image_local -eq $null) {
  try {
    (docker pull ftwkoopmans/msdap:$VERSION)
  } catch {
    Write-Host "An error occurred @ docker pull:"
    Write-Host $_
    pause
	Exit 1
  }

  if($LastExitCode -ne 0) {
    pause
    Exit $LastExitCode
  }
}


### stop and remove all previous instances by image name
$container_id_remove = (docker ps -a | findstr -i msdap | % {$_ -replace " .*", ""})
# docker ps -a | awk '{ print $1,$2 }' | grep msdap | awk '{print $1 }'
if($container_id_remove -ne $null) {
  Write-Host "$((Get-Date).ToString("HH:mm:ss")) - Stopping up-and-running MS-DAP containers (if any) ..."
  try {
    # https://stackoverflow.com/a/32074098
    (docker rm $(docker stop $container_id_remove)) 2>&1 | Out-Null
  } catch {
    Write-Host "An error occurred @ docker stop:"
    Write-Host $_
  }
}


### throw error if the port is already taken (since we previously closed down all msdap containers, who's using it?)
### this is optional since docker run will perform a similar check, but this approach enables a usecase-specific error message
try {
  # rstudio server port
  $tcp = new-object System.Net.Sockets.TcpClient
  $tcp.ReceiveTimeout = 100
  $tcp.SendTimeout = 100
  $tmp = ($tcp.Connect("localhost", 3839)) 2>&1 | Out-Null
  $tcp.close()

  # if we reach this without error, there is an open TCP port where we want to open MS-DAP. throw error
  Write-Host ""
  Write-Host "ERROR: port 3839 is already in use by some other software, cannot start MS-DAP (which would open on the same port)"
  Write-Host ""
  Write-Host "try to manually restart Docker. If the problem persists, restart your computer"
  pause
  exit 1
} catch {
}


### open browser after an N second delay
$url_rstudio = "http://localhost:3839"
Write-Host ""
Write-Host "$((Get-Date).ToString("HH:mm:ss")) - RStudio & MS-DAP will be available at: ${url_rstudio}"
Write-Host ""
Write-Host "use control+C (or close this powershell window) to stop the container"
Write-Host ""
Write-Host ""
# start job background, so we must pass the URL parameter to the new powershell process. ref; https://stackoverflow.com/a/26240244
Start-Job {param($path) Start-Sleep 3; Start-Process $path } -Arg $url_rstudio | Out-Null


Write-Host "$((Get-Date).ToString("HH:mm:ss")) - Starting MS-DAP container ..."

### docker run (mounting current directory on host as /data within container)
$timezone_utc_offset = [System.TimeZone]::CurrentTimeZone.GetUtcOffset([datetime]::Now).TotalHours
docker run -p 3839:8787 -e DISABLE_AUTH=true -v ${PWD}:/data -it ftwkoopmans/msdap:$VERSION /bin/bash -c "echo HOST_TIMEZONE_UTC_OFFSET=$timezone_utc_offset >> /usr/local/lib/R/etc/Renviron && /init"



# don't close until user hits enter key (eg; when double-clicking the script in Windows Explorer, we don't want the powershell window to automatically close)
pause
