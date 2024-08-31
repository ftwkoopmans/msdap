
# About Docker

The MS-DAP R package depends on many other tools / code libraries, so in
order to work with MS-DAP one needs to install all these dependencies as
well as the R programming language and RStudio (ref; MS-DAP R package
install guide). So there’s some effort involved in installing all
software and dependencies. Furthermore, suppose that one of the
dependencies (e.g. an R package for a DEA algorithm) has an impactful
change in version X and on one computer you use MS-DAP with version X of
that library and on another version Y; running the same data analysis on
both computers might not yield the exactly same results.

This is where “containers” can be useful; these are isolated
environments that include (in 1 big bundle) all code and dependencies
needed to run an application on any computing environment (i.e. the
container doesn’t even rely on any software/tools present in the
operating system, it’s completely independent). We provide a container
for MS-DAP that includes the MS-DAP R package, all dependencies, RStudio
and the R programming language all neatly packed together. On your
computer, only the Docker software is needed to use the MS-DAP
container.

Advantages of using the Docker container:

- perfect reproducibility of data analysis results, the entire software
  stack is always the same
- only tool that needs to be installed on your computer is Docker

Disadvantages of using the Docker container:

- if you’re already working with R/RStudio, getting up-and-running by
  just installing the MS-DAP R package is much easier
- there’s a learning curve to use Docker, even though we provide
  launcher scripts to make things easier you may need to consult with
  your local IT support if things don’t go smoothly
  - dealing with files / sharing files between computer and Docker
    container can be especially tricky
- performing additional analyses outside standard MS-DAP workflows might
  take some extra work if you want to use R packages that aren’t
  available within the Docker container
  - i.e. you don’t want to install new/extra R packages into the MS-DAP
    container as you might mistakingly update other packages that MS-DAP
    depends on, thereby breaking the reproducibility aspect of working
    with containers

*note; this is intended as a very brief introduction and only highlights
some common (dis)advantages from the perspective of MS-DAP users. This
is by no means a complete introduction to containers or Docker, nor a
thorough examination of the pros and cons of working with containers.*

# Using the Dockerized MS-DAP

## install Docker for Windows or macOS

- download and install Docker for Windows/macOS at:
  <https://www.docker.com/products/docker-desktop>
  - *windows; If WSL2 has not not been installed on your system yet, a
    popup notification with installation instructions may appear. Follow
    the presented link, or [click
    here](https://docs.microsoft.com/en-us/windows/wsl/wsl2-kernel) to
    “download the latest WSL2 Linux kernel”. Finally, install this file
    and reboot.*
- configure Docker
  - start Docker Desktop
  - open the settings screen
    - windows; right-click the Docker icon in the taskbar at
      bottom-right of the screen (looks like a whale, you may need to
      click the little arrow on the windows taskbar if icons are hidden)
      \> click settings
    - macOS; Docker \> menu bar \> Preferences
  - update settings
    - resources \> advanced: set number of CPU cores to use (if in
      doubt, set to 4 or more)
    - resources \> advanced: set amount of RAM to use (if in doubt, set
      to 4~8GB)
    - resources \> file sharing: click the + symbol and select the
      directory that holds all your experimental data (eg; C:/data)
- reboot your computer after installing Docker Desktop

## install and test Docker on Linux

Follow the official installation guide for your Linux distribution at:
<https://docs.docker.com/engine/install/>

Example of running the MS-DAP Docker container that was tested on a
Ubuntu live-USB;

1)  Installation following the guide on
    <https://docs.docker.com/engine/install/ubuntu/> and using the
    “convenience script” for installation provided by Docker;

- `curl -fsSL https://get.docker.com -o get-docker.sh`
- `sudo sh get-docker.sh`

2)  Setting up users following
    <https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user>

- `sudo groupadd docker` \# might get a warning that group already
  exists, which is fine
- `sudo usermod -aG docker $USER` \# or instead of `$USER` , enter your
  username
- `newgrp docker` \# update group
- Now we’ll need to log out (the current user) and log back in so the
  group membership is re-evaluated. (in the Ubuntu live-USB example this
  could be done by; `sudo -i` then `su ubuntu` )
- `docker run hello-world` \# test if docker works. Should see Docker
  downloading an image, print “Hello from Docker !” among other
  info/text

3)  run the MS-DAP docker launcher script (you can download it from this
    GitHub repository, there’s a versioned script attached to each
    release at <https://github.com/ftwkoopmans/msdap/releases> )

- `mkdir -m 777 ~/msdap_data` \# importantly, check the permission flags
  used here. Create a folder in your home directory where data is
  exchanged between host system and the Docker container
- `cd ~/msdap_data` \# we want to run the launcher script from within
  the directory where data is exchanged
- `~/Downloads/msdap_launcher_unix.sh` \# assume the launcher script is
  here and already made executable

## Start the Docker container

This section shows 2 different approaches to start the MS-DAP Docker
container; either use the all-in-one launcher script OR execute a few
commands on the commandline. For typical use-cases, we recommend using
the launcher script.

### Launcher script

A launcher script is provided for user convenience; it executes all
commands needed to run the MS-DAP docker container. To do this manually
instead, consult the next section for a step-by-step guide.

- download the script for the MS-DAP version you want to use at the
  [‘releases’ section of this github
  repository](https://github.com/ftwkoopmans/msdap/releases)
- copy/move the launcher script to the directory that holds your
  experiment data (eg; `C:/data/proteomics/`)
  - the file location is important: only data with the directory that
    contains the launcher script will be accessible within the Docker
    container. See further, subsection **Locating your mounted data
    within Docker containers**
- run the script to launch MS-DAP
  - windows; right-click `msdap_launcher_windows.ps1` and select “Run
    with PowerShell”
    - *note; upon first use, you may see a security warning because this
      script was downloaded from the internet (If it’s a “Execution
      Policy Change” warning, enter “Y” to proceed)*
  - macOS and Linux; execute the `msdap_launcher_unix.sh` script
    - *note; if you are unfamiliar with using scripts,* [here is a macOS
      guide](https://support.apple.com/guide/terminal/make-a-file-executable-apdd100908f-06b3-4e63-8a87-32e71241bab4/mac)
    - *note; on Linux, you need to make sure the Docker daemon is
      running prior to starting the script*

To use a different MS-DAP version, simply change the version number at
the top of the launcher script. To view all available MS-DAP Docker
containers versions, go to:
<https://hub.docker.com/r/ftwkoopmans/msdap/tags>

**what does the script do?**

- if Docker is not up-and-running, start Docker
- the launcher script will automatically download the MS-DAP container
  (if needed)
- if MS-DAP is already up-and-running, stop it
- start the MS-DAP docker container
- your default web browser is launched and navigates to a local
  webserver within the Docker container @ <http://localhost:3839>

### commandline

This section details each step involved in using the MS-DAP Docker
container. Instead of following this multi-step process, you may find it
easier to use the *launcher script* described in the previous section
(and consider the below guide as a fall-back in case of technical issues
with the launcher script).

1)  First, make sure that Docker is up and running. For Windows and
    macOS, simply start the Docker Desktop application.

- On Windows, you may hover the mouse cursor over the Docker tray icon
  that looks like a whale in the bottom-right section of the taskbar to
  confirm its status is “Docker Desktop is running” (if no icon is
  visible but you did start Docker, maybe the icon is hidden and you
  have to click the icon on the bottom-right of the Windows taskbar to
  “Show hidden icons”).

2)  Next, start a terminal. On Windows, start powershell (start \> type:
    `powershell` \> click `Windows PowerShell` icon). On macOS, open
    Terminal.

3)  Download the MS-DAP Docker container by executing the following
    command in the terminal (only have to do this ‘installation step’
    once):

`docker pull ftwkoopmans/msdap:1.0.3`

4)  Run the Docker container and give it access to the directory on your
    computer (the ‘host’ system) that contains the proteomics data.
    **importantly, Docker containers only have access to files on the
    host system that are explicitly mounted: update the highlighted path
    in the following command** (see further @ next subsection).

<code>docker run -p 3839:8787 -e PASSWORD=msdap -v
**C:/data/proteomics**:/data -it ftwkoopmans/msdap:1.0.3</code>

5)  Open a browser, such as Firefox or Chrome, and access MS-DAP at
    <http://localhost:3839> . You will be presented with a complete
    RStudio environment that contains MS-DAP!

note: to stop the container, go to the terminal used to launch MS-DAP
and press `control + C`.

note: to use a different MS-DAP version, simply change the version
number in both docker commands. To view all available MS-DAP Docker
containers, go to:
<https://hub.docker.com/repository/docker/ftwkoopmans/msdap>

## Locating your mounted data within Docker containers

Docker containers only have access to files on the host system that are
explicitly mounted. So if for example all proteomics datasets are
located in `C:/data/proteomics` on your computer (the ‘host’ system), we
will mount this directory in the Docker container to make its contents
available to the software running inside the container at file location
`/data` (as seen within the container).

**example: your data is located in C:/data/proteomics/…**

| Host system                   | as seen within container |
|-------------------------------|--------------------------|
| C:/data/proteomics            | /data                    |
| C:/data/proteomics/exp1       | /data/exp1               |
| C:/data/proteomics/exp1/fasta | /data/exp1/fasta         |

For reference, the syntax of this command:
`docker run ... -v <FULL PATH ON YOUR SYSTEM>:/data ...` (the `/data` is
the file location as seen from *within* the container and should not be
altered).

### Quickstart: an example R script

Use a snippet of R code that uses the MS-DAP R package to analyse the
Klaassen et al. (2018) dataset that was pre-processed with MetaMorpheus
and included with MS-DAP:

- To immediately view the expected results, [click here to download the
  PDF
  report](/examples/data/dataset_Klaassen2018_pmid26931375_report.pdf)
  we already prepared and uploaded.

Data analysis steps:

- assuming the MS-DAP Docker container is running and you accessed it
  using a web-browser at <http://localhost:3839>
- in the RStudio menu on the top-left; file \> New File \> R Script
- copy/paste the below code snippet to the panel in center of the screen
- run all lines of code: select all (by mouse drag or click anywhere
  then `control+A`) \> click the `run` button on top of the code panel
  to run all selected lines
- output log and result summary are shown at the bottom of the screen
- output files can be found in the *mapped directory* that you defined
  when launching MS-DAP

``` r
library(msdap)

# set working directory. use forward slashes in file paths
setwd("/exampledata/dataset_Klaassen2018_pmid26931375")
# load the dataset
dataset = import_dataset_metamorpheus(path = getwd(), protein_qval_threshold = 0.05)
# gather protein information, such as gene symbols, from fasta files
dataset = import_fasta(dataset, files = c("UP000000589_10090.fasta", "UP000000589_10090_additional.fasta"))
# optionally, we can apply filters to protein names. here we remove keratin and IGGs from the IP dataset
dataset = remove_proteins_by_name(dataset, regular_expression = "ig \\S+ chain|keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)")
# import sample metadata and define a contrast for DEA
dataset = import_sample_metadata(dataset, filename = "sample_metadata.xlsx")
dataset = setup_contrasts(dataset, contrast_list = list(c("wt", "ko")))
# run the pipeline. check the documentation by running the R command: ?analysis_quickstart
dataset = analysis_quickstart(dataset,
                              filter_min_detect = 1,
                              filter_min_quant = 3,
                              filter_by_contrast = TRUE,
                              norm_algorithm = c("vwmb", "modebetween_protein"),
                              dea_algorithm = c("deqms", "msempire", "msqrob"),
                              dea_qvalue_threshold = 0.05,
                              dea_log2foldchange_threshold = NA,
                              diffdetect_min_samples_observed = 2, # for differential detection only
                              output_qc_report = TRUE,
                              multiprocessing_maxcores = 4,
                              # note; within the Docker container, /data maps to the directory on your computer
                              # as defined when launching the MS-DAP Docker container
                              output_dir = "/data/dataset_Klaassen2018_pmid26931375",
                              output_within_timestamped_subdirectory = TRUE)
# optionally, print a brief summary to console
print_dataset_summary(dataset)
```

### troubleshooting

**documentation from Docker to verify your Docker installation is ok**

- <https://docs.docker.com/docker-for-windows/>
- <https://docs.docker.com/docker-for-mac/>

**permission denied while trying to connect to the Docker daemon**

The current user is not privileged to use docker. For instance, docker
daemon is running as user U from group G and either; a) user U is root,
or b) current user is not a member of group G (i.e. current user should
be added to user group G).

The docker installation guide shows how to run docker as non-root user;
<https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user>

**failed to create output director / reason ‘Permission denied’**

When using the MS-DAP Docker container, you need to make sure the
directory from which the launcher script is started is configured such
that the current user has read and write permissions. For example,
create a directory with world-writable permissions and start the
launcher script from within. See Ubuntu example above.

**cannot locate input data from within the Docker container**

1)  double-check that the *mapped directory* is correctly provided when
    starting the Docker container, as detailed in the above
    documentation.

2)  within the MS-DAP bundled RStudio at <http://localhost:3839> ; run
    the following R command in the console (on the bottom of the RStudio
    interface) to list all files in the *mapped directory* `/data` to
    verify whether your data is visible from within the container:
    `dir("/data")`

**timestamped MS-DAP output directories do not match the current time**

This is a known issue; the timezone of the Docker container (and R
runtime within) may not match the host computer. We try to deal with
this in the launcher scripts, but issues may persist.
