@echo off

REM This script applies an OpenMS workflow for label-free data analysis to all mzML files in the same directory using the fasta files on the same location.
REM workflow is the same as the iPRG2015 OpenMS tutorial by Hannes Rost
REM
REM requirements:
REM - OpenMS 2.5 or newer
REM - MSGFPlus installed as OpenMS thirdparty plugin
REM - Percolator installed as OpenMS thirdparty plugin
REM - java runtime is required for MSGFPlus, you probably have this but in case of problems/errors with MSGFPlus grab latest version from www.java.com
REM - after installation of the above, you may need to reboot first as we found on some test systems
REM
REM 1) place this .bat file in the directory where your input files are
REM 2) input: centroided mzML files  +  fasta files WITHOUT decoys (will be generated on-the-fly)
REM 3) output: mzTab output file compatible with MS-DAP  +  consensusXML for (optional) further processing by OpenMS
REM
REM
REM Version: 0.2
REM Licence: GPL-3.0-or-later
REM Author: Frank Koopmans
REM Email: frank.koopmans@vu.nl



REM settings below:
REM -----------------------------------------------------------------------------------------------------------------------

set NTHREADS=8
set "OPENMS=C:\Program Files\OpenMS-2.5.0"


REM end of settings
REM -----------------------------------------------------------------------------------------------------------------------

REM get timestamp. code from: https://stackoverflow.com/a/23476347
for /f "tokens=2 delims==" %%a in ('wmic OS Get localdatetime /value') do set "dt=%%a"
set "YY=%dt:~2,2%" & set "YYYY=%dt:~0,4%" & set "MM=%dt:~4,2%" & set "DD=%dt:~6,2%"
set "HH=%dt:~8,2%" & set "Min=%dt:~10,2%" & set "Sec=%dt:~12,2%"
set "LOGFILE=log_%YYYY%-%MM%-%DD%_%HH%-%Min%-%Sec%.txt"


echo start pipeline>%LOGFILE% 2>&1


set "MSGFPLUS=%OPENMS%\share\OpenMS\THIRDPARTY\MSGFPlus\MSGFPlus.jar"
set "PERCOLATOR=%OPENMS%\share\OpenMS\THIRDPARTY\Percolator\percolator.exe"

if not exist "%OPENMS%" (
  echo please update your OPENMS file path in the .bat file, directory not found; %OPENMS%
  pause
  exit /b
)

if not exist "%MSGFPLUS%" (
  echo cannot locate MSGFPLUS, which should be located at; %MSGFPLUS%
  pause
  exit /b
)

if not exist "%PERCOLATOR%" (
  echo cannot locate Percolator, which should be located at; %MSGFPLUS%
  pause
  exit /b
)

REM check if java is installed and available in current path (so downstream tools can find it). code from https://stackoverflow.com/a/33450403
java -version 1>nul 2>nul || (
   echo Not able to find Java executable or version. You can download and install from www.java.com
   exit /b 2
)

setlocal EnableDelayedExpansion


REM create a single target+decoy fasta file from all target fasta files in current dir
set /A HASFASTA=0
set CURDIR=%~dp0
set "FASTA_TARGET=!CURDIR!merged.fasta"
set "FASTA_TARGETDECOY=!CURDIR!merged_targetdecoy.fasta"

if exist !FASTA_TARGETDECOY! (
  set /A HASFASTA=1
  echo pre-existing "merged_targetdecoy.fasta" file found and used in this analysis
  echo pre-existing "merged_targetdecoy.fasta" file found and used in this analysis >>%LOGFILE% 2>&1

) else (

  REM collect all fasta files (full path, encapsulated by double-quotes) in space delimited string
  for %%G in (*.fasta) do (
    set /A HASFASTA=1
    REM from filename to full path
    set "FASTA=%%~dpG%%G"

    if defined FILES_FASTA (
      set FILES_FASTA=!FILES_FASTA! ^"!FASTA!^"
    ) else (
      set FILES_FASTA=^"!FASTA!^"
    )
  )

  REM merge fasta files and add decoys
  if defined FILES_FASTA (
    echo merge fasta files and add decoys; !FILES_FASTA!
    "%OPENMS%\bin\FileMerger.exe" -in !FILES_FASTA! -out "!FASTA_TARGET!" >>%LOGFILE% 2>&1
    "%OPENMS%\bin\DecoyDatabase.exe" -in "!FASTA_TARGET!" -out "!FASTA_TARGETDECOY!" >>%LOGFILE% 2>&1
  )

)

if HASFASTA==0 (
  echo no fasta files in current directory: please provide fasta files without decoys. this pipeline generates a target-decoy variant of each fasta file in current dir
  pause
  exit /b
)


REM experimental design
set EXPDESIGN="!%~dp0!mock_exp_design.tsv"
echo Fraction_Group	Fraction	Spectra_Filepath	Label	Sample>"!EXPDESIGN!"
set /a CNT = 0


REM create a target+decoy variant of each fasta file in current dir
set /A HASMZML=0
for %%G in (*.mzML) do (
  set /A HASMZML=1
  REM from filename to full path
  set "MZML=%%~dpG%%G"
  set "IDXML=%%~dpG%%~nG.idXML"

  echo processing file: %%G ...

  REM experimental design
  set /a CNT += 1
  echo !CNT!	1	!MZML!	1	!CNT!>>"!EXPDESIGN!"

  REM concatenate full filenames, each surrounded by quotes (in case the user's directories or filenames have spaces)
  if defined FILES_MZML (
    set FILES_MZML=!FILES_MZML! ^"!MZML!^"
  ) else (
    set FILES_MZML=^"!MZML!^"
  )

  if defined FILES_IDXML (
    set FILES_IDXML=!FILES_IDXML! ^"!IDXML!^"
  ) else (
    set FILES_IDXML=^"!IDXML!^"
  )

  REM apply OpenMS tools; search engine, percolator, fix confidence scores for compatability with downstream tools
  if not exist !IDXML! (
    "%OPENMS%\bin\MSGFPlusAdapter.exe" -in "!MZML!" -out "!IDXML!" -database "!FASTA_TARGETDECOY!" -executable "!MSGFPLUS!" -max_precursor_charge 5 -threads %NTHREADS% >>%LOGFILE% 2>&1
    "%OPENMS%\bin\PeptideIndexer.exe" -fasta "!FASTA_TARGETDECOY!" -in "!IDXML!" -out "!IDXML!" -enzyme:specificity none >>%LOGFILE% 2>&1
    "%OPENMS%\bin\PSMFeatureExtractor.exe" -in "!IDXML!" -out "!IDXML!" >>%LOGFILE% 2>&1
    "%OPENMS%\bin\PercolatorAdapter.exe" -in "!IDXML!" -out "!IDXML!" -percolator_executable "!PERCOLATOR!" -post-processing-tdc -subset-max-train 100000 >>%LOGFILE% 2>&1
    "%OPENMS%\bin\FalseDiscoveryRate.exe" -in "!IDXML!" -out "!IDXML!" -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins >>%LOGFILE% 2>&1
    "%OPENMS%\bin\IDFilter.exe" -in "!IDXML!" -out "!IDXML!" -score:pep 0.05 >>%LOGFILE% 2>&1
    "%OPENMS%\bin\IDScoreSwitcher.exe" -in "!IDXML!" -out "!IDXML!" -old_score q-value -new_score MS:1001493 -new_score_orientation lower_better -new_score_type "Posterior Error Probability" >>%LOGFILE% 2>&1
  ) else (
    echo !IDXML! already exists: skipping search engine for current mzML, assuming it was previously processed by this script
    echo !IDXML! already exists: skipping search engine for current mzML, assuming it was previously processed by this script >>%LOGFILE% 2>&1
  )
)

if HASMZML==0 (
  echo no mzML files in current directory: nothing to do...
  pause
  exit /b
)



REM apply OpenMS LFQ pipeline tool
echo running ProteomicsLFQ...
"%OPENMS%\bin\ProteomicsLFQ.exe" -in !FILES_MZML! -ids !FILES_IDXML! -Alignment:max_rt_shift 0.1 -fasta "!FASTA_TARGETDECOY!" -targeted_only "true" -transfer_ids "false" -mass_recalibration "false" -design "!EXPDESIGN!" -out_cxml "result.consensusXML" -out "result.mzTab" >>%LOGFILE% 2>&1


pause
exit /b
