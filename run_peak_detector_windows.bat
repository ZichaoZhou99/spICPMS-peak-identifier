@echo off
REM Short-dwell ICP-MS peak detector runner for Windows
REM Put this .bat file and spicpms_peak_detector.py in the same folder.

set /p INPUT=Drag and drop your CSV file here, then press Enter: 

python "%~dp0spicpms_peak_detector.py" %INPUT% --dwell-us 50 --baseline-points 1000 --background-method poisson --threshold-sigma 5 --min-bins 1 --merge-gap-bins 0

pause
