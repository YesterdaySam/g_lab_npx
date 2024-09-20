
:: You can call CatGT three ways:
::
:: 1) > CatGT cmd-line-parameters
:: 2) > runit.bat cmd-line-parameters
:: 3a) Edit parameters in runit.bat, then call it ...
:: 3b) > runit.bat
::
:: This script effectively says:
:: "If there are no parameters sent to runit.bat, call CatGT
:: with the parameters hard coded here, else, pass all of the
:: parameters through to CatGT."
::

@echo off
@setlocal enableextensions
@cd /d "%~dp0"

set LOCALARGS=-dir=D:\Data\Kelton\probe_data\KW006 ^
-run=KW006_08092024_rec_D8_CA1 ^
-g=0 ^
-t=0 ^
-prb_fld ^
-lf ^
-lffilter=butter,12,0,300 ^
-prb=0 ^
-save=2,0,0,384


if [%1]==[] (set ARGS=%LOCALARGS%) else (set ARGS=%*)

%~dp0CatGT %ARGS%

