@echo off
if "%~1"=="" (
  probioscore --help
) else (
  probioscore %*
)
