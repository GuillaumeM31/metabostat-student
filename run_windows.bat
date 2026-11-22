@echo off
REM ========================================================
REM   MetaboStat / MSnet-Stat - Student Edition (Windows)
REM ========================================================
cd /d "%~dp0"
echo Starting MetaboStat Application...

where Rscript >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    Rscript run_app.R
) else (
    echo.
    echo NOTE: Rscript was not found in your system PATH.
    echo Trying default R installation directories...
    
    REM Search for standard 64-bit R installation
    for /f "delims=" %%i in ('dir /b /ad /o-n "C:\Program Files\R\R-*" 2^>nul') do (
        if exist "C:\Program Files\R\%%i\bin\Rscript.exe" (
            "C:\Program Files\R\%%i\bin\Rscript.exe" run_app.R
            goto :done
        )
    )

    echo.
    echo ERROR: Could not automatically locate Rscript.exe.
    echo To run the application:
    echo  1. Open RStudio
    echo  2. Open the file 'app.R' or 'run_app.R'
    echo  3. Click 'Run App' or run: source("run_app.R")
    echo.
    pause
)

:done
