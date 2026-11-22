#!/bin/bash
# MetaboStat Quick Launcher for macOS
cd "$(dirname "$0")"
echo "========================================================"
echo "  MetaboStat / MSnet-Stat - Student Edition (macOS)"
echo "========================================================"
echo "Starting application..."

if command -v Rscript &> /dev/null; then
    Rscript run_app.R
else
    echo "ERROR: Rscript not found in PATH."
    echo "Please ensure R is installed from https://cloud.r-project.org"
    echo "Or launch from RStudio by opening app.R and clicking 'Run App'."
    read -p "Press Enter to exit..."
fi
