#!/bin/bash
# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Run all scripts in sequence
# ==============================================================================

set -e

echo "=== Project 06: Clinical Dataset Engineering ==="
echo ""

cd "$(dirname "$0")"

echo "Step 1/6: Setup and dependencies..."
Rscript R/00_setup.R

echo ""
echo "Step 2/6: Extracting SDTM domains..."
Rscript R/01_extract_sdtm.R

echo ""
echo "Step 3/6: Building NONMEM dataset..."
Rscript R/02_build_nm_dataset.R

echo ""
echo "Step 4/6: Building ADNCA dataset..."
Rscript R/03_build_adnca.R

echo ""
echo "Step 5/6: Running independent QC..."
Rscript R/04_qc_nm_dataset.R

echo ""
echo "Step 6/6: Assembling submission package..."
Rscript R/06_assemble_submission.R

echo ""
echo "=== All steps complete ==="
echo "Outputs in: outputs/"
echo "Submission package in: submission_package_mock/"
echo ""
echo "To render QC report:"
echo "  Rscript -e \"rmarkdown::render('R/05_render_qc_report.Rmd', output_dir='outputs/')\""
echo ""
echo "To launch QC dashboard:"
echo "  Rscript -e \"shiny::runApp('shiny_app/')\""
