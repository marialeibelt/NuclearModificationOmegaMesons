# Measuring the Nuclear Modification of Omega Meson Production in pp and OO Collisions with ALICE
## Description
Analysis framework for the measurement of nuclear modification (R<sub>OO</sub>) of ω → π⁺π⁻π⁰ production in pp and OO collisions at √s<sub>NN</sub> = 5.36 TeV with ALICE.

## About this project

This project was developed during my time as a CERN Summer Student 2025.

It contains the analysis for ω → π⁺π⁻π⁰ reconstruction and the measurement of the nuclear modification factor (R<sub>OO</sub>) in ALICE, including signal extraction, correction, and plotting.

### Overview of files

- **`omegarec.cpp`** – Extracts and processes histograms from raw data ROOT files. Performs background fitting and subtraction, then fits the ω signal in different p<sub>T</sub> bins.  
- **`plot_omega_compareflex.cpp`** – Corresponding plotting script. Compares signals extracted from different files and writes the resulting figures to a ROOT file.  
- **`omegagen.cpp`** – Extracts ω signals from generated minimum bias Monte Carlo (MC).  
- **`plot_omegagen.cpp`** – Plotting script for generated MC. Writes the corresponding histograms to a ROOT file.  
- **`omegaeff.cpp`** – Calculates the MC efficiency in chosen p<sub>T</sub> bins.  
- **`omegacorr.cpp`** – Produces the fully corrected ω spectra.  
- **`plot_omegacorr.cpp`** – Plotting script for corrected spectra, writing the figures to a ROOT file.  
- **`ROO.cpp`** – Calculates and writes the R<sub>OO</sub> values into a ROOT file.

## LICENCE
MIT License

Copyright (c) 2025 Maria Leibelt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
