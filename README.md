<h1><i>Molektra: A Graphical User Interface for Continuous Shape and Symmetry Measures</i></h1>


<p align="center">
  <img src="gui/molektra_banner.png" width="650" alt="Molektra">
</p>

---
## About
  `Molektra` is a molecular visualizer for continuous shape measures and symmetry analysis, featuring built-in tools for rendering coordination polyhedra and symmetry elements. <br> <br>
Below, you will find instructions on how to download the pre-compiled executables from this repository and a quickstart to using `Molektra`.

<p align="left">
  <img src="gui/banner_2.png" width="350" alt="Installation">
</p>

### Option 1: Standalone Executables (Recommended)
Since `Molektra` is provided as a standalone executable, **no Python environment or additional libraries are required**.
  1. Go to XXX (../../XXX) page on the right side of this repository. <br> <br>
  2. Download the version corresponding to your operating system:
     - **Windows:** Download `Molektra_Windows.exe`
     - **macOS (Apple Silicon/ M-series):** Download `Molektra_macOS_arm.app`(for M1, M2, M3, M4 chips)
     - **macOS (Intel):** Download `Molektra_macOs_intel.app` (for Intel-based Macs) <br> <br>
  3. Double-click the downloaded file to launch the application. <br>
*(Note: Depending on your system settings, Windows SmartScreen or macOS Gatekeeper might prompt a security warning. You can safely bypass this by clicking "More info > Run anyway" on Windows, or by right-clicking the app and selecting "Open" on macOS).*

### Option 2: Run from Source code
If you prefer to inspect the code or run `Molektra` directly via Python, you can clone the repository and install the required dependencies manually. 
```bash
#Clone the repository
git clone [https://github.com/LoraSM/Molektra.git]
cd Molektra

#Install dependencies
pip install -r requirements_molektra.txt

#Launch the application
python main.py
```

<p align="left">
  <img src="gui/banner_3.png" width="350" alt="Quick Start">
</p>
