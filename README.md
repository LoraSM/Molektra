<h1><i>Molektra: A Graphical User Interface for Continuous Shape and Symmetry Measures</i></h1>


<p align="center">
  <img src="gui/molektra_banner.png" width="650" alt="Molektra">
</p>

---
## About
  *Molektra* is a molecular visualizer for continuous shape measures and symmetry analysis, featuring built-in tools for rendering coordination polyhedra and symmetry elements. <br> <br>
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
git clone [https://github.com/LoraSM/Molektra.git](https://github.com/LoraSM/Molektra.git)
cd Molektra

#Install dependencies
pip install -r requirements_molektra.txt

#Launch the application
python main.py
```

<p align="left">
  <img src="gui/banner_3.png" width="350" alt="Quick Start">
</p>

Using `Molektra` is straightforward. Here is the basic workflow to analyze a molecule:
1. **Launch the app** go to `Load File > Open` to load your molecular structure (e.g., `.xyz`, `.pdb` or `.cif`), or you can simply drop the file into the visualization screen.<br><br>
2. **Select your atoms** by clicking on the central atom first, followed by the coordinating ligands. Alternatively, click the central atom and press **Cmd + E** (macOs) or **Ctrl+E** (Windows) to automatically select all neighboring atoms.<br><br>
3. **Run the continous shape measures (CShM) analysis**: Once your central atom is selected, press the `Shape Measures` button. `Molektra` automatically detects the coordinating neighbors to perform the calculations. The CShM values and their corresponding polyhedra will be displayed on a new screen. <br><br>
4. **Run the symmetry analysis**: First, detect the point group of your selection (mandatory step). Then, click the `Symmetry Measures` button to open a new screen detailing each symmetry element and the global symmetry analysis. <br><br>
5. **Export the results**: To save your work, press the `Export` button. You can choose to export the 3D visualization as an image or save the numerical data as a results file. <br><br>

### Full User Manual 
For a detailed, step-by-step explanation of all features and advanced visualization settings, please refer to the complete manual:
**[Read the Molektra User Manual](docs/Molektra_Manual.pdf)**
<br><br>

<p align="left">
  <img src="gui/additional_banner.png" width="650" alt="Additional Information">
</p>

### How to cite
If you use `Molektra` in your research, please cite our work: 

### Bug reports & Feedback
Found a bug or have a feature request? Please feel free to open an issue!
Go to the [Issues](../../issues) tab and describe the problem. Any feedback to improve `Molektra` is highly appreciated. 

### License
`Molektra` is open-source software licensed under the [MIT License](LICENSE).
