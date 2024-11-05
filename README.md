# Spectral Efficiency and Tandem Performance Calculators
* The Spectral Efficiency (SE) calculator takes a set of current vs voltage (I-V) and quantum efficiency (QE) data along with a spectrum input and plots the spectral efficiency of the cells represented by the datasets for the spectrum provided. 

* The Tandem Performance calculator uses the calculated spectral efficiency of top and bottom cells and the transmission data of the top cell to calculate performance metrics for a tandem device in either a two or four terminal architecture.

* Note: The calculators are designed for cell-level calculations. For best results, module-level data should be adjusted to correspond to a single constituent cell (i.e. module voltage should be divided by the number of cells to show average cell voltage).

The SE and Tandem Performance calculators are based on the following publication:
  Yu, Z., Leilaeioun, M. & Holman, Z. Selecting tandem partners for silicon solar cells. Nat Energy 1, 16137 (2016). https://doi.org/10.1038/nenergy.2016.137

# Download and Install

* Recommended: The SE and Tandem Calculators run most reliably in a conda virtual environment. To set up a virtual environment, first install conda from https://conda.io/projects/conda/en/latest/user-guide/install/index.html. Once conda is installed, use 'conda env create -n \<name>' replacing \<name> with your environment name to create a new virtual environment. Activate your environment with 'conda activate \<name>' and install the packages listed below using 
* 'conda install colorama ipywidgets  matplotlib numpy pandas traitlets cycler ipython scipy openpyxl' followed by 
* 'pip install iv-params os-sys textwrap3 tk pytest-warnings datetime'.

* Clone or download the code for the SE and Tandem Performance calculators here: https://github.com/NREL/SE-and-Tandems

* Jupyter Lab is required for running these calculators. Jupyter Lab can be installed using "pip install jupyterlab" and run using "jupyter lab" (documentation: https://jupyter.org/install). Using Jupyter Notebook in place of Lab may cause errors with data upload and visualization.

* Please follow the tutorial below to use the Spectral Efficiency and Tandem Performance calculators. Example files including proper formatting of these files can be found in "SE-and-Tandems/Examples."

Packages needed for Spectral Efficiency calculator 
--------------
* os
* sys
* datetime
* tkinter
* colorama
* ipywidgets 
* matplotlib
* numpy
* pandas
* traitlets
* cycler
* ipython
* scipy
* openpyxl
* iv-params
* textwrap
* warnings

Packages needed for Tandem Performance calculator
--------------
* os
* sys
* tkinter
* colorama
* ipywidgets 
* matplotlib
* numpy
* pandas
* traitlets
* cycler
* IPython
* scipy
* openpyxl
* warnings





Instructions for Spectral Efficiency calculator:
------------------

### Running Example Calculations
 
* To begin, open Jupyter Lab (enter "jupyter lab" into terminal or equivalent)

* Navigate to the folder "SE-and-Tandems/notebooks" and open SE\_Calculation\_GUI.ipynb

* Run each cell sequentially until the SE GUI appears (press shift+return while a cell is selected to run)

* To use the GUI, you will need a I-V file and an EQE file for each of the cells you intend to compare as well as a single spectrum file (i.e. AM1.5G). For this tutorial, example files are located in folders under "SE-and-Tandems/Examples/SE"

* Upload the I-V files together by clicking on the "Select Files" icon

* Repeat for the EQE files

* Upload a spectrum file in the "Select Spectral Irradiance File" section

* Click "Calculate." The graph of the input I-V and EQE and the output spectral efficiency will appear. The traditional single junction efficiency for the provided cells will also be displayed 

* Upon calculation, the calculated SE data saves in .csv format to the 'calculatedSE' folder and the figure will be saved in .png format to 'SE\_images.' The saved data files are named "\<input file header>\_calculatedSE\_\<index\>", where 'index' is incremented to avoid overwriting duplicate file names. The SE image files are named "\<input file header>\all_cells_calculatedSE\_\<index\>".

* The calculated SE (calculatedSE) data (.csv file) includes 8 data sets in columns:

    1) QE\_smooth (The provided EQE data will be interpolated with the same wavelength spacing of incident AM1.5 spectrum)
    
    2) wl (Wavelength in nm scale)
    
    3) I\_sun (AM1.5 irradiation)
    
    4) Jsc\_wl (spectrally-resolved Jsc: primarily QE smooth data multiplied by AM1.5 irradiation)
    
    5) v (Voltage column from the studied cell's I-V)
    
    6) i (Current column from the studied cell's I-V) 
    
    7) Eg (computed bandgap using EQE derivation in respect to wavelength) 
    
    8) SE (calculated Spectral Efficiency data)

* To run a new calculation, simply replace the necessary files by clicking the respective "Files Selected" in the GUI. You can also click "Clear Output Results" to remove calculations from the notebook output before running a new calculation, but this is not necessary. You do not need to rerun the python cells. 

### File format and data preparation

* The names of both I-V and EQE files for the same single cell must be identical, with the exception of the suffix

* As an example, I-V file Si\_H2M2W3\_IV corresponds to EQE file Si\_H2M2W3\_QE

* The I-V, EQE and spectrum files must be in .csv, .xls or .xlsx (Excel) format. 

/// For better understanding, please refer to the template I-V and EQE files in the 'Examples/SE' folder

### I-V data 

* The I-V data should be plotted in the graph's first quadrant

* For example, the voltage and current ranges should be positive at the maximum power point

* The voltage and current units should be V and mA, respectively

* The first header of I-V data should be v for voltage column and i for current column

### EQE data

* The wavelength should be in nm scale, and the EQE should be expressed in percentage

* The EQE data's first header should be defined by wavelength (nm) for the wavelength column and QE (%) for the EQE column

### Spectrum file

* If you use your own spectrum file, please prepare it in the same way as the provided AM1.5G file



Instructions for Tandem Performance calculator:
----------------------------------

* Open Jupyter Lab, navigate to "SE-and-Tandems/notebooks" and open Tandem\_Performance\_GUI.ipynb

* You will need the results (.csv) of a SE calculation (described above) or you can use the example input files in "SE-and-Tandems/Examples/Tandem Performance"

* Run the notebook cells sequentially until the tandems GUI appears

* Upload the top cell and bottom cell SE files using the respective file selection buttons. Several top cell files can be uploaded simultaneously with a single bottom cell file

* Next, choose a coupling method from the dropdown menu. Coupling can be 

    1) ideal: no files needed (100% of light below the top cell's bandgap will be transmitted),
    
    2) fixed: provide the percent transmission with the sliding bar, no additional files needed (you can specify how much light will be transmitted below the top cell's bandgap) or
    
    3) lambda coupling through a separate transmission file for the topmost cell by following "Select Transmission Files for lambda Coupling." 
    
  Note on bandgaps: For ideal and fixed coupling, check the band gap of the top cell SE file(s) to ensure Eg is physical and not derived from calculation noise. There is no need to include a bandgap for lambda coupling. 
  
  Note on lambda coupling: For accuracy when using lambda coupling, you will need to provide a transmission file for each top cell. For the script to match the proper transmission file to each top cell spectral efficiency, name the transmission file the same as the corresponding top cell file minus '\_calculatedSE'. To differentiate the files, you may optionally append '\_transmission', '\_trans', '\_tran', or '\_tr' somewhere in the transmission file filename. As an example, if the top cell spectral efficiency file is named 'top\_cell\_15\_calculatedSE\_6', two appropriate transmission file names would be 'top\_cell\_15\_transmission\_6' or 'top\_cell\_15\_6_tr'. An example transmission file can be found in "Examples/Tandem Performance."

* Once the proper files are selected, click "Calculate." The results will display below the GUI and save data to "tandems_data" and "tandems\_images" as a .csv and .png, respectively. You can click "Clear Output Results" to remove calculations from the notebook output before running a new calculation, but this is not necessary. 

* Save files are named: 'tandem data\_'/'tandem graph\_' + top cell name + bottom cell name + coupling method + transmitted light + _index (where 'index' is incremented to avoid overwriting duplicate file names)


### 2T tandem calculation based on top and bottom cell spectral efficiency 

* "Top Cell" max power point = Top Cell SE * I<sub>sun</sub> * ((min[J<sub>sc,top</sub>,J<sub>sc,bottom,filtered</sub>])/J<sub>sc,top</sub>)
* "Bottom Cell" max power point = Bottom Cell SE * I<sub>sun</sub> * ((min[J<sub>sc,top</sub>,J<sub>sc,bottom,filtered</sub>])/J<sub>sc,bottom,filtered</sub>)
*  SE is Spectral Efficiency 
* The "I<sub>sun</sub>" is AM1.5 spectral irradiance [W/m^2]
* Here, the "min" function selects the lower J<sub>sc</sub> from the Top and filtered Bottom cells
* 2T Tandem efficiency = "Top Cell" max power point + "Bottom Cell" max power point


### Accurate form of 2T tandem calculation by generating tandem IV from top and bottom cell IV

* Obtaining the tandem IV by adding the voltage of the top cell to the voltage of the filtered bottom cell at the same current of either the top or filtered bottom cell that has a lower J<sub>sc</sub>
* After that, taking the max power point from the obtained tandem IV curve, 2T Tandem efficiency = max power point of "tandem IV"



## Known Errors 

* Python windows for file selection open behind the browser window, not in front (and may raise a CATransaction warning, but should not influence use of the calculators otherwise).

* Separate python windows open for each notebook GUI once the first file is selected. These windows stay open even after calculation until the terminal session is closed. Once you are finished using the GUI, you can force quit the idle python instances.

* Cancelling a file dialog window clears any previous input.
