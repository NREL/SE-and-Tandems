#SE calculator with GUI interface

import os
import sys
from tkinter import Tk, filedialog, Toplevel

import colorama
from colorama import Fore

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import traitlets
from cycler import cycler
from IPython.display import clear_output, display, Markdown
from scipy.interpolate import interp1d
import textwrap

from iv_params.iv_params import IV_Params


#Defining the class for sending different error message
class WrongFileFormat(Exception):
    pass

class WrongIVandQEName(Exception):
    pass

class WrongDataPreparation(Exception):
    pass

class SpectrumDataMissing(Exception):
    pass


#get the constant factor of spectralJsc   
def SEJsc_factor():
    q      = 1.602176634*10**-19 #(C) # Elementary charge in SI unit
    c      = 2.998*10**8  #m/s
    h      = 6.626*10**-34 #J.s
    k = (q/(h*c*10**9))   #10**9 for m to nm conversion
    return k


# Clip the filename to get cell name and material name
def clip_name(filename):
    filenopath = os.path.basename(filename)
    tpl1 = filenopath.partition('_')
    cell = tpl1[0]
    tpl2 = filenopath.rpartition('_')
    name = tpl2[0]
    return cell, name


#input I-V, EQE and spectra files and then calculate spectral efficiency by processing other parameters: Voc, FF using iv_params and Jsc(λ) by multplying EQE(λ) with spectra Irr(λ)

def calc_se_eff(irrfile, IV_file_list, QE_file_list): #irrfile=spectra file, IV_file_list=I-V file, QE_file_list= EQE file
    

    # Create list for SE
    SE_list = []
    #Create list for single-junction efficiency
    SJEffi = []
    #Create list for filenames
    FN_list = []
    #Create list for bandgap
    Eg_list = []
    
    # input constant factor of spectralJsc from previously defined function: def SEJsc_factor() 
    spectralJscfactor = SEJsc_factor() 
    
    
    
    #set up plots for input I-V and EQE data
    fig1, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 4))
    ax1.set_ylabel("Current (mA)",color="red",fontsize=14) #was Current Density (mA cm$^{-2}$)
    ax1.set_xlabel("Voltage (V)",fontsize=14)
    ax2.set_ylabel("EQE (%)",color="blue",fontsize=14)
    ax2.set_xlabel("Wavelength (nm)",fontsize=14)
    fig1.suptitle('Input I-V and EQE data', fontsize=16)
    
    #set up plots for spectral efficiency data
    fig2, ax3 = plt.subplots(1, figsize=(7, 5))   #(6,4) sets height same as fig1, but slightly wider than one subplot
    ax3.set_xlim(280,1300)
    ax3.set_xlabel("Wavelength (nm)",fontsize=14)
    ax3.set_ylabel("Spectral Efficiency (%)",fontsize=14)
    fig2.suptitle('Spectral Efficiency', fontsize=16)

    #set up plots for individual spectral efficiency data plots to save
    fig3, ax4 = plt.subplots(1, figsize=(7, 5))   #(6,4) sets height same as fig1, but slightly wider than one subplot
    ax4.set_xlim(280,1300)
    ax4.set_xlabel("Wavelength (nm)",fontsize=14)
    ax4.set_ylabel("Spectral Efficiency (%)",fontsize=14)
    fig3.suptitle('Spectral Efficiency', fontsize=16)
        
    #reading spectra file
    for filename3 in irrfile:
            
        if filename3[-3:] == 'csv':
            Spectrafile = pd.read_csv(filename3)
            
        elif filename3[-4:] == 'xlsx':
            Spectrafile  = pd.read_excel(filename3)
            
        elif filename3[-3:] == 'xls':
            Spectrafile  = pd.read_excel(filename3)
                
    df = Spectrafile
    df = df.rename(columns={"Wavelength (nm)": "wl", "Spectral irradiance (W m-2 nm-1)": "I_sun","Cumulative photon flux (cm-2 s-1)":"cum_flux"})
    compiledSE = pd.DataFrame() 
    compiledSE['wl'] = df['wl']
    compiledSE['I_sun'] = df['I_sun']
    
    
    #reading inputted I-V and EQE files
    for filename1,filename2 in zip(IV_file_list, QE_file_list):
        
        if filename1[-3:] == 'csv':
            IV = pd.read_csv(filename1)
            
        elif filename1[-4:] == 'xlsx':
            IV = pd.read_excel(filename1)
            
        elif filename1[-3:] == 'xls':
            IV = pd.read_excel(filename1)
            
        else:
            raise WrongFileFormat
        
        # Set up names for I-V.
        IV_cell, IV_name = clip_name(filename1)
        
      
        if filename2[-3:] == 'csv':
            QE = pd.read_csv(filename2)
            
        elif filename2[-4:] == 'xlsx':
            QE = pd.read_excel(filename2)
            
        elif filename2[-3:] == 'xls':
            QE = pd.read_excel(filename2)
            
        else:
            raise WrongFileFormat
        
        # Set up names for EQE.
        QE_cell, QE_name = clip_name(filename2)
        
        
        
        # Calculate only when file names are same for I-V and EQE except suffix 
        if IV_name==QE_name:
            
            try:
            
                #remove junk from end of QE data
                if 'end' in QE['wavelength (nm)'].unique():
                    end_row = QE.index[QE['wavelength (nm)'] == 'end'][0]
                    QE = QE.iloc[:end_row].copy()
                QE = QE.rename(columns={"wavelength (nm)": "wl_raw", "QE (%)": "QE"})
                
            
                #redefining the columns name 
                QE["wl_raw"]= QE["wl_raw"]
                QE["QE"]= QE["QE"]
            
                #doing EQE derivation with respect to wavelength
                dydx = np.gradient(QE["QE"], QE["wl_raw"])
                
                #finding the minimum value of derivative EQE    
                min_value=min(dydx)
                #finding the index of minimum value of derivative EQE    
                min_index=list(dydx).index(min_value)
                #get the wavelength for the index of minimum value of derivative EQE    
                wavelength_min_deri=QE["wl_raw"][min_index]
                #convert the wavelength to eV
                Bandgap=1.24/(wavelength_min_deri*0.001)
                Bandgap= float("{:.2f}".format(Bandgap))

                #passing the bandgap information as a message
                msg3=str(QE_name) + ": " + '\033[1m' + str(Bandgap) + '\033[0m' + " eV"
                Eg_list.append(msg3)

                #Smooth QE data by interpolating 
                newQE = interp1d(QE['wl_raw'],QE['QE'], bounds_error = False )
                test = newQE(df['wl'])
                SE = pd.DataFrame(data=test, columns =['QE_smooth']) 
                SE['wl'] = df['wl']
                SE['I_sun'] = df['I_sun']
                
                #getting the values of voltage (Voc) and current (Jsc) from I-V curve using iv_params     
                ivFOM = IV_Params(IV.v, IV.i)
                ivFOM.calc_iv_params()
                ivFOM =ivFOM.calc_iv_params()
            
                #calculate Jsc(lambda).....spectrally-resolved Jsc
                SE['Jsc_wl'] = spectralJscfactor*SE['wl']*SE['QE_smooth']* SE['I_sun']
                SE['Jsc_wl'] = SE['Jsc_wl'].replace(np.nan,0)

                #Integrated Jsc from Spectrally-resolved Jsc as well as normalized from mA to A by dividing the electric current value by 1000
                SE_Integrated_Jsc= np.trapz(SE['Jsc_wl'])/1000
        
                #Adjusting the integrated Jsc from EQE to the Jsc from I-V - taking the Jsc from I-V curve into consideration as a benchmark
                Factor_Jsc=(ivFOM['isc'] - SE_Integrated_Jsc) / SE_Integrated_Jsc
                Factor_Jsc=1+Factor_Jsc
                SE['Jsc_wl']=SE['Jsc_wl']*Factor_Jsc
                
                #exporting data in columns of SE_data file 
                SE['QE_smooth']=SE['QE_smooth']*Factor_Jsc
                SE['v']=IV.v
                SE['i']=IV.i
                SE['Eg']=Bandgap
                
                #input power
                Pin=1  ## considering the standard Am1.5G input power for efficiency calculations (1 kW/m^2)
            
                #Efficiency using Integrated spectally-resolved Jsc from provided EQE curve
                Effi_SE = (ivFOM['voc'] * SE_Integrated_Jsc * Factor_Jsc * ivFOM['ff']) / Pin
                Effi_SE = float("{:.2f}".format(Effi_SE))
            
                #Efficiency using provided I-V curve
                Effi_IV = ivFOM['pmp'] / Pin   ## pmp is maximum power point of a given I-V curve 
                Effi_IV = float("{:.2f}".format(Effi_IV))  ## taking values up to two decimal 
            
                #passing message for the calculated efficiency from I-V curve
                msg= str(IV_name) +": " + '\033[1m' + str(Effi_IV) + '\033[0m' + "%"
            
                #Single-junction efficiency data in list
                SJEffi.append(msg)

                # Calculate SE efficiency
                SE['SE'] = ivFOM['voc'] * ivFOM['ff'] * SE['Jsc_wl'] / (SE['I_sun'].replace(0,np.nan))
                
            except:
                raise WrongDataPreparation
            

            exclude = ['QE_smooth', 'wl', 'I_sun', 'Jsc_wl', 'v', 'i', 'Eg']
            
            for col in SE.columns:
                if col not in exclude:
                    ax1.plot(IV['v'], IV['i'], label = f"{IV_name}")
                    ax2.plot(SE['wl'], SE['QE_smooth'], markerfacecolor='none', label= f"{QE_name}")
                    ax3.plot(SE['wl'], SE['SE'], label= f"{IV_name}")
                    ax4.plot(SE['wl'], SE['SE'], label= f"{IV_name}")

            #save SE data including all other processed data: Jsc(lamda), Smooth QE
            SE_line = [SE]
            SE_list.append(SE_line)
           
            
            #Saving data in folder
            i = 1
            while os.path.exists('../SE_data/{}{:d}.csv'.format(IV_name + '_calculatedSE' + '_', i)):
                i += 1
            SE.to_csv('../SE_data/{}{:d}.csv'.format(IV_name + '_calculatedSE' + '_', i))
            
            #saving file name in list
            msg2 = str(IV_name) + '_calculatedSE' + '_' + str(i)
            FN_list.append(msg2)

            ax1.legend(loc=1)
            ax2.legend(loc=0)
            ax3.legend(loc=1)
            
            #Save Figure
            plt.savefig('../SE_images/{}{:d}.png'.format('all_cells_calculatedSE' + '_', i, dpi=300))
            ax4.legend(loc=1)
        else:
            raise WrongIVandQEName

    ax1.set_ylim(0,)
    ax1.set_xlim(0,)
    ax3.set_ylim(0,)
    ax4.set_ylim(0,)

    #Saving cumulative SE image in folder
    i = 1
    while os.path.exists('../SE_images/{}{:d}.png'.format('all_cells_calculatedSE' + '_', i)):
        i += 1
    plt.savefig('../SE_images/{}{:d}.png'.format('all_cells_calculatedSE' + '_', i, dpi=300))
        
    
#commenting next return line prevents printing of figures to GUI, but not to log
    return SE, fig1, fig2, msg, FN_list, SJEffi, Eg_list


#GUI-part
class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""
    files = traitlets.List([], help="List of file paths").tag(sync=True)
    _files = traitlets.List([], help="List of file paths").tag(sync=True)

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(SelectFilesButton, self).__init__(*args, **kwargs)
        # Create the button.
        self.icon = "folder-open"
        self.style.button_color = "orange"
        # Set on click behavior.
        self.on_click(self.select_files)

        

    @staticmethod
    def select_files(b):
        """Generate instance of tkinter.filedialog.
        Parameters
        ----------
        b : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        b._files = filedialog.askopenfilename(multiple=True)
        b._files.sort()
        b.files = [os.path.basename(file) for file in b._files]
        b.description = "Files Selected"
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"
        
def se_calculation():
    
    def clear_plots(b):
        with outbox:
            clear_output()
    
    def on_button_click(b):
        #I-V files
        IV_file_list = topselect._files
        #EQE files
        QE_file_list = botselect._files
        #spectra file
        irrfile = spectra._files
        
        #recalling function calc_se_eff(IV_file_list, QE_file_list) where all the calculations done
        try:
            SE, fig1, fig2, msg, FN_list, SJEffi, Eg_list = calc_se_eff(irrfile, IV_file_list, QE_file_list) #tried semicolon here, does nothing
            
            with outbox: 
                clear_output() 
                print("-----------------------------------------")
                print(Fore.RED +'Traditional single-junction efficiency')
                print(Fore.BLACK +"-----------------------------------------")
                print(Fore.BLACK +'\n'.join(SJEffi))
                print(Fore.BLACK +"-----------------------------------------")
                print(Fore.RED +"Bandgap (Eg) of sub-cell absorber materials")
                print(Fore.BLACK +"-----------------------------------------")
                print(Fore.BLACK +'\n'.join(Eg_list))
                print(Fore.BLACK +"-----------------------------------------")
                print(Fore.RED +'Data and figures are saved in folders\n "SE_data" and "SE_images," respectively.')
                print(' ')
                print(Fore.RED +'The output filenames are: ')
                print(Fore.BLACK +',\n'.join(FN_list))
                print("-----------------------------------------")
                print(' ')
                display(fig1)
                display(fig2)   
    
                
        #Sending error message       
        except WrongFileFormat:
            with outbox:
                clear_output()
                print(Fore.RED +'Files format are not correct,\n the files should be csv or Excel format')
                
        except WrongIVandQEName:
            with outbox:
                clear_output()
                print(Fore.RED +'The name of files: IV and QE are not matched,\n please see the instructions')
                
        except SpectrumDataMissing:
            with outbox:
                clear_output()
                print(Fore.RED +'AM1.5G Spectrum data missing, this file must be located in the folder of "SE_inputdata"')
                
        except WrongDataPreparation:
            with outbox:
                clear_output()
                print(Fore.RED +'There is an error in data preparation; please double-check the instructions! Tips: Shift I-V files to Quadrant I (positive Voc and Jsc). Check column headers for typos. Eliminate any metadata aside from column names and values if needed.)')
                    
    '''
    use interactive_input for GUI in IPython
    '''
        
    cell_layout = widgets.Layout(display='inline_flex',
        flex_flow='row',
        justify_content='flex-end',
        width='320px')
        #width='100%')  
        
    #Input and Output Buttons
    
    # Column 1
    s1label = widgets.Label(value='Step 1: Select File(s) ')
    toplabel = widgets.Label(value='Select I-V file(s)')

    # Button for selecting IV file (s)
    topselect = SelectFilesButton(description = 'Select Files', layout=cell_layout)
    IV_file_list = topselect.files
    toptext = widgets.Select(options = [], layout=cell_layout)
    toplink = widgets.link((topselect, 'files'), (toptext, 'options'))                    
    botlabel = widgets.Label(value='Select EQE file(s)')

    # Button for selecting QE file(s)
    botselect = SelectFilesButton(description = 'Select Files', layout=cell_layout)

    QE_file_list = botselect.files

    bottext = widgets.Select(options = [], layout=cell_layout)
    botlink = widgets.link((botselect, 'files'), (bottext, 'options'))
    
    # Column 2
    s2label = widgets.Label(value='Step 2: Select Spectral Irradiance File')
    #Button for importing spectra file
    spectra = SelectFilesButton(
        description="Upload Spectral Irradiance File", layout=cell_layout
    )
    irrfile = spectra.files

    spectratext = widgets.Select(options = [], layout=cell_layout)
    spectralink = widgets.link((spectra, 'files'), (spectratext, 'options'))
    
    # Column 3
    #Buttons for calculating results and clearing the results
    s3label = widgets.Label(value='Step 3: Calculate SE and Display Results')
    outbox = widgets.Output()

    
    # Actions       
    
    # button to start calculations
    calcbutton = widgets.Button(
        description="Calculate", 
        layout=cell_layout,
        button_color='lightgreen')
    
    clr_plots = widgets.Button(
            description="Clear Output Results",
            layout=cell_layout)

    leftcntrls = [s1label, toplabel, topselect, toptext, botlabel, botselect, bottext]
    
    ctrcntrls = [s2label, spectra, spectratext]
    
    rtcntrls = [s3label, calcbutton, clr_plots]

    dncntrls = [outbox]
    
    # call function to start calculations 
    calcbutton.on_click(on_button_click)
    clr_plots.on_click(clear_plots)
    
    # user interface        
    box_layout = widgets.Layout(display='flex',
        flex_flow='column',
        align_items='center',
        width='33.33%',                        
        height = '380px')

    #for output results
    box_layout_2 = widgets.Layout(display='flex',
        flex_flow='column',
        align_items='center',
        width='100%',
        height = '1400px')
    
    leftbox = widgets.VBox(leftcntrls, layout=box_layout)
    
    ctrbox = widgets.VBox(ctrcntrls, layout=box_layout)

    rtbox = widgets.VBox(rtcntrls, layout=box_layout)

    dnbox = widgets.VBox(dncntrls, layout=box_layout_2)

    ui = widgets.HBox([leftbox, ctrbox, rtbox])
    ui2 = widgets.HBox([dnbox])
        
    display(ui, ui2)

if __name__=="__main__":
    se_calculation()
