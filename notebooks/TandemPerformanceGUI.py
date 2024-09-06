#Tandem Performance calculator with GUI interface

import os
import sys 
from tkinter import Tk, filedialog

import ipywidgets as widgets
from IPython.display import clear_output, display
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import traitlets
from cycler import cycler
from scipy.interpolate  import interp1d

import colorama
from colorama import Fore

from iv_params.iv_params import IV_Params

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["#1F77B4", "#AEC7E8", "#FF7F0E","#FFBB78", "#2CA02C", "#98DF8A","#D62728", "#FF9896", "#8C564B","#C49C94", "#E377C2", "#F7B6D2","#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5", "#9467BD", "#C5B0D5"])


#Defining the class for sending different error message
class NoSelectionMade(Exception):
    pass

class WrongFileFormat(Exception):
    pass

class WrongLambdaFileFormat(Exception):
    pass

class FileMissing(Exception):
    pass

class BandgapMissing(Exception):
    pass

# Clip the filename to get cell name
def clip_name(filename):
    filenopath = os.path.basename(filename)
    tpl1 = filenopath.rpartition('.')
    cell_name = tpl1[0].replace('_calculatedSE', '')
    return cell_name

#Doing all 4T calculations here   
def calc_tandem_4T_eff(top_file_list, bot_file_list, lambdafile, coupling, transmission): #top_file_list= SE file of top cells, bot_file_list=SE file of bottom cells, lambdafile= Transmission file of top cell, coupling= different optical coupling methods such as ideal, fixed and lambda, and transmission= controlling light transmission using fixed coupling onto bottom cells
    
    
    # Create list for tandems
    tandem_list = []
    # Create list for QE
    qe_list = []
    # Create list for 4T tandem efficiency data using SE 
    tandem_4T = []
    # Create list for 4T tandem efficiency data using J-V curves
    tandem_iv = []
    
    #plotting tandem power density
    fig, (ax1) = plt.subplots(1, figsize=(6, 4))
    ax1.set_ylabel("Power Density (mW cm$^{-2}$ nm$^{-1}$)",color="red",fontsize=14)
    ax1.set_xlabel("Wavelength (nm)",fontsize=14)
    fig.suptitle('Tandem power density', fontsize=16)
    if coupling == 'fixed':
        fig.suptitle('4T Tandem power density \n'+coupling+' coupling, '+str(transmission)+'% transmission', fontsize=16)
    else:
        fig.suptitle('4T Tandem power density \n'+coupling+' coupling', fontsize=16)
    fig.tight_layout()
    
    #plotting QE
    fig2, (ax2) = plt.subplots(1, figsize=(6, 4))
    ax2.set_ylabel("QE (%)",color="red",fontsize=14)
    ax2.set_xlabel("Wavelength (nm)",fontsize=14)
    fig2.suptitle('EQE', fontsize=16)
    fig2.tight_layout()
    
    #plotting SE 4T
    fig3, (ax3) = plt.subplots(1, figsize=(6, 4))
    ax3.set_ylabel("SE (%)",color="red",fontsize=14)
    ax3.set_xlabel("Wavelength (nm)",fontsize=14)
    if coupling == 'fixed':
        fig3.suptitle('Spectral Efficiency of 4T Tandem \n'+coupling+' coupling, '+str(transmission)+'% transmission', fontsize=16)
    else:
        fig3.suptitle('Spectral Efficiency of 4T Tandem \n'+coupling+' coupling', fontsize=16)
    fig3.tight_layout()

      
    #plotting SE SJ
    fig4, (ax4) = plt.subplots(1, figsize=(6, 4))
    ax4.set_ylabel("SE (%)",color="red",fontsize=14)
    ax4.set_xlabel("Wavelength (nm)",fontsize=14)
    fig4.suptitle('Spectral Efficiency of Single-junction', fontsize=16)
    fig4.tight_layout()
    
    #4T tandem using I-V curves
    fig5, (ax5) = plt.subplots(1, figsize=(6, 4))
    ax5.set_ylabel("Current (mA/cm^2)",color="red",fontsize=14)
    ax5.set_xlabel("Voltage (V)",fontsize=14)
    fig5.suptitle('Cells I-V', fontsize=16)
    fig5.tight_layout()
    

    #define list of cell names for loop to prevent double counting of cells
    cellNameList = []

    ## Loop through bottom cells first, and then top cells.
    for filename1 in bot_file_list:

        if filename1[-3:] == 'csv':
            bot_se = pd.read_csv(filename1)
            
        elif filename1[-4:] == 'xlsx':
            bot_se = pd.read_excel(filename1)
            
        elif filename1[-3:] == 'xls':
            bot_se = pd.read_excel(filename1)
            
        else:
            raise WrongFileFormat
    
        
        ## Set up names for bottom cell.
        bot_name = clip_name(filename1)

        for filename2 in top_file_list:
            
            if filename2[-3:] == 'csv':
                top_se = pd.read_csv(filename2)
            
            elif filename2[-4:] == 'xlsx':
                top_se = pd.read_excel(filename2)
            
            elif filename2[-3:] == 'xls':
                top_se = pd.read_excel(filename2)
            
            else:
                raise WrongFileFormat
    
            
            # Set up names for top cell.
            top_name = clip_name(filename2)
            
            # Create f columns in each dataframe to define the ideal and fixed coupling
            
            #top cell flag (1 above Eg (top) and 0 below Eg(top))
            top_se['ftop'] = np.NaN
            #bottom cell flag holds transmission value based on coupling method
            bot_se['fbot'] = np.NaN

            # Define coupling
            
            # conditions for lambda coupling 
            if coupling == 'lambda':
                #Using topcell transmission file
                try:
                    top_se['ftop'] = 1
                    
                    for filename3 in lambdafile:

                        transmission_name = clip_name(filename3).replace('_transmission','').replace('_trans','').replace('_tran','').replace('_tr','')
                        if top_name == transmission_name:

                            if filename3[-3:] == 'csv':
                                Transm = pd.read_csv(filename3)
            
                            elif filename3[-4:] == 'xlsx':
                                Transm = pd.read_excel(filename3)
            
                            elif filename3[-3:] == 'xls':
                                Transm = pd.read_excel(filename3)
                            else:
                                raise WrongLambdaFileFormat
                        
                    smoothTr = interp1d(Transm['wavelength (nm)'],Transm['transmission'], bounds_error = False ) 
                    topTr = pd.DataFrame(data=smoothTr(bot_se['wl']), columns = ['Tr_smooth'])
                    bot_se['fbot'] = topTr['Tr_smooth']
                    bot_se['fbot'] = bot_se['fbot'].replace(np.nan, 0)


                except:
                    raise FileMissing
                    
                
            # conditions for ideal coupling 
            elif coupling == 'ideal':
                transmission = 100
                # Extracting bandgap for ideal and fixed coupling from top cell calculatedSE file
                try:
                    bandgap = top_se['Eg'][0]
                    top_bg_nm = round (1.24/bandgap * 1000)  
                except:
                    raise BandgapMissing
    
                # Define ftop for the top cell.
                #top_se['wl']=photon wavelength (nm)
                top_se.loc[top_se['wl'] < top_bg_nm, 'ftop'] = 1
                top_se.loc[top_se['wl'] >= top_bg_nm, 'ftop'] = 0
          
                # Assuming perfect transmission.
                bot_se.loc[bot_se['wl'] < top_bg_nm, 'fbot'] = 0
                bot_se.loc[bot_se['wl'] >= top_bg_nm, 'fbot'] = 1
                
            # conditions for fixed coupling 
            elif coupling == 'fixed':
            
                # Extracting bandgap for ideal and fixed coupling from top cell calculatedSE file
                try:
                    bandgap = top_se['Eg'][0]
                    top_bg_nm = round (1.24/bandgap * 1000)  
                except:
                    raise BandgapMissing
    
                # Define ftop for the top cell.
                top_se.loc[top_se['wl'] < top_bg_nm, 'ftop'] = 1
                top_se.loc[top_se['wl'] >= top_bg_nm, 'ftop'] = 0
                
                # Assuming a fixed percentage of transmission.
                bot_se.loc[bot_se['wl'] < top_bg_nm, 'fbot'] = 0
                bot_se.loc[bot_se['wl'] >= top_bg_nm, 'fbot'] = transmission/100

  
            else:
                raise NoSelectionMade
        
            # FUNCTION OPPORTUNITY HERE? Need to numerically integrate each subcell.
            
            
            ##### 4T calculations 
            
            # Top cell
            #computing topcell power
            
            #from the provided top cell file - top_se['SE'] column for spectral efficiency (percent) and top_se['I_sun'] column for solar irradiance (W/m^2)  
            #top_se['ftop'] is top cell flag (1 above Eg (top) and 0 below Eg(top))
            top_se['effic'] = top_se['SE']*top_se['ftop']*top_se['I_sun']
            #top_se['effic'] = power/nm (W/nm*m"2)
            top_se['effic'] = top_se['effic'].replace(np.nan,0)
            
            
            #computing topcell efficiency = computing topcell power/standard incident specturm power (1000 W/m^2)
            top_eff = np.trapz(top_se['effic'])/1000
            #result until two decimal points
            top_eff = float("{:.2f}".format(top_eff))
            
    
    
            # Bottom cell
            #computing bottom cell power
            
             #from the provided top cell file - bot_se['SE'] column for speactral efficiency (percent) and bot_se['I_sun'] column for solar irradiance (W/m^2)  
            #bot_se['fbot'] is top cell flag (0 above Eg (top) and 1 below Eg(top))
            bot_se['effic'] = bot_se['SE']*bot_se['fbot']*bot_se['I_sun']
            
            #bot_se['effic'] = power/nm (W/nm*m"2)
            bot_se['effic'] = bot_se['effic'].replace(np.nan,0)
            
            #computing bottom cell efficiency: bottom cell power/standard incident spectrum power (1000 W/m^2)
            bot_eff = np.trapz(bot_se['effic'])/1000
            bot_eff = float("{:.2f}".format(bot_eff))
            
            # Calculate tandem efficiency using SE
            tandem_eff = top_eff + bot_eff
            tandem_eff = "{:.2f}".format(tandem_eff)
            
            #constant factor of spectrally_Jsc
            q = 1.602176634*10**-19  # Elementary charge in SI unit
            c = 2.998*10**8
            h = 6.626*10**-34
            spectralJscfactor = (q/(h*c*10**9))   # conversion from m to nm
            
            # top QE
            top_se['QE'] = top_se['QE_smooth']*top_se['ftop']
            top_se['QE'] = top_se['QE'].replace(np.nan,0)
            
            #calculate Jsc from top QE
            EQE_1 = spectralJscfactor*top_se['wl']*top_se['QE']* top_se['I_sun']
            EQE_1= EQE_1.replace(np.nan,0)
            
            #Integrated Jsc from Spectally-resolved Jsc (top cell).....normalized from mA to A by dividing the electric current value by 1000
            QE_Integrated_top_Jsc= np.trapz(EQE_1)/1000
            QE_Integrated_top_Jsc = float("{:.2f}".format(QE_Integrated_top_Jsc))
            
            # bottom cell QE
            bot_se['QE']=bot_se['QE_smooth']*bot_se['fbot']
            bot_se['QE'] = bot_se['QE'].replace(np.nan,0)
            
            #calculate Jsc from bottom cell QE
            EQE_2 = spectralJscfactor*bot_se['wl']*bot_se['QE']*bot_se['I_sun']
            EQE_2= EQE_2.replace(np.nan,0)
            
            #Integrated Jsc from Spectally-resolved Jsc (bottom cell).....normalized from mA to A by dividing the electric current value by 1000
            QE_Integrated_bot_Jsc= np.trapz(EQE_2)/1000
            QE_Integrated_bot_Jsc = float("{:.2f}".format(QE_Integrated_bot_Jsc))
           
            #passing the information of calculated Jsc from QE
            msg2="Top cell Jsc from " + str(top_name) + "'s QE" + " is " + '\033[1m' + str(QE_Integrated_top_Jsc) + '\033[0m' + " mA/cm^2" + "\nBottom cell Jsc from filtered " + str(bot_name) + "'s QE" + " is " + '\033[1m' + str(QE_Integrated_bot_Jsc) + '\033[0m' + " mA/cm^2"
            
            #passing QE information
            qe_list.append(msg2)
            
            
            #plotting QE for 4T tandem
            if top_name not in cellNameList:
                ax2.plot(top_se['wl'], top_se['QE'], label="top QE " + top_name) 
            if coupling == 'lambda':
                ax2.plot(bot_se['wl'], bot_se['QE'], label=bot_name + " filtered through " + top_name)
            elif bot_name not in cellNameList:
                ax2.plot(bot_se['wl'], bot_se['QE'], label=bot_name + ", transmission = " + str(transmission) + "%")

            #spectral efficiency for single-junction 
            top_se['1JSE']= top_se['SE']
            bot_se['1JSE']= bot_se['SE']
            #spectral efficiency for 4T tandem
            top_se['4TSE']= top_se['SE']*top_se['ftop']
            bot_se['4TSE']= bot_se['SE']*bot_se['fbot']  
            
            
            ##calculating max power point and Jsc from J-V curve using iv_param
            
            #top cell
            #reading i and v from calculatedSE files
            ivFOM_top = IV_Params(top_se.v, top_se.i)
            ivFOM_top.calc_iv_params()
            ivFOM_top =ivFOM_top.calc_iv_params()
            #max power point
            Topcell_Max_Power_Point=ivFOM_top['pmp']
            #Jsc from J-V curve
            Jsc_top=ivFOM_top['isc']
           
            
            #bottom cell
            #reading i and v from calculatedSE files
            ivFOM_bot_raw = IV_Params(bot_se.v, bot_se.i)
            ivFOM_bot_raw.calc_iv_params()
            ivFOM_bot_raw =ivFOM_bot_raw.calc_iv_params()
            # original Jsc from J-V curve
            Jsc_bot_raw=ivFOM_bot_raw['isc']
            
            ##filtering bottom cell 
            
            #subtracting Jsc (integrated bottom cell QE in 4T tandem) from Jsc (J-V curve of bottom cells)
            delta_Jsc=Jsc_bot_raw - QE_Integrated_bot_Jsc
            #subtracting delta_Jsc from column i (current) of original bottom cell file
            filtered_bot_i=bot_se.i-delta_Jsc
            #form a new I-V for bottom cell
            ivFOM_bot = IV_Params(bot_se.v, filtered_bot_i)
            ivFOM_bot.calc_iv_params()
            ivFOM_bot =ivFOM_bot.calc_iv_params()
            #calculate max power point from new I-V of bottom cell
            Botcell_Max_Power_Point=ivFOM_bot['pmp']
            
            #tandem cell using J-V curves
            pmp_2T=(Topcell_Max_Power_Point+Botcell_Max_Power_Point)*10.0 #10.0 to convert from mW/cm^2 to W/m^2
            pmp_2T=float("{:.2f}".format(pmp_2T))
            
            #plotting 
            if top_name not in cellNameList:
                ax5.plot(top_se.v.replace(0,), top_se.i.replace(0,), label=("top I-V " + top_name))   #doesn't seem to need zero replacement - maybe it does?
            if bot_name not in cellNameList:
                ax5.plot(bot_se.v, bot_se.i, label=("bottom I-V " + bot_name))
            if coupling == 'lambda':
                ax5.plot(bot_se.v, filtered_bot_i, label=(bot_name + " filtered through " + top_name))
            elif bot_name not in cellNameList:
                ax5.plot(bot_se.v, filtered_bot_i, label=(bot_name + ", transmission = " + str(transmission) + "%"))
            
            #passing information of 4T tandem efficiency 
            #msg5 = "Power Density at Max Power Point for 4T Tandem of " + str(top_name) + " and " + "filtered_" + str(bot_name) + ": " + str(pmp_2T) + " W/m^2"
            #tandem_iv.append(msg5)
            
            # MAKE THE FLAGS BOOLEAN SO TEST IS EASY!
            #Yes = 1, No = 0 flags. Set these for your desired output.
            #outputpath = '../tandems' # name of directory where we will save tandem data
            #Print data in this notebook
            print_data = 1
            #Save data to individual files
            save_ind = 1
            #Save all data to a single file
            save_all = 1
            #Print graphs in this notebook
            print_graphs = 1
                
            ##########
            if save_ind == 1:
            #create a compiled parameter dataframe
                    
                if coupling == 'ideal':
                    transmission = 100
                elif coupling == 'lambda':
                    transmission = "light transmission as a function of wavelength" 

                tandem_data = [[top_name, top_eff, coupling, transmission],[bot_name, bot_eff],["4T_tandem", tandem_eff]]
                columns = ['cell name', 'efficiency', 'coupling', 'transmission']    
                tandemSEtemp = pd.DataFrame(tandem_data, columns=columns) 
                tandem_wl_data = [top_se['wl'],top_se['effic'],top_se['4TSE'],top_se['1JSE'],top_se['QE'],bot_se['wl'],bot_se['effic'],bot_se['4TSE'],bot_se['1JSE'],bot_se['QE'],top_se.v,top_se.i,bot_se.v, filtered_bot_i, bot_se.i,]
                wl_columns = ['top cell wavelength (nm)','top cell efficiency','top cell 4T SE', 'top cell 1J SE','top cell QE','bottom cell wavelength (nm)','bottom cell efficiency','bottom cell 4T SE', 'bottom cell 1J SE','bottom cell QE','top cell voltage','top cell current','bottom cell voltage','filtered bottom cell current','bottom cell current']
                wlResolvedData = pd.DataFrame({col: data for col, data in zip(wl_columns, tandem_wl_data)})
                tandemSE = pd.concat([tandemSEtemp,wlResolvedData],axis=1)
                
            
            #Data saving
            if coupling=='fixed':
                i = 1
                while os.path.exists('../tandems_data/{}{:d}.csv'.format('4T_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
                    i += 1
                tandemSE.to_csv('../tandems_data/{}{:d}.csv'.format('4T_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i))
            else:
                i = 1
                while os.path.exists('../tandems_data/{}{:d}.csv'.format('4T_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
                    i += 1
                tandemSE.to_csv('../tandems_data/{}{:d}.csv'.format('4T_' + top_name + '+' + bot_name + '_' + coupling + '_', i))
            
            #Data Message Printing
            if print_data ==1:
                if coupling == 'lambda':
                    msg= "\033[1m" + "4T Tandem" + "\033[0m" + "\nCoupling Method: " + str(coupling) + "\nlight transmitted onto bottom cell: " + str(transmission) + "\nTop Cell is " + str(top_name) + "\n" + str(top_name) + " limiting efficiency = " + str(top_eff) + " %" + "\nBottom Cell is " + str(bot_name) + "\n" + str(bot_name) + " limiting efficiency = " + str(bot_eff) + " %" + "\nThe tandem efficiency for " + str(top_name) + " and " + str(bot_name) + " is "+ str(tandem_eff) + " %" + "\nThe results are saved in the 'tandems_data' folder and the file name is: 4T_tandem data_" + str(top_name) + "+" + str(bot_name) + "+" + str(coupling) + "+" + str(transmission) + "_" + str(i)
                    tandem_4T.append(msg)   
                else:
                    msg= "\033[1m" + "4T Tandem" + "\033[0m" + "\nCoupling Method: " + str(coupling) + "\nlight transmitted onto bottom cell: " + str(transmission) + " %" + "\nTop Cell is " + str(top_name) + "\nBandgap of Top Cell: " + str(bandgap) + "eV" + "\n" + str(top_name) + " limiting efficiency = " + str(top_eff) + " %" + "\nBottom Cell is " + str(bot_name) + "\n" + str(bot_name) + " limiting efficiency = " + str(bot_eff) + " %" + "\nThe tandem efficiency for " + str(top_name) + " and " + str(bot_name) + " is "+ str(tandem_eff) + " %" + "\nThe results are saved in the 'tandems_data' folder and the file name is: 4T_tandem data_" + str(top_name) + "+" + str(bot_name) + "+" + str(coupling) + "+" + str(transmission) + "_" + str(i)
                    tandem_4T.append(msg)
                
            # Data plotting for power density, 4T SE and SJ SE       
            if print_graphs == 1:
                if top_name not in cellNameList:
                    ax1.plot(top_se['wl'], top_se['effic'], label=top_name)
                    ax3.plot(top_se['wl'], top_se['4TSE'], '--', label=top_name)
                    ax4.plot(top_se['wl'], top_se['1JSE'], '--', label=top_name)
                if bot_name not in cellNameList:
                    ax4.plot(bot_se['wl'], bot_se['1JSE'], '--', label=bot_name)
                if coupling == 'lambda':
                    ax1.plot(bot_se['wl'], bot_se['effic'], label=bot_name + " filtered through " + top_name)
                    ax3.plot(bot_se['wl'], bot_se['4TSE'], '--', label=bot_name + " filtered through " + top_name)
                elif bot_name not in cellNameList:
                    ax1.plot(bot_se['wl'], bot_se['effic'], label=bot_name + ", transmission = " + str(transmission) + "%")
                    ax3.plot(bot_se['wl'], bot_se['4TSE'], '--', label=bot_name + ", transmission = " + str(transmission) + "%")                

            #Define which top and bottom cells in list have already been processed in for loops (duplicates are fine)
            cellNameList.append(top_name)
            cellNameList.append(bot_name) 

    #Plotting axes limits and legends provided outside of loop
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax1.set_ylim(0,)

    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 

    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax3.set_ylim(0,)

    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax5.set_xlim(0,)
    ax5.set_ylim(0,)

    #Saving images
    if coupling=='fixed':
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig.savefig('../tandems_images/{}{:d}.png'.format('4T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig2.savefig('../tandems_images/{}{:d}.png'.format('4T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig3.savefig('../tandems_images/{}{:d}.png'.format('4T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig4.savefig('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig5.savefig('../tandems_images/{}{:d}.png'.format('4T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
    else:
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig.savefig('../tandems_images/{}{:d}.png'.format('4T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig2.savefig('../tandems_images/{}{:d}.png'.format('4T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig3.savefig('../tandems_images/{}{:d}.png'.format('4T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig4.savefig('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('4T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig5.savefig('../tandems_images/{}{:d}.png'.format('4T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)


    return tandem_eff, tandem_data, fig, fig2, fig3, fig4, fig5, tandem_list, tandem_4T, qe_list, tandem_iv
        
    

    

#Doing all 2T calculations here
    
def calc_tandem_2T_eff(top_file_list, bot_file_list, lambdafile, coupling, transmission): #top_file_list= SE file of top cells, bot_file_list=SE file of bottom cells, lambdafile= Transmission file of top cell, coupling= different optical coupling methods such as ideal, fixed and lambda, and transmission= controlling light transmission using fixed coupling onto bottom cells
    
    # Create list for tandems
    tandem_list = []
    # Create list for tandems
    qe_list = []
    # Create list for 2T tandem efficiency data using SE 
    tandem_2T = []
    # Create list for 2T tandem I-V 
    tandem_v = []
    
    #plotting tandem power density
    fig, (ax1) = plt.subplots(1, figsize=(6, 4))
    ax1.set_ylabel("Power Density (mW cm$^{-2}$ nm$^{-1}$)",color="red",fontsize=14)
    ax1.set_xlabel("Wavelength (nm)",fontsize=14)
    if coupling == 'fixed':
        fig.suptitle('2T Tandem power density \n'+coupling+' coupling, '+str(transmission)+'% transmission', fontsize=16)
    else:
        fig.suptitle('2T Tandem power density \n'+coupling+' coupling', fontsize=16)
    fig.tight_layout()
    
    #plotting QE
    fig2, (ax2) = plt.subplots(1, figsize=(6, 4))
    ax2.set_ylabel("QE (%)",color="red",fontsize=14)
    ax2.set_xlabel("Wavelength (nm)",fontsize=14)
    fig2.suptitle('EQE', fontsize=16)
    fig2.tight_layout()
    
    #plotting SE
    fig3, (ax3) = plt.subplots(1, figsize=(6, 4))
    ax3.set_ylabel("SE (%)",color="red",fontsize=14)
    ax3.set_xlabel("Wavelength (nm)",fontsize=14)
    if coupling == 'fixed':
        fig3.suptitle('Spectral Efficiency of 2T Tandem \n'+coupling+' coupling, '+str(transmission)+'% transmission', fontsize=16)
    else:
        fig3.suptitle('Spectral Efficiency of 2T Tandem \n'+coupling+' coupling', fontsize=16)
    fig3.tight_layout()
    
    #plotting SE SJ
    fig4, (ax4) = plt.subplots(1, figsize=(6, 4))
    ax4.set_ylabel("SE (%)",color="red",fontsize=14)
    ax4.set_xlabel("Wavelength (nm)",fontsize=14)
    fig4.suptitle('Spectral Efficiency of Single-junction', fontsize=16)
    fig4.tight_layout()
    
    #2T tandem by adding the voltage of top and bottom cell at the same current
    fig5, (ax5) = plt.subplots(1, figsize=(6, 4))
    ax5.set_ylabel("Current (mA/cm^2)",color="red",fontsize=14)
    ax5.set_xlabel("Voltage (V)",fontsize=14)
    fig5.suptitle('Tandem I-V', fontsize=16)
    fig5.tight_layout()
    
    #2T tandem by adding the voltage of top and bottom cell at the same current
    fig6, (ax6) = plt.subplots(1, figsize=(6, 4))
    ax6.set_ylabel("Current (mA/cm^2)",color="red",fontsize=14)
    ax6.set_xlabel("Voltage (V)",fontsize=14)
    fig6.suptitle('Top and Bottom cells raw I-V', fontsize=16)
    fig6.tight_layout()
    
    #list of cell names to prevent double plotting
    cellNameList = []

    ## Loop through bottom cells first, and then top cells.
    for filename1 in bot_file_list:

        if filename1[-3:] == 'csv':
            bot_se = pd.read_csv(filename1)
            
        elif filename1[-4:] == 'xlsx':
            bot_se = pd.read_excel(filename1)
            
        elif filename1[-3:] == 'xls':
            bot_se = pd.read_excel(filename1)
            
        else:
            raise WrongFileFormat
    
        ## Set up names for bottom cell.
        bot_name = clip_name(filename1)

        for filename2 in top_file_list:
            
            if filename2[-3:] == 'csv':
                top_se = pd.read_csv(filename2)
            
            elif filename2[-4:] == 'xlsx':
                top_se = pd.read_excel(filename2)
            
            elif filename2[-3:] == 'xls':
                top_se = pd.read_excel(filename2)
            
            else:
                raise WrongFileFormat
    
            # Set up names for top cell.
            top_name = clip_name(filename2)
 
            # Create f columns in each dataframe.
            
            #top cell flag (1 above Eg (top) and 0 below Eg(top))
            top_se['ftop'] = np.NaN
            #bottom cell flag (0 above Eg (top) and 1 below Eg(top))
            bot_se['fbot'] = np.NaN
            
            
            # Define coupling
        
            if coupling == 'lambda':
                #Using topcell transmission file
                try:
                    top_se['ftop'] = 1
                    
                    for filename3 in lambdafile:
                        transmission_name = clip_name(filename3).replace('_transmission','').replace('_trans','').replace('_tran','').replace('_tr','')
                        if top_name == transmission_name:
                            if filename3[-3:] == 'csv':
                                Transm = pd.read_csv(filename3)
            
                            elif filename3[-4:] == 'xlsx':
                                Transm = pd.read_excel(filename3)
            
                            elif filename3[-3:] == 'xls':
                                Transm = pd.read_excel(filename3)
                            else:
                                raise WrongLambdaFileFormat
                        
                    smoothTr = interp1d(Transm['wavelength (nm)'],Transm['transmission'], bounds_error = False )
                    topTr = pd.DataFrame(data=smoothTr(bot_se['wl']), columns = ['Tr_smooth'])
                    bot_se['fbot'] = topTr['Tr_smooth']
                    bot_se['fbot'] = bot_se['fbot'].replace(np.nan, 0)
                except:
                    raise FileMissing
                    
            #conditions for ideal coupling
            elif coupling == 'ideal':
                transmission = 100
                # Extracting bandgap for ideal and fixed coupling from top cell calculatedSE file
                try:
                    bandgap = top_se['Eg'][0]
                    top_bg_nm = round (1.24/bandgap * 1000)  
                except:
                    raise BandgapMissing
    
                # Define ftop for the top cell.
                #top_se['wl']=photon wavelength (nm)
                top_se.loc[top_se['wl'] < top_bg_nm, 'ftop'] = 1
                top_se.loc[top_se['wl'] >= top_bg_nm, 'ftop'] = 0
          
                # Assuming perfect transmission.
                bot_se.loc[bot_se['wl'] < top_bg_nm, 'fbot'] = 0
                bot_se.loc[bot_se['wl'] >= top_bg_nm, 'fbot'] = 1
                
            #conditions for fixed coupling
            elif coupling == 'fixed':
            
                # Extracting bandgap for ideal and fixed coupling from top cell calculatedSE file
                try:
                    bandgap = top_se['Eg'][0]
                    top_bg_nm = round (1.24/bandgap * 1000)  
                except:
                    raise BandgapMissing
    
                # Define ftop for the top cell.
                top_se.loc[top_se['wl'] < top_bg_nm, 'ftop'] = 1
                top_se.loc[top_se['wl'] >= top_bg_nm, 'ftop'] = 0
                
                # Assuming a fixed percentage of transmission.
                bot_se.loc[bot_se['wl'] < top_bg_nm, 'fbot'] = 0
                bot_se.loc[bot_se['wl'] >= top_bg_nm, 'fbot'] = transmission/100

  
            else:
                raise NoSelectionMade
        
            # FUNCTION OPPORTUNITY HERE? Need to numerically integrate each subcell.
            
            
            ##### 2T calculations 
            
            
            
            #constant spectralJscfactor
            q = 1.602176634*10**-19  # Elementary charge in SI unit
            c = 2.998*10**8
            h = 6.626*10**-34
            spectralJscfactor = (q/(h*c*10**9))   # 10**9 for conversion from m to nm
            
            
            #top cell for 2T tandem
            
            #### top QE
            top_se['QE'] = top_se['QE_smooth']*top_se['ftop']
            top_se['QE'] = top_se['QE'].replace(np.nan,0)
            
            
            #calculate Jsc from top QE
            EQE_1 = spectralJscfactor*top_se['wl']*top_se['QE']* top_se['I_sun']
            EQE_1= EQE_1.replace(np.nan,0)
            
            
            #Integrated Jsc from Spectally-resolved Jsc (top).....normalized from mA to A by dividing the electric current value by 1000
            QE_Integrated_top_Jsc= np.trapz(EQE_1)/1000
            QE_Integrated_top_Jsc = float("{:.2f}".format(QE_Integrated_top_Jsc))
            
            
            #bottom cell for 2T tandem
            
            
            #### bot QE
            bot_se['QE']=bot_se['QE_smooth']*bot_se['fbot']
            bot_se['QE'] = bot_se['QE'].replace(np.nan,0)
            
            
            #calculate Jsc from bottom cell
            EQE_2 = spectralJscfactor*bot_se['wl']*bot_se['QE']*bot_se['I_sun']
            EQE_2= EQE_2.replace(np.nan,0)
            
            
            #Integrated Jsc from Spectally-resolved Jsc (bottom).....normalized from mA to A by dividing the electric current value by 1000
            QE_Integrated_bot_Jsc= np.trapz(EQE_2)/1000
            QE_Integrated_bot_Jsc = float("{:.2f}".format(QE_Integrated_bot_Jsc))
           
            #passing information for calculated Jsc from EQE
            msg2="Top cell Jsc from " + str(top_name) + "'s QE" + " is " + '\033[1m' + str(QE_Integrated_top_Jsc) + '\033[0m' + " mA/cm^2" + "\nBottom cell Jsc from filtered " + str(bot_name) + "'s QE" + " is " + '\033[1m' + str(QE_Integrated_bot_Jsc) + '\033[0m' + " mA/cm^2"
            
            #passing QE information
            qe_list.append(msg2)
           
            
            #plotting QE
            if top_name not in cellNameList:
                ax2.plot(top_se['wl'], top_se['QE'], label="top QE " + top_name) 
            #this section needs to change if multiple transmission files are allowed -> get rid of if/elif and leave one line of ax2
            if coupling == 'lambda':
                ax2.plot(bot_se['wl'], bot_se['QE'], label=bot_name + " filtered through " + top_name)
            elif bot_name not in cellNameList:
                ax2.plot(bot_se['wl'], bot_se['QE'], label=bot_name + ", transmission = " + str(transmission) + "%")
            
            ######
            #single-junction SE
            top_se['1JSE']= top_se['SE']
            bot_se['1JSE']= bot_se['SE']
            
            #2T tandem SE factored by lower Jsc determined from QE of either top or bottom cells
            
            #top_se['ftop'] is top cell flag (1 above Eg (top) and 0 below Eg(top))
            top_se['2TSE']= top_se['SE']*top_se['ftop'] * (min(QE_Integrated_top_Jsc, QE_Integrated_bot_Jsc)/QE_Integrated_top_Jsc)
            #bot_se['ftop'] is bottom cell flag (0 above Eg (top) and 1 below Eg(top))
            bot_se['2TSE']= bot_se['SE']*bot_se['fbot'] * (min(QE_Integrated_top_Jsc, QE_Integrated_bot_Jsc)/QE_Integrated_top_Jsc)
            
            ####
            
            #2T tandem power factored by lower Jsc determined from QE of either top or bottom cells
            
            #top cell power
            top_se['2Teffic'] = top_se['SE']*top_se['ftop']*top_se['I_sun']
            top_se['2Teffic'] = top_se['2Teffic'] * (min(QE_Integrated_top_Jsc, QE_Integrated_bot_Jsc)/QE_Integrated_top_Jsc) #######is this eqn correct?
            #top_se['2Teffic'] = power/nm (W/nm*m"2)
            top_se['2Teffic'] = top_se['2Teffic'].replace(np.nan,0)
            #computing top cell efficiency = computing topcell power/standard incident spectrum power (1000 W/m^2)
            top_2T_eff = np.trapz(top_se['2Teffic'])/1000
            top_2T_eff = float("{:.2f}".format(top_2T_eff))
            
            
            ######
            #bottom cell power
            bot_se['2Teffic'] = bot_se['SE']*bot_se['fbot']*bot_se['I_sun']
            bot_se['2Teffic'] = bot_se['2Teffic'] * (min(QE_Integrated_top_Jsc, QE_Integrated_bot_Jsc)/QE_Integrated_top_Jsc)
            #bot_se['2Teffic'] = power/nm (W/nm*m"2)
            bot_se['2Teffic'] = bot_se['2Teffic'].replace(np.nan,0)
            #computing bottom cell efficiency = computing topcell power/standard incident spectrum power (1000 W/m^2)
            bot_2T_eff = np.trapz(bot_se['2Teffic'])/1000
            bot_2T_eff = float("{:.2f}".format(bot_2T_eff))
            
            
            ######
            #2T tandem efficiency based on SE
            Effi_2T = top_2T_eff+bot_2T_eff
            Effi_2T = float("{:.2f}".format(Effi_2T))
           
            
            #2T tandem efficiency calculations based on J-V curves of top and bottom cells
            
            #top cell
            #reading top cell current (i) and voltage (v)
            ivFOM_top = IV_Params(top_se.v, top_se.i)
            ivFOM_top.calc_iv_params()
            ivFOM_top =ivFOM_top.calc_iv_params()
            #top cell current from J-V curve
            Jsc_top=ivFOM_top['isc']
            
            
            #bottom cell
            #reading bottom cell current (i) and voltage (v)
            ivFOM_bot = IV_Params(bot_se.v, bot_se.i)
            ivFOM_bot.calc_iv_params()
            ivFOM_bot =ivFOM_bot.calc_iv_params()
            #bottom cell original current from J-V curve
            Jsc_bot=ivFOM_bot['isc']
            
            ##filtering bottom cell 
            
            #subtracting Jsc (integrated bottom cell QE in 2T tandem) from Jsc (J-V curve of bottom cells)
            delta_Jsc=Jsc_bot - QE_Integrated_bot_Jsc
            filtered_bot_i=bot_se.i-delta_Jsc
            
            
            #comparing top cell and filtered bottom cell Jsc 
            top_i = top_se.i
            bot_i = filtered_bot_i

            # Determine the call with the max current
            if np.max(top_i) > np.max(bot_i):
                added_current = bot_i
            else:
                added_current = top_i

            # equal spacing 
            new_i = np.linspace(0, added_current.max(), 5000)  #######5000 is probably more than enough

            # Interpolate each dataset
            bot_interp = interp1d(filtered_bot_i, bot_se.v, bounds_error=False) 
            top_interp = interp1d(top_se.i.replace(0,), top_se.v.replace(0,), bounds_error=False) #.replace(0,) needed to reduce stray lines in plotting

            #new voltage of top and bot
            new_bot_v = bot_interp(new_i) 
            new_top_v = top_interp(new_i)

            # Added voltage
            added_voltage = new_top_v + new_bot_v
            
            ####
            
            ##2T tandem efficiency calculations 
            
            #from 2T tandem I-V
            ivFOM_tandem = IV_Params(added_voltage, new_i)
            #ivFOM_tandem.calc_iv_params()
            ivFOM_tandem = ivFOM_tandem.calc_iv_params()  # why is this variable overwritten twice? Does the previous line do anything?
            #calculate max power point from 2T tandem I-V
            Tandem_Max_Power_Point=ivFOM_tandem['pmp']*10.0 #10.0 to convert from mW/cm^2 to W/m^2
            Tandem_Max_Power_Point=float("{:.2f}".format(Tandem_Max_Power_Point))
            
            
            #message for 2T tandem efficiency based on J-V curves 
            #msg5 = "Power Density at Max Power Point for 2T Tandem of " + str(top_name) + " on " + str(bot_name) + ": " + str(Tandem_Max_Power_Point) + " W/m^2"
            #msg5 = "" #optional edit - list pmp, voc, vmp, imp - sum top cell and bottom cell for each - need to check limiting current for 2T though 4/8/24
            #passing 2T tandem efficiency information
            #tandem_v.append(msg5)
          
            
            #plotting

            ax5.plot(added_voltage, new_i, label=("tandem I-V: " + top_name + " and " + bot_name))
            
            if bot_name not in cellNameList:
                ax6.plot(bot_se.v, bot_se.i, label=("raw bottom " + bot_name))  #doesn't seem to need zero replacement
                #.replace(0,np.nan,inplace=True) not allowed, None value not allowed   same for .replace(0,np.max(top_se.i),inplace=True)
#                ax6.plot(bot_se.v, bot_se.i, label=("raw bottom " + bot_name)) #.replace(0,np.nan,inplace=True) not allowed, None value not allowed   same for .replace(0,np.max(top_se.i),inplace=True) # .replace(0,'') throws error
            if coupling == 'lambda':
                ax6.plot(new_bot_v, new_i, label=(bot_name + " filtered through " + top_name))  
            elif bot_name not in cellNameList:
                ax6.plot(new_bot_v, new_i, label=("bottom " + bot_name + ", transmission = " + str(transmission) + "%"))
            if top_name not in cellNameList:
                ax6.plot(top_se.v.replace(0,), top_se.i.replace(0,), label=("raw top " + top_name)) 
            
            
            
            
            # MAKE THE FLAGS BOOLEAN SO TEST IS EASY!
            #Yes = 1, No = 0 flags. Set these for your desired output.
            #outputpath = '../tandems' # name of directory where we will save tandem data
            #Print data in this notebook
            print_data = 1
            #Save data to individual files
            save_ind = 1
            #Save all data to a single file
            save_all = 1
            #Print graphs in this notebook
            print_graphs = 1
                
            ##########
            if save_ind == 1:
            #create a compiled eff dataframe
                    
                if coupling == 'ideal':
                    transmission = 100
                elif coupling == 'lambda':
                    transmission = "light transmission as a function of wavelength" 

                tandem_data = [[top_name, top_2T_eff, coupling, transmission],[bot_name, bot_2T_eff],["2T_tandem", Effi_2T]]
                columns = ['cell name', 'efficiency', 'coupling', 'transmission']    
                tandemSEtemp = pd.DataFrame(tandem_data, columns=columns) 
                tandem_wl_data = [top_se['wl'],top_se['2Teffic'],top_se['2TSE'],top_se['1JSE'],top_se['QE'],bot_se['wl'],bot_se['2Teffic'],bot_se['2TSE'],bot_se['1JSE'],bot_se['QE'],top_se.v,top_se.i,bot_se.v, filtered_bot_i, bot_se.i,]
                wl_columns = ['top cell wavelength (nm)','top cell efficiency','top cell 2T SE', 'top cell 1J SE','top cell QE','bottom cell wavelength (nm)','bottom cell efficiency','bottom cell 2T SE', 'bottom cell 1J SE','bottom cell QE','top cell voltage','top cell current','bottom cell voltage','filtered bottom cell current','bottom cell current']
                wlResolvedData = pd.DataFrame({col: data for col, data in zip(wl_columns, tandem_wl_data)})
                wl_IV = [new_top_v,new_bot_v,added_voltage,new_i]
                wl_IV_columns = ['new top cell voltage','new bottom cell voltage','added voltage','new current']
                wlResolvedIV = pd.DataFrame({col: data for col, data in zip(wl_IV_columns, wl_IV)})
                tandemSE = pd.concat([tandemSEtemp,wlResolvedData,wlResolvedIV],axis=1)


                
            #Data saving
            if coupling=='fixed':
                i = 1
                while os.path.exists('../tandems_data/{}{:d}.csv'.format('2T_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
                    i += 1
                tandemSE.to_csv('../tandems_data/{}{:d}.csv'.format('2T_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i))
            else:
                i = 1
                while os.path.exists('../tandems_data/{}{:d}.csv'.format('2T_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
                    i += 1
                tandemSE.to_csv('../tandems_data/{}{:d}.csv'.format('2T_' + top_name + '+' + bot_name + '_' + coupling + '_', i))
            
            
            #Data Printing
            if print_data ==1:
                if coupling == 'fixed':
                    
                    msg3= '\033[1m' + "2T Tandem" + '\033[0m' + "\nCoupling Method: " + str(coupling) + "\nlight transmitted onto bottom cell: " + str(transmission) + "\nTop Cell is " + str(top_name) + "\n" + str(top_name) + " limiting efficiency = " + str(top_2T_eff) + "%" + "\nBottom Cell is " + str(bot_name) + "\n" + str(bot_name) + " limiting efficiency = " + str(bot_2T_eff) + "%" + "\nThe tandem efficiency for " + str(top_name) + " and " + str(bot_name) + " is "+ str(Effi_2T) + "%" + "\nThe results are saved in the 'tandems_data' folder and the file name is: 2T_tandem data_" + str(top_name) + '+' + str(bot_name) + '+' + str(coupling) + '+' + str(transmission) + '_' + str(i)
                    
                    
                    tandem_2T.append(msg3)
                    
                else:
                    
                    msg3= '\033[1m' + "2T Tandem" + '\033[0m' + "\nCoupling Method: " + str(coupling) + "\nlight transmitted onto bottom cell: " + str(transmission) + "%" + "\nTop Cell is " + str(top_name) + "\n" + str(top_name) + " limiting efficiency = " + str(top_2T_eff) + "%" + "\nBottom Cell is " + str(bot_name) + "\n" + str(bot_name) + " limiting efficiency = " + str(bot_2T_eff) + "%" + "\nThe tandem efficiency for " + str(top_name) + " and " + str(bot_name) + " is "+ str(Effi_2T) + "%" + "\nThe results are saved in the 'tandems_data' folder and the file name is: 2T_tandem data_" + str(top_name) + '+' + str(bot_name) + '+' + str(coupling) + '+' + '_' + str(i)
                    
                    tandem_2T.append(msg3)
                              
                
            # Data plotting        
            if print_graphs == 1:

                if coupling == 'lambda':
                    ax1.plot(bot_se['wl'], bot_se['2Teffic'], label=bot_name+' filtered through '+top_name)
                    ax3.plot(bot_se['wl'], bot_se['2TSE'], '--', label=bot_name+' filtered through '+top_name)
                elif bot_name not in cellNameList:
                    ax1.plot(bot_se['wl'], bot_se['2Teffic'], label=bot_name + ", transmission = " + str(transmission) + "%")
                    ax3.plot(bot_se['wl'], bot_se['2TSE'], '--', label=bot_name + ", transmission = " + str(transmission) + "%")

                if top_name not in cellNameList:
                    ax1.plot(top_se['wl'], top_se['2Teffic'], label=top_name)
                    ax3.plot(top_se['wl'], top_se['2TSE'], '--', label=top_name)
                    ax4.plot(top_se['wl'], top_se['1JSE'], '--', label=top_name)

                if bot_name not in cellNameList:
                    ax4.plot(bot_se['wl'], bot_se['1JSE'], '--', label=bot_name)

                
                      
            cellNameList.append(top_name)
            cellNameList.append(bot_name)  

    #Plotting axes limits and legends provided outside of loop
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax1.set_ylim(0,)

    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 

    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax3.set_ylim(0,)

    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax5.set_xlim(0,)
    ax5.set_ylim(0,)
    
    ax6.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax6.set_xlim(0,)
    ax6.set_ylim(0,)

    #Saving images
    if coupling=='fixed':
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig.savefig('../tandems_images/{}{:d}.png'.format('2T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig2.savefig('../tandems_images/{}{:d}.png'.format('2T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig3.savefig('../tandems_images/{}{:d}.png'.format('2T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig4.savefig('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig5.savefig('../tandems_images/{}{:d}.png'.format('2T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i)):
            i += 1
        fig6.savefig('../tandems_images/{}{:d}.png'.format('SJ_IV_' + top_name + '+' + bot_name + '_' + coupling + '_' + str(transmission) + '_', i), format="png", bbox_inches = "tight", dpi=500)
    else:
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig.savefig('../tandems_images/{}{:d}.png'.format('2T_power_density_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig2.savefig('../tandems_images/{}{:d}.png'.format('2T_EQE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig3.savefig('../tandems_images/{}{:d}.png'.format('2T_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig4.savefig('../tandems_images/{}{:d}.png'.format('SJ_SE_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('2T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig5.savefig('../tandems_images/{}{:d}.png'.format('2T_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
        i = 1
        while os.path.exists('../tandems_images/{}{:d}.png'.format('SJ_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i)):
            i += 1
        fig6.savefig('../tandems_images/{}{:d}.png'.format('SJ_IV_' + top_name + '+' + bot_name + '_' + coupling + '_', i), format="png", bbox_inches = "tight", dpi=500)
       
    return  Effi_2T, fig, fig2, fig3, fig4, fig5,  fig6, tandem_data, tandem_list, qe_list, tandem_2T, tandem_v
       



     
#Widgets part
        
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
        # List of selected fileswill be set to b.value
        b._files = filedialog.askopenfilename(multiple=True)
        #b._files.sort()
        b.files = [os.path.basename(file) for file in b._files]
        b.description = "Files Selected"
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"

    
  
    
    
def te_controls():
    
    
    def clear_plots(b):
        with outbox:
            clear_output()
        with outbox2:
            clear_output()

    #By recalling values from calc_tandem_4T_eff doing calculations using widgets buttons
    ### 4T Tandem
    def on_button_click(b):
        
        #top and bottom cells files
        top_file_list = topselect._files
        bot_file_list = botselect._files
        
        #top transmission file for lambda coupling
        lambdafile = trnsel3._files 
        
        
        ##activating coupling methods
        cpl3=coupl3.value
        
        #activating slider for fixed coupling
        tr3=trnsm3.value
        
        #using try except method, it does calculate when everything is right otherwise send an error message to the user
        try:
            tandem_eff, tandem_data, fig, fig2, fig3, fig4, fig5, tandem_list, tandem_4T, qe_list, tandem_iv = calc_tandem_4T_eff(top_file_list, bot_file_list, lambdafile, cpl3, tr3)
            
            with outbox:
                display(fig2)
                print(" ")
                #print(qe_list)
                print(Fore.BLACK +',\n'.join(qe_list))
                print(" ")
                display(fig4)
                print(" ")
                display(fig3)
                print("Note: Coupling methods are used to estimate the spectral efficiency of top and bottom cells in a 4T Tandem")
                print(" ")
                print(" ")
                #if  self.print_data_flag == True:
                print(Fore.BLACK +',\n'.join(tandem_4T))
                print(" ")
                display(fig)
                print(" ")
                display(fig5)
                print(" ")
                print(Fore.BLACK +',\n'.join(tandem_iv))
                print("------------------------------------------------------------------")
                
            with outbox2:
                clear_output() 
                print(Fore.RED +'See the 4T Tandem results below') 
                print("-----------------------------")
            
        #Sending error message
        except NoSelectionMade:
            with outbox2:
                clear_output() 
                print(Fore.RED +'No selection was made!\n Please select a coupling method')
                print("-----------------------------")
        
        except WrongFileFormat:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Files format are not correct,\n please see the instructions')
                print("-----------------------------")
                
        except WrongLambdaFileFormat:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Top transmission file format is not correct,\n please see the instructions')
                print("-----------------------------")
                
        except FileMissing:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Top transmission file missing,\n please see the instructions')
                print("-----------------------------")
                
        except BandgapMissing:
            with outbox2:
                clear_output() 
                print(Fore.RED +"Bandgap of top cell missing!\nPlease manually enter the bandgap of top cell \nby designating a column's header as Eg \nand entering the value in the column's second row")
                print("-----------------------------")
   
    
    #By recalling values from calc_tandem_4T_eff doing calculations using widgets buttons
    ### 2T Tandem
    
    def on_button_click_2T(b):
        
        #top and bottom cells files
        top_file_list = topselect._files
        bot_file_list = botselect._files
        
        #top transmission file for lambda coupling
        lambdafile = trnsel3._files 
        
        
        ##activating coupling methods
        cpl3=coupl3.value
        
        #activating slider for fixed coupling
        tr3=trnsm3.value
        
        #using try except method, it does calculate when everything is right otherwise send an error message to the user
        try:
            Effi_2T, fig, fig2, fig3, fig4, fig5, fig6, tandem_data, tandem_list, qe_list, tandem_2T, tandem_v = calc_tandem_2T_eff(top_file_list, bot_file_list, lambdafile, cpl3, tr3)
            
            with outbox:
                display(fig2)
                print(" ")
                print(Fore.BLACK +',\n'.join(qe_list))
                print(" ")
                display(fig4)
                print(" ")
                display(fig3)
                print("Note: Coupling methods and the lower Jsc of either cell are used to estimate the spectral efficiency of top and bottom cells in a 2T Tandem")
                print(" ")
                print(" ")
                #if  self.print_data_flag == True:
                print(Fore.BLACK +',\n'.join(tandem_2T))
                print(" ")
                display(fig)
                print(" ")
                display(fig6)
                print(" ")
                display(fig5)
                print(" ")
                print(Fore.BLACK +',\n'.join(tandem_v))
                print(" ")
                print("------------------------------------------------------------------")
                
            with outbox2:
                clear_output() 
                print(Fore.RED +'See the 2T Tandem results below') 
                print("-----------------------------")
            
        #Sending error message
        except NoSelectionMade:
            with outbox2:
                clear_output() 
                print(Fore.RED +'No selection was made!\n Please select a coupling method')
                print("-----------------------------")
        
        except WrongFileFormat:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Files format are not correct,\n please see the instructions')
                print("-----------------------------")
                
        except WrongLambdaFileFormat:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Top transmission file format is not correct,\n please see the instructions')
                print("-----------------------------")
                
        except FileMissing:
            with outbox2:
                clear_output() 
                print(Fore.RED +'Top transmission file missing,\n please see the instructions')
                print("-----------------------------")
                
        except BandgapMissing:
            with outbox2:
                clear_output() 
                print(Fore.RED +"Bandgap of top cell missing!\nPlease manually enter the bandgap of top cell \nby designating a column's header as Eg \nand entering the value in the column's second row")
                print("-----------------------------")
    
    
    
   
#    def updateWidget(*args)
        
    
    
    '''
    use interactive_input for GUI in IPython
    '''
        
    cell_layout = widgets.Layout(display='inline_flex',
        flex_flow='row',
        justify_content='flex-end',
        width='320px' )                        

    # controls 

    # Column 1
    s1label = widgets.Label(value='Step 1: Select cells')
    
    toplabel = widgets.Label(value='Select Top Cell(s)')
    # Button for selecting top cell(s)
    topselect = SelectFilesButton(
        description = 'Select Files', layout=cell_layout
    )
    top_file_list = topselect.files

    toptext = widgets.Select(options = [], layout=cell_layout)
    toplink = widgets.link((topselect, 'files'), (toptext, 'options'))
                               
    botlabel = widgets.Label(value='Select Bottom Cell(s)')
    # Button for selecting top cell(s)
    botselect = SelectFilesButton(
        description = 'Select Files', layout=cell_layout
    )

    bot_file_list = botselect.files

    bottext = widgets.Select(options = [], layout=cell_layout)
    botlink = widgets.link((botselect, 'files'), (bottext, 'options'))
    
    
    # Column 2
    s2label = widgets.Label(value='Step 2: Select 1 from 3 coupling methods')
    # Set up coupling options
    cplopt = [('Ideal coupling, cutoff at top cell bandgap', 'ideal'),
              ('Fixed transmission percent to bottom cell', 'fixed'),
              ('Use top cell transmission file', 'lambda')]
    cploptN = [('No selection made', 'none')] + cplopt


    c3label = widgets.Label(value='Coupling methods:')    
    coupl3 = widgets.Dropdown(
        options=cploptN, 
        value='lambda', disabled=False,
        layout=cell_layout
    )
    
    
    #Slider for fixed coupling transmission
    trnsm3 = widgets.IntSlider(
        value=90, min=0, max=100, step=1, 
        description='Fixed T %:', disabled=False, 
        continuous_update=False, orientation='horizontal', 
        readout=True, readout_format='d',
        layout=cell_layout
    )
        
    #Button for importing lambda file
    trnsel3 = SelectFilesButton(
        description="Select Transmission Files for lambda Coupling", layout=cell_layout
    )
    lambdafile = trnsel3.files

    lambdatext = widgets.Select(options = [], layout=cell_layout)
    lambdalink = widgets.link((trnsel3, 'files'), (lambdatext, 'options'))

    # Column 3
    #s3label = widgets.Label(value='Step 3: Select output options')
    
    s3label = widgets.Label(value='Step 3: Select calculate options')
    
    
    
    # Actions   
  
    # button to start calculations
    calc = widgets.Button(
        description="Calculate 4T Tandem",
        layout=cell_layout)
    
    calc_2 = widgets.Button(
        description="Calculate 2T Tandem",
        layout=cell_layout)
    
    clr = widgets.Button(
        description="Clear Output Results",
        layout=cell_layout)


    #displaying message using .Output method in jupyter lab
    outbox = widgets.Output()
    outbox2 = widgets.Output()
    
    leftcntrls = [s1label, toplabel, topselect, toptext,
                  botlabel, botselect, bottext
                 ]
                
    ctrcntrls = [s2label, c3label, coupl3, trnsm3, trnsel3, lambdatext
                ]
    
    
    rtcntrls = [s3label, calc, calc_2, clr, outbox2
               ]

    dncntrls = [outbox]
    
    # call function to start calculations
    calc.on_click(on_button_click)
    calc_2.on_click(on_button_click_2T)
    clr.on_click(clear_plots)
    
        
     # user interface        
    box_layout = widgets.Layout(display='flex',
        flex_flow='column',
        align_items='center',
        width='33.33%',
        height = '400px')

    
    #for output results
    box_layout_2 = widgets.Layout(display='flex',
        flex_flow='column',
        align_items='center',
        width='100%',
        height = '3000px')
                           
    leftbox = widgets.VBox(leftcntrls,layout=box_layout)
    ctrbox = widgets.VBox(ctrcntrls, layout=box_layout)
    rtbox = widgets.VBox(rtcntrls, layout=box_layout)
    dnbox = widgets.VBox(dncntrls, layout=box_layout_2)
    
    
    
    ui = widgets.HBox([leftbox, ctrbox, rtbox])
    ui2 = widgets.HBox([dnbox])
    
    
    display(ui, ui2)

if __name__=="__main__":
    te_controls()
