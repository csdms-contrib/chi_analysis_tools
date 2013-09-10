
## chi_plot_multiple_mn.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This function creates a chi-elevation plot (see Perron and Royden, 2012; An 
## integral approach to bedrock river profile analysis), for river profiles
## stored in the *.tree file. The trunk channel is highlighted distinctly from
## the tributary channels.  The user needs to specify the FileName,
## OutputFigureName and OutputFIgureFormat in the code.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## DTM 07/11/2012
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## 
##
## This modified script now creates multiple plots for all channel tree files contained
## within the directory. Place all the *.tree files to be plotted in the same 
## directory with the python script and run from there.
## Change the data variables (l.68) to plot different paramaters
##
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## DAV 06/08/2013
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Modified to plot results of a sensitivity analysis done for the Child 
## model runs
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 25/08/2013
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



#import modules
import numpy as np, matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
#import re

from glob import glob
from matplotlib import rcParams


#==============================================================================
# numbers = re.compile(r'(\d+)')               
# def numericalSort(value):                    
#     parts = numbers.split(value)             
#     parts[1::2] = map(int, parts[1::2])      
#     return parts
# This function enables numerical sorting of files
#==============================================================================

def make_plots():
    
    
    DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\Child_runs\\'
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\PA\\'
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Child_runs_mk4\\'
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Child_runs_mk4\\sensitivity_skip\\' 
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Child_runs_mk4\\sensitivity_tn\\' 
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Child_runs_mk4\\sensitivity_msl\\'     
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Child_runs_mk4\\sensitivity_sigma\\' 
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Test_data\\'    
    
    # set sizes for fonts. These sizes are designed to be visible on a powerpoint
    # slide (10 inches x 7.5 inches) and at the same time be legible if this
    # slide is shrunk to 1 column width (85mm)
    label_size = 20
    title_size = 15
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    #########################
    #                       #
    #   READ IN THE DATA    #
    #                       #
    #########################

    for fname in glob(DataDirectory+"*_fullProfileMC_forced_*.tree"): # assigns a number to each iteration (i.e. for every .tree file in the directory)
    #for fname in glob(DataDirectory+"*_mover_n_*"): # assigns a number to each iteration (i.e. for every .tree file in the directory)
        #print idx, fname 
        
        split_fname = fname.split('\\')
        no_tree_levs = len(split_fname) 
        this_fname = split_fname[no_tree_levs-1]
        print "fname: "+fname+" and this fname: "+this_fname
        
        fname_prefix = this_fname.split('.')[0]
       
        # open file
        f = open(fname,'r')
        lines = f.readlines()[0:]   # read in the data           
        no_lines = len(lines)-1   # get the number of lines (=number of data)
        
        # get the parameter values
        sigma_value = this_fname.split('_')[5]    #get the value of sigma used
        skip_value = this_fname.split('_')[6]    #get the value of skip used 
        msl_value = this_fname.split('_')[7]    #get the value of minimum_segment_length used 
        tn_value = this_fname.split('_')[8]    #get the value of target_nodes used 
    
        # get the A_0 and m/n values
        line = lines[0].strip().split(" ")
        A_0 = float(line[0]) 
        m_over_n = float(line[1]) 
        
        print "Running a file!"
        print "A_0: " +str(A_0)+" m_over_n: "+str(m_over_n)

        # set the name of the output figure
        OutputFigureName = DataDirectory+fname_prefix
        OutputFigureFormat = 'pdf'
        print "The output figure name is: "
        print OutputFigureName
    
        # data variables
        channel_id = np.zeros(no_lines, dtype=np.int)     # ID number for channel segment
    
        chi = np.zeros(no_lines)            # chi
        elev = np.zeros(no_lines)           # elevation
        fitted_elevation_mean = np.zeros(no_lines) #fitted elevation
        m_mean = np.zeros(no_lines)
    
        for i in range (0,no_lines):
            line = lines[i+1].strip().split(" ")
            #print line
            channel_id[i] = int(float(line[0]))
            chi[i] = float(line[7])
            elev[i] = float(line[8])
            fitted_elevation_mean[i] = float(line[20])  
            m_mean[i] = float(line[11])         # the M_chi
        f.close()
  
        # get the segments  
        n_data = no_lines
        n_channel_segments = channel_id[n_data-1]+1
        segment_lengths = np.zeros(n_channel_segments, dtype=np.int)
        for i in range (0,n_data):
            segment_lengths[channel_id[i]] = segment_lengths[channel_id[i]] + 1   
 
        # SET UP COLOURMAPS
        jet = plt.get_cmap('jet')
    
        # channel ID
        Channel_ID_MIN = np.min(channel_id)
        Channel_ID_MAX = np.max(channel_id)
        cNorm_channel_ID  = colors.Normalize(vmin=Channel_ID_MIN, vmax=Channel_ID_MAX)  # the max number of channel segs is the 'top' colour
        scalarMap_channel_ID = cmx.ScalarMappable(norm=cNorm_channel_ID, cmap=jet)
           
        #########################
        #                       #
        #   MAKE CHI-PLOTS      #
        #                       #
        #########################    
        
        # set the size to that of a powerpoint slide
        plt.figure(1, facecolor='white',figsize=(10,7.5))
        ax = plt.subplot(1,1,1)
        # Set up colourmap for plotting channel segments
        jet = plt.get_cmap('jet')
        Channel_ID_MIN = np.min(channel_id)
        Channel_ID_MAX = np.max(channel_id)
        cNorm_channel_ID  = colors.Normalize(vmin=Channel_ID_MIN, vmax=Channel_ID_MAX)  # the max number of channel segs is the 'top' colour
        scalarMap_channel_ID = cmx.ScalarMappable(norm=cNorm_channel_ID, cmap=jet)
    
        
        # loop through the data, extracting each profile segment in turn   
        # this plots the data in chi-elevation space
        data_pointer = 0    # points to data element in chi and elev vectors
        for i in range (0,len(segment_lengths)):
            chi_seg = np.zeros(segment_lengths[i])
            elev_seg = np.zeros(segment_lengths[i])
        
            for j in range (0,segment_lengths[i]):
                chi_seg[j] = chi[data_pointer]
                elev_seg[j] = elev[data_pointer]
                data_pointer = data_pointer + 1

            if i == 0:
                # plot trunk stream in black, with thicker line 
                l1, = ax.plot(chi_seg, elev_seg, "k-", linewidth=4)

      
        # Configure final plot
        plt.xlabel('$\chi$ (m)', fontsize= axis_size)
        plt.ylabel('Elevation (m)', fontsize =axis_size)
        plt.title('$A_0$: '+str(A_0)+' m$^2$,  $m/n$: '+str(m_over_n)+ " sigma: "+ sigma_value+" skip: " +skip_value+" msl: " +msl_value+" tn: "+tn_value, fontsize = title_size)  
        
        ax.spines['top'].set_linewidth(3)
        ax.spines['left'].set_linewidth(3)
        ax.spines['right'].set_linewidth(3)
        ax.spines['bottom'].set_linewidth(3) 
        ax.tick_params(axis='both', width=3)     
        
        plt.savefig(OutputFigureName+'_profiles.' + OutputFigureFormat, format=OutputFigureFormat)
        plt.clf()
 


   
if __name__ == "__main__":
    make_plots()