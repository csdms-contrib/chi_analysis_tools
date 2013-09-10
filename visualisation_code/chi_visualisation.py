## chi_visualisation.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This function creates a chi-elevation plot (see Perron and Royden, 2012; An 
## integral approach to bedrock river profile analysis), for river profiles
## stored in the *.tree file. The trunk channel is highlighted distinctly from
## the tributary channels.  The user needs to specify the FileName,
## OutputFigureName and OutputFIgureFormat in the code.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## DTM 07/11/2012
## Modified by SMM 07/08/2013
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#import modules
import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as colors
import matplotlib.cm as cmx

def make_plots():
    

    label_size = 20
    title_size = 18
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

    # open file
    #DataDirectory = 'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\India\\'
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\PA\\'
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Test_data\\'  
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Test_data\\'
    DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Apennines\\'
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Child_runs_mk4\\'   
    #DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\PA\\'  
    #DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\Apennines\\' 
    #FileName = 'pa_basin_fullProfileMC_mainstem_1189.tree'
    #FileName = 'pa_basin_fullProfileMC_mainstem_3124.tree'    
    #FileName = 'rio_torto_fullProfileMC_colinear_633.tree'
    #FileName = 'rio_torto_fullProfileMC_mainstem_114.tree'
    #FileName = 'rio_torto_fullProfileMC_colinear_114.tree'  
    #FileName = 'rio_torto_fullProfileMC_colinear_110.tree' 
    #FileName = 'rio_torto_fullProfileMC_mainstem_110.tree' 
    #FileName = 'rio_torto_fullProfileMC_forced_0.4_110.tree' 
    #FileName = 'UniUplift09b_mk4_t6_fullProfileMC_colinear.tree' 
    #FileName = 'pa_for_chi_fullProfileMC_forced_0.85_3_2_15_100_384.tree'
    #FileName = 'rio_torto_fullProfileMC_forced_0.7_20_2_15_100_110.tree'
    FileName = 'rio_torto_fullProfileMC_forced_0.6_20_1_12_100_110.tree'
  
    #OutputFigureName = 'chi_plot'
    #OutputFigureFormat = 'eps'
    fname = DataDirectory + FileName
    f = open(DataDirectory + FileName,'r')  # open file
    lines = f.readlines()   # read in the data
    n_lines = len(lines)   # get the number of lines (=number of data)
    n_data = n_lines -1
    # data variables
    channel_id = np.zeros(n_data)     # ID number for channel segment
    receiver_channel = np.zeros(n_data)
    #node_on_receiver_channel = np.zeros(n_data, dtype=np.int)
    node = np.zeros(n_data)           # node number
    row = np.zeros(n_data)            # row
    col = np.zeros(n_data)            # column
    flow_dist = np.zeros(n_data)      # flow distance
    chi = np.zeros(n_data)            # chi
    elevation = np.zeros(n_data)           # elevation
    drainage_area = np.zeros(n_data)  # drainage area
    n_data_points_uic = np.zeros(n_data)
    m_mean = np.zeros(n_data)
    #m_standard_deviation = np.zeros(n_data)
    m_standard_error = np.zeros(n_data)
    b_mean = np.zeros(n_data)
    #b_standard_deviation = np.zeros(n_data)
    b_standard_error = np.zeros(n_data)
    DW_mean = np.zeros(n_data)
    #DW_standard_deviation = np.zeros(n_data)
    DW_standard_error = np.zeros(n_data)
    fitted_elevation_mean = np.zeros(n_data)
    #fitted_elevation_standard_deviation = np.zeros(n_data)
    fitted_elevation_standard_error = np.zeros(n_data)
    
    print "Reading " + FileName

    # a few lines for splitting up the filename and getting the parameter values
    # out. This is slightly redundant since you actually have to write the filename
    # in in the first place but these data extraction methods are included so 
    # that any derivative functions can loop through these files
    # and extract the relevant parameter values. 
    # first seperate the filename from the path
    split_fname = fname.split('\\')
    no_tree_levs = len(split_fname) 
    this_fname = split_fname[no_tree_levs-1]
    print "fname: "+fname+" and this fname: "+this_fname
   
    # now seperate the fname prefix from the extension     
    split_fname = fname.split('.')    
    no_tree_levs = len(split_fname) 
    fname_prefix = split_fname[no_tree_levs-2]
    
    print fname_prefix
    
    # now split once more based on the junction, etc. 
    split_fname = fname_prefix.split('_')  
    no_tree_levs = len(split_fname) 
    
    jn_str = split_fname[no_tree_levs-1]
    tn_str = split_fname[no_tree_levs-2]
    msl_str = split_fname[no_tree_levs-3]
    skip_str = split_fname[no_tree_levs-4]
    sigma_str = split_fname[no_tree_levs-5]
    
    param_str = "junction: " + jn_str + " $\sigma$: " + sigma_str + " skip: " + skip_str+ " msl: " + msl_str + " tn: " + tn_str

    
    # get the A_0 and m/n values
    line = lines[0].strip().split(" ")
    A_0 = float(line[0]) 
    m_over_n = float(line[1]) 
    
    for i in range (0,n_data):
        line = lines[i+1].strip().split(" ")
        #print line
        channel_id[i] = int(float(line[0]))
        receiver_channel[i] = int(float(line[1])) 
        #node_on_receiver_channel = int(float(line[2]))
        node[i] = int(float(line[3]))
        row[i] = int(line[4])
        col[i] = int(line[5])
        flow_dist[i] = float(line[6])
        chi[i] = float(line[7])
        elevation[i] = float(line[8])
        drainage_area[i] = float(line[9])
        n_data_points_uic[i] = float(line[10])
        m_mean[i] = float(line[11])
        #m_standard_deviation[i] = float(line[12])
        m_standard_error[i] = float(line[13])
        b_mean[i] = float(line[14])
        #b_standard_deviation[i] = float(line[15])
        b_standard_error[i] = float(line[16])
        DW_mean[i] = float(line[17])
        #DW_standard_deviation[i] = float(line[18])
        DW_standard_error[i] = float(line[19])
        fitted_elevation_mean[i] = float(line[20])
        #fitted_elevation_standard_deviation = float(line[21])
        fitted_elevation_standard_error[i] = float(line[22])
    f.close()
    
    # Determine number of segments, and their respective lengths
    n_channel_segments = channel_id[n_data-1]+1
    segment_lengths = np.zeros(n_channel_segments, dtype=np.int)
    for i in range (0,n_data):
        segment_lengths[channel_id[i]] = segment_lengths[channel_id[i]] + 1  
        
    #########################
    #                       #
    #   MAKE CHI-PLOTS      #
    #                       #
    #########################    
    print "Producing figures..."
    
    # SET UP COLOURMAPS
    jet = plt.get_cmap('jet')
    hot = plt.get_cmap('RdYlBu_r')
    
    
    # m-values
    m_MIN = np.min(m_mean)
    m_MAX = np.max(m_mean)
    cNorm_m_values  = colors.Normalize(vmin=m_MIN, vmax=m_MAX)  # the max number of channel segs is the 'top' colour
    scalarMap_m_values = cmx.ScalarMappable(norm=cNorm_m_values, cmap=hot)

    # channel ID
    Channel_ID_MIN = np.min(channel_id)
    Channel_ID_MAX = np.max(channel_id)
    cNorm_channel_ID  = colors.Normalize(vmin=Channel_ID_MIN, vmax=Channel_ID_MAX)  # the max number of channel segs is the 'top' colour
    scalarMap_channel_ID = cmx.ScalarMappable(norm=cNorm_channel_ID, cmap=jet)

    minfd = min(flow_dist)
    print "minimum flow distance: " +str(minfd)
    
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # PLOT OF CHI m-VALUE (+ERROR) AGAINST CHI WITH EACH CHANNEL LABELLED WITH DISTINCT COLOUR
    plt.figure(2, facecolor='white',figsize=(10,7.5))    
    ax = plt.subplot(1,1,1)
    for i in range (0,len(channel_id)):
        
        colorVal = scalarMap_channel_ID.to_rgba(channel_id[i])  
        if channel_id[i]==0:
            plt.plot(chi[i], m_mean[i], "o", markersize=10, color='black', markeredgecolor = 'black',alpha=0.7)  
        else:
            plt.plot(chi[i], m_mean[i], "o", markersize=8, color=colorVal, markeredgecolor = colorVal,alpha = 0.7)
                 
           # plt.plot(chi[i], m_mean[i] + m_standard_error[i], ".", linewidth=0.25, color='k')
           # plt.plot(chi[i], m_mean[i] - m_standard_error[i], ".", linewidth=0.25, color='k')
            
    # Configure final plot
    ax.spines['top'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['right'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5) 
    ax.tick_params(axis='both', width=2.5)    
    plt.xlabel('$\chi$ (m)', fontsize = axis_size)
    plt.ylabel('Gradient in $\chi$ space', fontsize = axis_size)
    plt.title('$A_0$: '+str(A_0)+' m$^2$, and $m/n$: '+str(m_over_n)+ ' '+param_str, fontsize = title_size)

    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # LONGITUDINAL PROFILES WITH COLOUR SCALE GIVING CHI-m VALUE
    plt.figure(3, facecolor='white',figsize=(10,7.5))  
    ax = plt.subplot(1,1,1)    
    for i in range (0,len(channel_id)):
            colorVal = scalarMap_m_values.to_rgba(m_mean[i])
            plt.plot((flow_dist[i]-minfd)/1000, elevation[i], "o", markersize=6, color=colorVal, markeredgecolor = colorVal,alpha = 0.5)       
                  
    # Configure final plot
    sm = plt.cm.ScalarMappable(cmap=hot,norm=plt.normalize(vmin=np.min(m_mean), vmax=np.max(m_mean)))
    sm._A = []
    
    cbar = plt.colorbar(sm,orientation='horizontal',use_gridspec=True)
    cbar.set_label('Gradient in $\chi$ space', fontsize = axis_size)
    plt.xlabel('Distance upstream (km)', fontsize = axis_size)
    plt.ylabel('Elevation (m)', fontsize = axis_size)
    plt.title('$A_0$: '+str(A_0)+' m$^2$, and $m/n$: '+str(m_over_n)+' '+param_str, fontsize = title_size)
    
    ax.spines['top'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['right'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5) 
    ax.tick_params(axis='both', width=2.5)   
    
    cbar.ax.spines['top'].set_linewidth(2.5)
    cbar.ax.spines['left'].set_linewidth(2.5)
    cbar.ax.spines['right'].set_linewidth(2.5)
    cbar.ax.spines['bottom'].set_linewidth(2.5) 
    cbar.ax.tick_params(axis='both', width=2.5)       
    
    plt.tight_layout()

    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # BASIC CHI-PLOT WITH EACH CHANNEL LABELLED WITH A DIFFERENT COLOUR
    plt.figure(1, facecolor='white',figsize=(10,7.5))
    ax = plt.subplot(1,1,1)

    
    data_pointer = 0    # points to data element in chi and elev vectors
    for i in range (0,len(segment_lengths)):
        chi_seg = np.zeros(segment_lengths[i])
        elev_seg = np.zeros(segment_lengths[i])
        
        for j in range (0,segment_lengths[i]):
            chi_seg[j] = chi[data_pointer]
            elev_seg[j] = elevation[data_pointer]
            data_pointer = data_pointer + 1

        if i == 0:
            # plot trunk stream in black, with thicker line 
            l1, = ax.plot(chi_seg, elev_seg, "k-", linewidth=4, alpha = 0.3)
        else:
            # plot other stream segments plot
            colorVal = scalarMap_channel_ID.to_rgba(i) # this gets the distinct colour for this segment
            ax.plot(chi_seg, elev_seg, "-", linewidth=4, color=colorVal, alpha = 0.3) 


   # now loop again plotting the fitted segments    
    data_pointer = 0    # points to data element in chi and elev vectors
    for i in range (0,len(segment_lengths)):
        chi_seg = np.zeros(segment_lengths[i])
        elev_seg = np.zeros(segment_lengths[i])
        
        for j in range (0,segment_lengths[i]):
            chi_seg[j] = chi[data_pointer]
            elev_seg[j] = fitted_elevation_mean[data_pointer]
            data_pointer = data_pointer + 1

        if i == 0:
            # plot trunk stream in black, with thicker line 
            l2, = plt.plot(chi_seg, elev_seg, "k", linewidth=3, dashes=(10,2))
        else:
            # plot other stream segments plot
            colorVal = scalarMap_channel_ID.to_rgba(i) # this gets the distinct colour for this segment
            ax.plot(chi_seg, elev_seg, linewidth=3, color=colorVal, dashes=(10,2)) 
    
    # Configure final plot
    plt.xlabel('$\chi$ (m)', fontsize= axis_size)
    plt.ylabel('Elevation (m)', fontsize =axis_size)
    plt.title('$A_0$: '+str(A_0)+' m$^2$, and $m/n$: '+str(m_over_n)+' '+param_str, fontsize = title_size)  
    plt.legend((l1, l2), ('Data', 'Best fit segments'), 'lower right',prop={'size':label_size}) 
    
    ax.spines['top'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['right'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5) 
    ax.tick_params(axis='both', width=2.5)  


    #plt.savefig(OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat)


    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#    # CHI-PLOT USING FITTED PROFILES (+ERROR) RATHER THAN RAW DATA, LABELLED WITH DISTINCT COLOUR
#    plt.figure(4, facecolor='white')    
#    data_pointer = 0    # points to data element in chi and elev vectors
#    for i in range (0,len(segment_lengths)):
#        chi_seg = np.zeros(segment_lengths[i])
#        fitted_elev_seg = np.zeros(segment_lengths[i])
#        fitted_elev_ulim = np.zeros(segment_lengths[i])        
#        fitted_elev_llim = np.zeros(segment_lengths[i])
#        
#        for j in range (0,segment_lengths[i]):
#            chi_seg[j] = chi[data_pointer]
#            fitted_elev_seg[j] = fitted_elevation_mean[data_pointer]            
#            fitted_elev_ulim[j] = fitted_elevation_mean[data_pointer] + fitted_elevation_standard_error[data_pointer]            
#            fitted_elev_llim[j] = fitted_elevation_mean[data_pointer] - fitted_elevation_standard_error[data_pointer]
#            data_pointer = data_pointer + 1
#            # plot other stream segments plot
#        colorVal = scalarMap_channel_ID.to_rgba(i)
#        plt.plot(chi_seg, fitted_elev_ulim, "-k", linewidth=0.2)
#        plt.plot(chi_seg, fitted_elev_llim, "-k", linewidth=0.2)
#        plt.plot(chi_seg, fitted_elev_seg, "-", linewidth=0.5, color=colorVal) 
#            
#    # Configure final plot
#    plt.xlabel('$\chi$ / m')
#    plt.ylabel('fitted elevation / m')
#
#    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#    # CHI-PLOT DISPLAYING BOTH RAW DATA AND FITTED PROFILES
#    plt.figure(5, facecolor='white')    
# 
#    data_pointer = 0    # points to data element in chi and elev vectors
#    for i in range (0,len(segment_lengths)):
#        chi_seg = np.zeros(segment_lengths[i])
#        elev_seg = np.zeros(segment_lengths[i])
#        fitted_elev_seg = np.zeros(segment_lengths[i])
#        
#        for j in range (0,segment_lengths[i]):
#            chi_seg[j] = chi[data_pointer]
#            elev_seg[j] = elevation[data_pointer]
#            fitted_elev_seg[j] = fitted_elevation_mean[data_pointer]
#            data_pointer = data_pointer + 1
#
#        # plot other stream segments plot
#        colorVal = scalarMap_channel_ID.to_rgba(i) # this gets the distinct colour for this segment
#        plt.plot(chi_seg, elev_seg, "-k", linewidth=1)
#        plt.plot(chi_seg, fitted_elev_seg, "-", linewidth=1, color=colorVal)
#        
#     # Configure final plot
#    plt.xlabel('$\chi$ / m')
#    plt.ylabel('elevation / m')
#
#    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#    # CHI-PLOT USING ONLY MAIN TRIBUTARY WITH FITTED PROFILES (+ERROR)
#    plt.figure(6, facecolor='white')    
#    for i in range (0,len(channel_id)):   
#        plt.plot(chi[channel_id==0], fitted_elevation_mean[channel_id==0] + fitted_elevation_standard_error[channel_id==0], "-k", linewidth=0.2)
#        plt.plot(chi[channel_id==0], fitted_elevation_mean[channel_id==0] - fitted_elevation_standard_error[channel_id==0], "-k", linewidth=0.2)
#        plt.plot(chi[channel_id==0], elevation[channel_id==0], "-b", linewidth=0.5)            
#        plt.plot(chi[channel_id==0], fitted_elevation_mean[channel_id==0], "-r", linewidth=0.5)
#
#
#
#     # Configure final plot
#    plt.xlabel('$\chi$ / m')
#    plt.ylabel('elevation / m')
#
#    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#    # PLOT OF CHI b-VALUE (+ERROR) AGAINST CHI WITH EACH CHANNEL LABELLED WITH DISTINCT COLOUR
#    plt.figure(7, facecolor='white')    
#    for i in range (0,len(channel_id)):
#            colorVal = scalarMap_channel_ID.to_rgba(channel_id[i])
#            plt.plot(chi[i], b_mean[i], ".", linewidth=0.2, color=colorVal)            
#            #plt.plot(chi[i], b_mean[i] + b_standard_error[i], ".", linewidth=0.25, color='k')
#            #plt.plot(chi[i], b_mean[i] - b_standard_error[i], ".", linewidth=0.25, color='k')
#            
#    # Configure final plot
#    plt.xlabel('$\chi$ / m')
#    plt.ylabel('b-value / m')
#    
#    
#    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#    # PLOT OF DW STATISTIC (+ERROR) AGAINST CHI WITH EACH CHANNEL LABELLED WITH DISTINCT COLOUR
#    plt.figure(8, facecolor='white')    
#    for i in range (0,len(channel_id)):
#            colorVal = scalarMap_channel_ID.to_rgba(channel_id[i])
#            plt.plot(chi[i], DW_mean[i], ".", linewidth=0.2, color=colorVal)            
#            #plt.plot(chi[i], DW_mean[i] + DW_standard_error[i], ".", linewidth=0.25, color='k')
#            #plt.plot(chi[i], DW_mean[i] - DW_standard_error[i], ".", linewidth=0.25, color='k')
#            
#    # Configure final plot
#    plt.xlabel('$\chi$ / m')
#    plt.ylabel('DW statistic')
    
    plt.show()
    
if __name__ == "__main__":
    make_plots()