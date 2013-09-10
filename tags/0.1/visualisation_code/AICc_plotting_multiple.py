## AICc_plotting_multiple.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This function creates plots that allow inspection of the best fit
## AICc values and their standard deviations for the collinearity test
## as well as for individual channels
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 04/08/2013
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#import modules
import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as colors
import matplotlib.cm as cmx
from glob import glob

def make_plots():
 

    label_size = 20
    title_size = 18
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size   
    
    OutputFigureFormat = 'png'

    #########################
    #                       #
    #   READ IN THE DATA    #
    #                       #
    #########################

    # open file
    #DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\PA\\'    
    #DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\Apennines\\torrente_lapa\\'
    DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\Apennines\\'
    # DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\PA\\'
    # DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Apennines\\'
    #DataDirectory =  'm:\\topographic_tools\\LSDRaster_chi_package\\Test_data\\'       
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Test_data\\'
    #DataDirectory =  'c:\\code\\topographic_analysis\\LSDRaster_chi_package\\Child_runs_mk4\\'


    for FileName in glob(DataDirectory+"*BF_movern_*.movern"): # assigns a number to each iteration (i.e. for every .tree file in the directory)

        split_fname = FileName.split('\\')
        no_tree_levs = len(split_fname) 
        this_fname = split_fname[no_tree_levs-1]
        print "fname: "+FileName+" and this fname: "+this_fname
        
        fname_prefix = this_fname.split('.')[0]
  
        fname = FileName
        f = open(FileName,'r')  # open file
        lines = f.readlines()   # read in the data
        n_lines = len(lines)   # get the number of lines (=number of data)
        #n_data = n_lines -1
        
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
        split_fname = this_fname.split('.')    
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
 
        # set the name of the output figure
        OutputFigureName = DataDirectory+fname_prefix
        print "The output figure name is: "
        print OutputFigureName
   
        # Okay, now we move on to extract the data. 
        # first get the m/n values
        line = lines[0].strip().split(" ")
        
        n_movern = len(line)-1
        
        movern = np.zeros(n_movern)
        collinear_AICc_mean = np.zeros(n_movern)
        collinear_AICc_stdd = np.zeros(n_movern)
        
        for i in range(0,n_movern):
            movern[i] = float(line[i+1])
            
        # now get the AICc colinearity means
        line = lines[1].strip().split(" ")
        
        for i in range(0,n_movern):
            collinear_AICc_mean[i] = float(line[i+1])  
      
        # now get the AICc colinearity standard deviations
        line = lines[2].strip().split(" ")
        
        for i in range(0,n_movern):
            collinear_AICc_stdd[i] = float(line[i+1])
            
        minimum_collinear = min(collinear_AICc_mean)
      
        minimum_element = 0
        for i in range(0,n_movern):
            if collinear_AICc_mean[i] == minimum_collinear:
                minimum_element = i
                        
        minimum_mn_collinear = movern[minimum_element]
        minimum_stdd_AIC_collinear = collinear_AICc_stdd[minimum_element]
        
        print 'min AICc: '+str(minimum_collinear)+ ' and min stdd: ' +str(minimum_stdd_AIC_collinear)   
          
        line_for_variability = np.zeros(n_movern)
        for i in range(0,n_movern):
            line_for_variability[i] =  minimum_stdd_AIC_collinear+minimum_collinear
            
        
        # now figure out how many channels there are
        n_channels = (n_lines-3)/2
        
        AICc_mean_vecvec = np.zeros((n_channels,n_movern))
        AICc_stdd_vecvec = np.zeros((n_channels,n_movern))
         
        print 'n_channels: '+ str(n_channels)
        
        # loop through the channels
        for i in range(0,n_channels):
            line = lines[2*i+3].strip().split(" ")
            line2 = lines[2*i+4].strip().split(" ")
            
            for j in range(0,n_movern):
                AICc_mean_vecvec[i,j] = float(line[j+1])  
                AICc_stdd_vecvec[i,j] = float(line2[j+1])
                
        #print AICc_mean_vecvec   
        #print AICc_stdd_vecvec 
        
        #print AICc_mean_vecvec[0]
        #print AICc_mean_vecvec[0,0] 
         #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        # AICc plot for the collinear test
        fig = plt.figure(1, facecolor='white',figsize=(10,7.5))
        ax =  fig.add_subplot(1,1,1)
        
        dmn = movern[1]-movern[0]    
        plt.xlim(movern.min()-dmn*0.5, movern.max()+dmn*0.5)
    
        
    
        #plt.plot(movern, collinear_AICc_mean, movern, collinear_AICc_mean+collinear_AICc_stdd, movern, collinear_AICc_mean-collinear_AICc_stdd, "-", linewidth=2) 
        plt.errorbar(movern, collinear_AICc_mean, collinear_AICc_stdd,color='blue',lw=3 )
        plt.plot(movern, line_for_variability,'k',lw=3, dashes=(10,2), alpha = 0.5)
        plt.fill_between(movern, collinear_AICc_mean+collinear_AICc_stdd,  collinear_AICc_mean-collinear_AICc_stdd, color='grey', alpha=0.3,lw=3)
        
        plt.annotate('Minimum $AICc$ at $m/n$ = '+str(minimum_mn_collinear), xy=(minimum_mn_collinear, minimum_collinear+minimum_stdd_AIC_collinear),  xycoords='data',
                    xytext=(0.8, 0.90), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle="fancy", #linestyle="dashed",
                                color="0.5",
                                shrinkB=5,
                                connectionstyle="arc3,rad=0.1",
                                ),
                    horizontalalignment='right', verticalalignment='top',fontsize = label_size)   
     
        #plt.text((movern[0]+movern[0])/2, minimum_collinear,'Minimum $AICc$ plus 1 standard deviation',verticalalignment='top',color='black')  
       
        # Configure final plot
        ax.spines['top'].set_linewidth(2.5)
        ax.spines['left'].set_linewidth(2.5)
        ax.spines['right'].set_linewidth(2.5)
        ax.spines['bottom'].set_linewidth(2.5) 
        ax.tick_params(axis='both', width=2.5)    
    
    
        plt.xlabel('$m/n$',fontsize = axis_size)
        plt.ylabel('$AICc$',fontsize = axis_size)
        plt.title('Results, collinearity test, '+ param_str,fontsize = title_size)
        plt.savefig(OutputFigureName+'_cl.' + OutputFigureFormat, format=OutputFigureFormat)
        plt.clf()
    
         #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        # AICc plot for the individual channels
        for i in range(0,n_channels):
            thisfig = plt.figure(i+2, facecolor='white',figsize=(10,7.5))
            ax = thisfig.add_subplot(1,1,1)
            
            plt.errorbar(movern, AICc_mean_vecvec[i], AICc_stdd_vecvec[i],lw=3, color='blue')
            
            plt.xlim(movern.min()-dmn*0.5, movern.max()+dmn*0.5)
                    
            thischanAICc = AICc_mean_vecvec[i]
            #print i
            #print thischanAICc
            minimum_chan = min(AICc_mean_vecvec[i])
     
            minimum_element = 0
            for j in range(0,n_movern):
                if AICc_mean_vecvec[i,j] == minimum_chan:
                    minimum_element = j
                        
            minimum_mn_chan = movern[minimum_element]
            minimum_stdd_AIC_chan = AICc_stdd_vecvec[i,minimum_element]
     
            for j in range(0,n_movern):
               line_for_variability[j] =  minimum_stdd_AIC_chan+minimum_chan 
        
            plt.plot(movern, line_for_variability,'k',lw=3, dashes=(10,2), alpha = 0.5)
            plt.fill_between(movern, AICc_mean_vecvec[i]+AICc_stdd_vecvec[i],  AICc_mean_vecvec[i]-AICc_stdd_vecvec[i], lw=3,color='grey', alpha=0.3)
        
            plt.annotate('Minimum $AICc$ at $m/n$ = '+str(minimum_mn_chan), xy=(minimum_mn_chan, minimum_chan+minimum_stdd_AIC_chan),  xycoords='data',
                    xytext=(0.8, 0.90), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle="fancy", #linestyle="dashed",
                                color="0.5",
                                shrinkB=5,
                                connectionstyle="arc3,rad=0.1",
                                ),
                    horizontalalignment='right', verticalalignment='top',fontsize = label_size)   
     
            #plt.text((movern[0]+movern[0])/2, minimum_chan,'Minimum $AICc$ plus 1 standard deviation',verticalalignment='top',color='black',fontsize = label_size)  
       
            # Configure final plot
            plt.xlabel('$m/n$',fontsize = axis_size)
            plt.ylabel('$AICc$',fontsize = axis_size)
            plt.title('Results, channel '+str(i)+ ', '+param_str ,fontsize = title_size)
            
            ax.spines['top'].set_linewidth(3)
            ax.spines['left'].set_linewidth(3)
            ax.spines['right'].set_linewidth(3)
            ax.spines['bottom'].set_linewidth(3) 
            ax.tick_params(axis='both', width=3) 
            plt.savefig(OutputFigureName+'_'+str(i) + '.' + OutputFigureFormat, format=OutputFigureFormat)
            plt.clf()

    
    
if __name__ == "__main__":
    make_plots()