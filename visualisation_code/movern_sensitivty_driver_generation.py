## movern_sensitivty_driver_generation.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This function takes a single driver file and then spawns
##  many driver files with different mean skip, minimum segement lengths, 
## , total nodes and sigma values in order to perform a sensitivy analysis
##  of the best fit m over n ratio as a batch run 
## This allows one to run the code on many CPUs, saving some time as 
## opposed to running the entire analysis in serial
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 04/09/2013
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def movern_sensitivty_driver_generation():
    
    # these are the number of different parameter values you want to use
    n_skip = 2
    n_sigma = 1
    n_msl = 2
    n_tn = 2

    # this is the starting value of the parameter values    
    start_skip = 1
    start_sigma = 3.0
    start_msl = 10
    start_tn = 80
 
    # these are the change to the parameter value each time you iterate   
    d_skip = 1
    d_sigma = 3
    d_msl = 5
    d_tn = 10

    #########################
    #                       #
    #   READ IN THE DATA    #
    #                       #
    #########################

    # det the directory and filename
    DataDirectory =  'm:\\papers\\Segment_fitting\\Data_and_code_repository\\Data\\Apennines\\test_sensitivity\\' 
    DriverFileName = 'rio_torto.driver'
  
    # combine these
    FileName = DataDirectory+DriverFileName
    
    print "FileName is"+FileName
  
  
    # first, get the prefix of the file
    split_fname = DriverFileName.split('.')
    no_tree_levs = len(split_fname) 
    fname_prefix  = split_fname[0]  
    if (no_tree_levs > 2):
        for i in range (1,no_tree_levs-1):
            fname_prefix+= "."+split_fname[i]
    
    print "The file prefix is: " + fname_prefix

    f = open(FileName,'r')  # open file
    lines = f.readlines()   # read in the data
    f.close()
    
    # overwrite the lines for sigma, msl, target nodes and mean skip
    counter = 0
    for sk in range(0,n_skip):
        for sg in range(0,n_sigma):
            for minsl in range(0,n_msl):
                for tarn in range(0,n_tn):
                    
                    counter+=1
                    
                    skip = start_skip+sk*d_skip
                    sigma = start_sigma+sg*d_sigma
                    msl = start_msl+minsl*d_msl
                    tn = start_tn+tarn*d_tn
                    
                    print "skip: " + str(skip)+ " sigma: " + str(sigma)+ " msl: " + str(msl) + " tn: " + str(tn)
                    
                    lines[6] = str(msl)+'\n'
                    lines[7] = str(sigma)+'\n'
                    lines[11]= str(tn)+'\n'
                    lines[17]= str(skip)+'\n'
                    
                    print lines
                    
                    this_fname = fname_prefix+"."+str(counter)+".driver"
                    f = open(DataDirectory + this_fname,'w')  # open file
                    f.writelines(lines)
                    f.close()
                    
                    
   

if __name__ == "__main__":
    movern_sensitivty_driver_generation()
    