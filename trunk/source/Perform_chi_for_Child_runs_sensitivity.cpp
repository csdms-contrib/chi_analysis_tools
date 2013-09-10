//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Perform_chi_analysis_for_Child_runs.cpp
//
// This program crunches the data from child runs. It is different from
// other Perform_chi_analysis_v3.cpp in that it has no reference to
// LSDRaster objects (since the data is based on a TIN)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2013 Simon M. Mudd 2013
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include "LSDStatsTools.hpp"
#include "LSDChiNetwork.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}

	string DEM_name;
	string fill_ext = "_fill";
	file_info_in >> DEM_name;
	int junction_number;
	double pruning_threshold;
	int threshold;
	double A_0;
	int minimum_segment_length;
	double sigma;
	double start_movern;
	double d_movern;
	double Minimum_Slope;
	int n_movern;
	int target_nodes;
	int n_iterations;
	double fraction_dchi_for_variation;
	double vertical_interval;
	double horizontal_interval;
	double area_thin_frac;
	int target_skip;

	file_info_in >> Minimum_Slope >> threshold >> junction_number
				>> pruning_threshold >> A_0 >> minimum_segment_length >> sigma >> start_movern
				>> d_movern >> n_movern >> target_nodes >> n_iterations >> fraction_dchi_for_variation
				>> vertical_interval >> horizontal_interval >> area_thin_frac >> target_skip;


	cout << "Paramters of this run: " << endl
	     << "junction number: " << junction_number << endl
	     << "pruning_threshold: " << pruning_threshold << endl
	     << "threshold: " << threshold << endl
	     << "A_0: " << A_0 << endl
	     << "minimum_segment_length: " << minimum_segment_length << endl
	     << "sigma: " << sigma << endl
	     << "start_movern " << start_movern << endl
	     << "d_movern: " << d_movern << endl
	     << "n_movern: " << n_movern << endl
	     << "target_nodes: " << target_nodes << endl
	     << "n_iterarions: " << n_iterations << endl
	     << "fraction_dchi_for_variation: " << fraction_dchi_for_variation << endl
	     << "vertical interval: " << vertical_interval << endl
	     << "horizontal interval: " << horizontal_interval << endl
	     << "area thinning fraction for SA analysis: " << area_thin_frac
       << "target_skip is: " << target_skip << endl;

	file_info_in.close();

	string Chan_ext = ".chan";
	string Chan_for_chi_ingestion_fname = path_name+DEM_name+Chan_ext;
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// create the chi network
	LSDChiNetwork ChiNetwork(Chan_for_chi_ingestion_fname);
	LSDChiNetwork ChiNetwork_extended(Chan_for_chi_ingestion_fname);
	ChiNetwork_extended.extend_tributaries_to_outlet();
	//cout << "N channels: " << ChiNetwork.get_n_channels() << endl;

	// get the best fit m over n ratio for this basin
	double bf_cum_movn_ms, bf_colinear_movn_breaks;
	string bfchi_fname = "BF_mover_n_";
	string fpt_ext = ".movern";

	vector<double> m_over_n_values;
	vector<double> AICc_mean_breaks;
	vector<double> AICc_stdd_breaks;
	vector< vector<double> > AICc_vals;
	vector< vector<double> > AICc_stddev;

  // initiate parameters for sensitivity analysis
  int n_sigma = 5;
  double d_sigma = 2.0;
  double start_sigma = 1.0;
  int n_skip = 3;
  int d_skip = 1;
  int start_skip = 2;
  int n_minimum_segment_length = 2;
  int d_minimum_segment_length = 10;
  int start_minimum_segment_length = 10;
  int n_target_nodes = 4;
  int d_target_nodes = 10;
  int start_target_nodes = 80;
  string sigma_str;
  string skip_str;
  string msl_str;
  string tn_str;
  string space = "_";
  string param_str;

  for (int sig = 0; sig<n_sigma; sig++)
  {
    for (int skp = 0; skp <n_skip; skp++)
    {
      for (int msl = 0; msl<n_minimum_segment_length; msl++)
      {
        for (int tn = 0; tn<n_target_nodes; tn++)
        {
          sigma = double(sig)*d_sigma+start_sigma;
          target_skip = skp*d_skip+start_skip;
          minimum_segment_length = msl*d_minimum_segment_length+start_minimum_segment_length;
          target_nodes = tn*d_target_nodes+start_target_nodes;

          cout << endl << endl << "now parameters: " <<endl;
          cout << "sigma: " << sigma << endl;
          cout << "target_skip: " << target_skip << endl;
          cout << "minimum_segment_length: " << minimum_segment_length << endl;
          cout << "target_nodes: " << target_nodes << endl;

          // convert the parameters to a string for the output filename
		      sigma_str = static_cast<ostringstream*>( &(ostringstream() << sigma) )->str();
          skip_str = static_cast<ostringstream*>( &(ostringstream() << target_skip) )->str();
          msl_str = static_cast<ostringstream*>( &(ostringstream() << minimum_segment_length) )->str();
          tn_str = static_cast<ostringstream*>( &(ostringstream() << target_nodes) )->str();

          param_str = sigma_str+space+skip_str+space+msl_str+space+tn_str;

          // get the best fit m over n in two different ways
        	int Monte_Carlo_switch = 1;
        	bf_colinear_movn_breaks = ChiNetwork.search_for_best_fit_m_over_n_colinearity_test_with_breaks(A_0,
        						 n_movern, d_movern,start_movern,
        						 minimum_segment_length, sigma, target_skip, target_nodes, n_iterations,
        						 m_over_n_values, AICc_mean_breaks, AICc_stdd_breaks,
        						 Monte_Carlo_switch);

        	bf_cum_movn_ms = ChiNetwork_extended.search_for_best_fit_m_over_n_individual_channels_with_breaks_monte_carlo(A_0,
        	                     n_movern, d_movern,start_movern,
        						 minimum_segment_length, sigma, target_skip, target_nodes, n_iterations,
        						 m_over_n_values, AICc_vals, AICc_stddev);

        	// now print the results from the m_over_n analysis
        	ofstream mn_analysis_out;
        	string mn_fo_fname = (path_name+DEM_name+bfchi_fname+param_str+fpt_ext);
        	mn_analysis_out.open(mn_fo_fname.c_str());
        	mn_analysis_out << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << m_over_n_values[i] << " ";
        	}
        	mn_analysis_out << endl << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << AICc_mean_breaks[i] << " ";
        	}
        	mn_analysis_out << endl << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << AICc_stdd_breaks[i] << " ";
        	}
        	mn_analysis_out << endl;



        	for(int chan = 0; chan<  ChiNetwork.get_n_channels(); chan++)
        	{
        		mn_analysis_out << chan << " ";
        		for (int i = 0; i< n_movern; i++)
        		{
        			mn_analysis_out << AICc_vals[chan][i] << " ";
        		}
        		mn_analysis_out << endl;
        		mn_analysis_out << chan << " ";
        		for (int i = 0; i< n_movern; i++)
        		{
        			mn_analysis_out << AICc_stddev[chan][i] << " ";
        		}
        		mn_analysis_out << endl;
        	}
        	mn_analysis_out.close();
        	cout << "best fit m over n mainstem: " << bf_cum_movn_ms
        	     << " and collinear: " << bf_colinear_movn_breaks << endl;

        	// now do monte carlo sampling to get the different steepnesses of the channels
        	fpt_ext = ".tree";

        	// get the breaks of all the channels
        	ChiNetwork_extended.split_all_channels(A_0, bf_colinear_movn_breaks, n_iterations,
        						target_skip, target_nodes, minimum_segment_length, sigma);

        	// monte carlo sample all channels
        	ChiNetwork_extended.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, bf_colinear_movn_breaks, n_iterations,
        							target_skip, minimum_segment_length, sigma);

        	string fpt_mc = "_fullProfileMC_colinear_";

        	ChiNetwork_extended.print_channel_details_to_file_full_fitted((path_name+DEM_name+fpt_mc+param_str+fpt_ext));

        	//=-=-=-=-=-=-
        	// NOW do the best fit m/n of the mainstem
        	// get the breaks of all the channels
        	ChiNetwork_extended.split_all_channels(A_0, bf_cum_movn_ms, n_iterations,
        						target_skip, target_nodes, minimum_segment_length, sigma);


        	// monte carlo sample all channels
        	ChiNetwork_extended.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, bf_cum_movn_ms, n_iterations,
        							target_skip, minimum_segment_length, sigma);

        	fpt_mc = "_fullProfileMC_mainstem_";

        	ChiNetwork_extended.print_channel_details_to_file_full_fitted((path_name+DEM_name+fpt_mc+param_str+fpt_ext));
        	}
        }
      }
    }




}
