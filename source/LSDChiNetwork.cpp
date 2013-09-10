//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiNetwork
// Land Surface Dynamics ChiNetwork
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for analysing channels using the integral method of channel
//  analysis
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


// source code for the LSDChiNetwork object
// this object organizes the chi analysis for several
// channel segments. It can be used in conjunction with the LSD topographic
// analysis package and also used as a standalone package
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		5/4/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <fstream>
#include "TNT/tnt.h"
#include "LSDChiNetwork.hpp"
#include "LSDMostLikelyPartitionsFinder.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChiNetwork_CPP
#define LSDChiNetwork_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChiNetwork::create(string channel_network_fname)
{
	ifstream channel_data_in;
	channel_data_in.open(channel_network_fname.c_str());

	if( channel_data_in.fail() )
	{
		cout << "\nFATAL ERROR: the channel network file \"" << channel_network_fname
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}

	int channel_number;
	int receiver_cnumber;
	int recevier_cnode;

	int node;
	int row;
	int col;

	double flow_dist;
	double elev;
	double drain_area;

	int last_cn = 0;		// this is 1 if this is the first node in a channel
	int last_receiver_node = -1;
	int last_receiver_channel = -1;

	vector<int> empty_int;
	vector<double> empty_double;

	vector<int> node_vec;
	vector<int> row_vec;
	vector<int> col_vec;
	vector<double> flow_dist_vec;
	vector<double> elev_vec;
	vector<double> drain_area_vec;

	channel_data_in >> NRows >> NCols >> XMinimum >> YMinimum >> DataResolution >> NoDataValue;

	while( channel_data_in >> channel_number >> receiver_cnumber >> recevier_cnode
	                       >> node >> row >> col >> flow_dist >> elev >> drain_area)
	{

		// get the receiver_channel and receiver node for the first channel (these will be recursive)
		if (last_receiver_node == -1)
		{
			last_receiver_node = recevier_cnode;
			last_receiver_channel = receiver_cnumber;
		}

		// if this is a new channel add the old channel data to the data members and reset the
		// vectors for assimilating data
		if (channel_number != last_cn)
		{

			cout << "new channel: " << channel_number << " last cn: " << last_cn << endl;

			node_indices.push_back(node_vec);
			row_indices.push_back(row_vec);
			col_indices.push_back(col_vec);

			elevations.push_back(elev_vec);
			flow_distances.push_back(flow_dist_vec);
			drainage_areas.push_back(drain_area_vec);

			node_on_receiver_channel.push_back(last_receiver_node);
			receiver_channel.push_back(last_receiver_channel);

			node_vec = empty_int;
			row_vec = empty_int;
			col_vec = empty_int;
			flow_dist_vec = empty_double;
			elev_vec = empty_double;
			drain_area_vec = empty_double;

			// reset the receiver nodde and channel
			last_receiver_node = recevier_cnode;
			last_receiver_channel = receiver_cnumber;

			last_cn = channel_number;
		}

		// now push back the data
		node_vec.push_back(node);
		row_vec.push_back(row);
		col_vec.push_back(col);
		flow_dist_vec.push_back(flow_dist);
		elev_vec.push_back(elev);
		drain_area_vec.push_back(drain_area);
	}


	// push back the data for the final channel
	node_indices.push_back(node_vec);
	row_indices.push_back(row_vec);
	col_indices.push_back(col_vec);
	elevations.push_back(elev_vec);
	flow_distances.push_back(flow_dist_vec);
	drainage_areas.push_back(drain_area_vec);
	node_on_receiver_channel.push_back(last_receiver_node);
	receiver_channel.push_back(last_receiver_channel);

	// now initiate the chi values
	int n_channels = int(elevations.size());
	for (int i = 0; i< n_channels; i++)
	{
		int n_nodes_in_channel = (node_indices[i].size());
		vector<double> empty_chi(n_nodes_in_channel,0.0);
		chis.push_back(empty_chi);
	}

	// close the infile
	channel_data_in.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// this function extends the tributaries so all tributaries start from their source and then
// run all the way down to the outlet. That is, downstream nodes are reapeated in each
// tributary.
//
// The purpose of this is to pick up segments that may only have propagated
// a very short distance upstream from the tributary junction and are therefore
// shorter than the minimum segment length.
//
// The function is constucted so all the member functions should continue to work on the channel
// network once the tributaries have been recalucalted
//
// IMPORTANT: This only works if all the tributaries drain to the mainstem
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChiNetwork::extend_tributaries_to_outlet()
{
	int n_channels = elevations.size();
	int n_nodes_on_mainstem = elevations[0].size();

	if(n_channels>1)
	{
		for(int chan = 1; chan<n_channels; chan++)
		{
			int this_receiver_channel = receiver_channel[chan];
			int this_node_on_receiver_channel = node_on_receiver_channel[chan];
			for(int node = this_node_on_receiver_channel+1; node< n_nodes_on_mainstem; node++)
			{
				node_indices[chan].push_back( node_indices[this_receiver_channel][node] );
				row_indices[chan].push_back ( row_indices[this_receiver_channel][node] );
				col_indices[chan].push_back ( col_indices[this_receiver_channel][node] );
				elevations[chan].push_back ( elevations[this_receiver_channel][node] );
				flow_distances[chan].push_back ( flow_distances[this_receiver_channel][node] );
				drainage_areas[chan].push_back ( drainage_areas[this_receiver_channel][node] );
				chis[chan].push_back ( chis[this_receiver_channel][node] );
			}

			node_on_receiver_channel[chan] = int( node_indices[chan].size() )-1;
			receiver_channel[chan] = chan;
		}
	}
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function prints the details of an individual channel to screen for bug checking
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChiNetwork::print_channel_details_to_screen(int channel_number)
{
  if (channel_number >= int(elevations.size()))
    {
      cout << "Channel number is in reference to a channel that doesn't exist" << endl;
      cout << " using the last possible channel " << endl;
      channel_number = int( elevations.size())-1;
    }

	vector<int> node = node_indices[channel_number];
	vector<int> row = row_indices[channel_number];
	vector<int> col = col_indices[channel_number];
	vector<double> elevation = elevations[channel_number];
	vector<double> flow_distance = flow_distances[channel_number];
	vector<double> drainage_area = drainage_areas[channel_number];
	vector<double> chi = chis[channel_number];

	int n_nodes = node.size();
	for (int i = 0; i< n_nodes; i++)
	  {
	    cout << channel_number << " " << receiver_channel[channel_number] << " "
                 << node_on_receiver_channel[channel_number] << " "
                 << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
		         << chi[i] << " " << elevation[i] << " " << drainage_area[i] << endl;
	  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function prints the details of all channels to a file
// format is
// A_0 m_over_n
// channel_number node_on_receiver_channel node_index row col flow_distance chi elevation darainage_area
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::print_channel_details_to_file(string fname, double A_0, double m_over_n)
{

	ofstream channel_profile_out;
	channel_profile_out.open(fname.c_str());
	calculate_chi(A_0, m_over_n);

	channel_profile_out << A_0 << " " << m_over_n << endl;
	int n_channels = chis.size();
	for (int channel_number = 0; channel_number< n_channels; channel_number++)
	{
		vector<int> node = node_indices[channel_number];
		vector<int> row = row_indices[channel_number];
		vector<int> col = col_indices[channel_number];
		vector<double> elevation = elevations[channel_number];
		vector<double> flow_distance = flow_distances[channel_number];
		vector<double> drainage_area = drainage_areas[channel_number];
		vector<double> chi = chis[channel_number];

		int n_nodes = node.size();
		for (int i = 0; i< n_nodes; i++)
		{
		    channel_profile_out << channel_number << " " << receiver_channel[channel_number] << " "
    	             << node_on_receiver_channel[channel_number] << " "
    	             << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
			         << chi[i] << " " << elevation[i] << " " << drainage_area[i] << endl;
		}
	}
	channel_profile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function prints the details of all channels to a file
// it includes data from monte carlo fittin
// format is
// A_0 m_over_n
// channel_number node_on_receiver_channel node_index row col flow_distance chi elevation darainage_area...
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::print_channel_details_to_file_full_fitted(string fname)
{

	double m_over_n = m_over_n_for_fitted_data;
	double A_0 = A_0_for_fitted_data;

	ofstream channel_profile_out;
	channel_profile_out.open(fname.c_str());
	calculate_chi(A_0, m_over_n);

	channel_profile_out << A_0 << " " << m_over_n << endl;
	int n_channels = chis.size();

	if (chi_b_means.size() == chis.size())
	{

		for (int channel_number = 0; channel_number< n_channels; channel_number++)
		{
			vector<int> node = node_indices[channel_number];
			vector<int> row = row_indices[channel_number];
			vector<int> col = col_indices[channel_number];
			vector<double> elevation = elevations[channel_number];
			vector<double> flow_distance = flow_distances[channel_number];
			vector<double> drainage_area = drainage_areas[channel_number];
			vector<double> chi = chis[channel_number];

			vector<double> m_mean = chi_m_means[channel_number];
			vector<double> m_standard_deviation = chi_m_standard_deviations[channel_number];
			vector<double> m_standard_error = chi_m_standard_errors[channel_number];

			vector<double> b_mean = chi_b_means[channel_number];
			vector<double> b_standard_deviation = chi_b_standard_deviations[channel_number];
			vector<double> b_standard_error = chi_b_standard_errors[channel_number];

			vector<double> DW_mean = chi_DW_means[channel_number];
			vector<double> DW_standard_deviation = chi_DW_standard_deviations[channel_number];
			vector<double> DW_standard_error = chi_DW_standard_errors[channel_number];

			vector<double> fitted_elev_mean = all_fitted_elev_means[channel_number];
			vector<double> fitted_elev_standard_deviation = all_fitted_elev_standard_deviations[channel_number];
			vector<double> fitted_elev_standard_error = all_fitted_elev_standard_errors[channel_number];

			vector<int> n_data_points_uic = n_data_points_used_in_stats[channel_number];

			int n_nodes = node.size();
			for (int i = 0; i< n_nodes; i++)
			{
				channel_profile_out << channel_number << " " << receiver_channel[channel_number] << " "
						 << node_on_receiver_channel[channel_number] << " "
						 << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
						 << chi[i] << " " << elevation[i] << " " << drainage_area[i] << " "
						 << n_data_points_uic[i] << " "
						 << m_mean[i] << " " << m_standard_deviation[i] << " " << m_standard_error[i] << " "
						 << b_mean[i] << " " << b_standard_deviation[i] << " " << b_standard_error[i] << " "
						 << DW_mean[i] << " " << DW_standard_deviation[i] << " " << DW_standard_error[i] << " "
						 << fitted_elev_mean[i] << " " << fitted_elev_standard_deviation[i] << " "
						 << fitted_elev_standard_error[i] << " " <<  endl;
			}
		}
	}
	else
	{
		cout << "LSDChiNetwork Line 276 you don't seem to have run the monte carlo fitting routine" << endl;
	}
	channel_profile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function prints the details of all channels to a file
// it includes data from monte carlo fittin
// format is
// A_0 m_over_n
// channel_number node_on_receiver_channel node_index row col flow_distance chi elevation darainage_area...
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::print_channel_details_to_file_full_fitted(string fname, int target_nodes, int minimum_segment_length)
{

	double m_over_n = m_over_n_for_fitted_data;
	double A_0 = A_0_for_fitted_data;

	ofstream channel_profile_out;
	channel_profile_out.open(fname.c_str());
	calculate_chi(A_0, m_over_n);

	channel_profile_out << A_0 << " " << m_over_n << endl;
	int n_channels = chis.size();

	if (chi_b_means.size() == chis.size())
	{
		int N = calculate_skip(target_nodes);
		is_channel_long_enough_test(minimum_segment_length,N);

		for (int channel_number = 0; channel_number< n_channels; channel_number++)
		{
			if(is_tributary_long_enough[channel_number] == 1)
			{
				vector<int> node = node_indices[channel_number];
				vector<int> row = row_indices[channel_number];
				vector<int> col = col_indices[channel_number];
				vector<double> elevation = elevations[channel_number];
				vector<double> flow_distance = flow_distances[channel_number];
				vector<double> drainage_area = drainage_areas[channel_number];
				vector<double> chi = chis[channel_number];

				vector<double> m_mean = chi_m_means[channel_number];
				vector<double> m_standard_deviation = chi_m_standard_deviations[channel_number];
				vector<double> m_standard_error = chi_m_standard_errors[channel_number];

				vector<double> b_mean = chi_b_means[channel_number];
				vector<double> b_standard_deviation = chi_b_standard_deviations[channel_number];
				vector<double> b_standard_error = chi_b_standard_errors[channel_number];

				vector<double> DW_mean = chi_DW_means[channel_number];
				vector<double> DW_standard_deviation = chi_DW_standard_deviations[channel_number];
				vector<double> DW_standard_error = chi_DW_standard_errors[channel_number];

				vector<double> fitted_elev_mean = all_fitted_elev_means[channel_number];
				vector<double> fitted_elev_standard_deviation = all_fitted_elev_standard_deviations[channel_number];
				vector<double> fitted_elev_standard_error = all_fitted_elev_standard_errors[channel_number];

				vector<int> n_data_points_uic = n_data_points_used_in_stats[channel_number];

				int n_nodes = node.size();
				for (int i = 0; i< n_nodes; i++)
				{
					channel_profile_out << channel_number << " " << receiver_channel[channel_number] << " "
							 << node_on_receiver_channel[channel_number] << " "
							 << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
							 << chi[i] << " " << elevation[i] << " " << drainage_area[i] << " "
							 << n_data_points_uic[i] << " "
							 << m_mean[i] << " " << m_standard_deviation[i] << " " << m_standard_error[i] << " "
							 << b_mean[i] << " " << b_standard_deviation[i] << " " << b_standard_error[i] << " "
							 << DW_mean[i] << " " << DW_standard_deviation[i] << " " << DW_standard_error[i] << " "
							 << fitted_elev_mean[i] << " " << fitted_elev_standard_deviation[i] << " "
							 << fitted_elev_standard_error[i] << " " <<  endl;
				}
			}
		}
	}
	else
	{
		cout << "LSDChiNetwork Line 276 you don't seem to have run the monte carlo fitting routine" << endl;
	}
	channel_profile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDChiNetwork::slope_area_extraction_vertical_intervals
//
// this one is the vertical intervals version: it measures slope over fixed vertical
// intervals as reccomended by Wobus et al 2006.
// This function gets slope data and area data for use in making slope area plots.
//  It generates several data elements, which are written to the file with name fname (passed to
//  function). The file format is for each row
//   SA_file << chan << " " << start_row << " " << mp_row << " " << end_row << " "
//		<< start_col << " " << mp_col << " " << end_row << " "
//		<< start_interval_elevations << " "
//		<< mp_interval_elevations << " " << end_interval_elevations << " "
//		<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
//		<< end_interval_flowdistance << " "
//		<< start_area << " " << mp_area << " " << end_area << " " << slope
//		<< " " << log10(mp_area) << " " << log10(slope) << endl;
//
//  where start, mp and end denote the start of the interval over which slope is measured, the midpoint
//    and the end.
//
// The area thin fraction is used to thin the data so that segments with large changes in drainage
// area are not used in the regression (because these will affect the mean slope)
// the fraction is determined by (downslope_area-upslope_area)/midpoint_area.
// So if the fraction is 1 it means that the change is area is equal to the area at the midpoint
// a restictive value is 0.05, you will eliminate major tributaries with a 0.2, and
// 1 will catch almost all of the data.
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::slope_area_extraction_vertical_intervals(double interval, double area_thin_fraction,
															string fname)
{
	int n_channels = elevations.size();
	double half_interval = interval/2;
	double target_end_interval_elevations;
	double target_mp_interval_elevations;
	double start_interval_elevations;
	double end_interval_elevations;
	double mp_interval_elevations = -9999;
	double start_interval_flowdistance;
	double end_interval_flowdistance;
	double mp_interval_flowdistance = -9999;
	int start_row, start_col;
	int mp_row = -9999;
	int mp_col = -9999;
	int end_row, end_col;
	double start_area;
	double mp_area = -9999;
	double end_area;
	double area_thin_frac_for_test;

	cout << "n channels is: " << n_channels << endl;

	ofstream SA_file;
	SA_file.open(fname.c_str());

	// loop through the channels
	for (int chan = 0; chan<n_channels; chan++)
	{
		int n_nodes_this_channel = (elevations[chan].size());
		int final_node_flag = 0;
		int end_flag;
		int mp_flag;
		int n = 0;

		//cout << "n: " << n << " final_node_flag: " << final_node_flag << " nnodes: " << n_nodes_this_channel << endl;
		// now, start at the top node (node 0) of each channel and
		// work down. The algorithm looks downstream until it hits
		// the midpoint, and then continues until it hits
		// the end point
		// if it encounters the end of the channel before it hits
		// the end point the loop exits
		while (final_node_flag == 0 && n < n_nodes_this_channel)
		{
			start_interval_elevations = elevations[chan][n];
			start_interval_flowdistance = flow_distances[chan][n];
			start_area = drainage_areas[chan][n];
			start_row = row_indices[chan][n];
			start_col = col_indices[chan][n];

			target_end_interval_elevations = 	start_interval_elevations-interval;
			target_mp_interval_elevations = start_interval_elevations-half_interval;

			int search_node = n+1;

			// reset midpoint and end flags
			mp_flag = 0;
			end_flag = 0;
			//cout << "n is: " << n << " and end_flag is: " << end_flag << endl;

			// now work downstream
			while (search_node < n_nodes_this_channel && end_flag == 0)
			{
				//cout << "search_node: " << search_node << " elev: " << elevations[chan][search_node]
				//     << " and target mp, end: " << target_mp_interval_elevations << " " << target_end_interval_elevations << endl;

				// see if search node is the midpoint node
				if ( elevations[chan][search_node] <= target_mp_interval_elevations && mp_flag  == 0)
				{
					mp_interval_elevations = elevations[chan][search_node];
					mp_interval_flowdistance = flow_distances[chan][search_node];
					mp_area = drainage_areas[chan][search_node];
					mp_row = row_indices[chan][search_node];
					mp_col = col_indices[chan][search_node];

					// set midpoint flag so it doens't collect downstream nodes
					mp_flag = 1;
				}

				// see if the search node is the end node
				if (elevations[chan][search_node] <= target_end_interval_elevations)
				{
					end_interval_elevations = elevations[chan][search_node];
					end_interval_flowdistance = flow_distances[chan][search_node];
					end_area = drainage_areas[chan][search_node];
					end_row = row_indices[chan][search_node];
					end_col = col_indices[chan][search_node];

					// set end flag so the search node is reset
					end_flag = 1;
				}

				search_node++;
			}

			// if the end flag == 0 (that means it found the end interval) then print the information
			// to file. If it didn't reach the end flag that means the previous node was the final
			// flag
			if (end_flag == 1)
			{
				double slope = (start_interval_elevations-end_interval_elevations)/
				               (start_interval_flowdistance-end_interval_flowdistance);

				// the data take a log of the slope so it is necessary to have this statement
				// the the case of a negative or zero slope
				if(slope <=0)
				{
					slope = 0.0000000001;
				}

				area_thin_frac_for_test = (end_area-start_area)/mp_area;
				if (area_thin_frac_for_test < area_thin_fraction)
				{
					SA_file << chan << " " << start_row << " " << mp_row << " " << end_row << " "
					    << start_col << " " << mp_col << " " << end_col << " "
					    << start_interval_elevations << " "
				     	<< mp_interval_elevations << " " << end_interval_elevations << " "
				     	<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
				     	<< end_interval_flowdistance << " "
				     	<< start_area << " " << mp_area << " " << end_area << " " << slope
				     	<< " " << log10(mp_area) << " " << log10(slope) << endl;
				}
			}
			else
			{
				final_node_flag = 1;
			}

			n++;

		}
	}
	SA_file.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDChiNetwork::slope_area_extraction_horizontal_intervals
//
// this one is the horizontal intervals version: it measures slope over fixed flow distance
// as used by many authors including DiBiasie et al 2010 and Ouimet et al 2009
// This function gets slope data and area data for use in making slope area plots.
//  It generates several data elements, which are written to the file with name fname (passed to
//  function). The file format is for each row
//   SA_file << chan << " " << start_row << " " << mp_row << " " << end_row << " "
//		<< start_col << " " << mp_col << " " << end_row << " "
//		<< start_interval_elevations << " "
//		<< mp_interval_elevations << " " << end_interval_elevations << " "
//		<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
//		<< end_interval_flowdistance << " "
//		<< start_area << " " << mp_area << " " << end_area << " " << slope
//		<< " " << log10(mp_area) << " " << log10(slope) << endl;
//
//  where start, mp and end denote the start of the interval over which slope is measured, the midpoint
//    and the end.
//
// The area thin fraction is used to thin the data so that segments with large changes in drainage
// area are not used in the regression (because these will affect the mean slope)
// the fraction is determined by (downslope_area-upslope_area)/midpoint_area.
// So if the fraction is 1 it means that the change is area is equal to the area at the midpoint
// a restictive value is 0.05, you will eliminate major tributaries with a 0.2, and
// 1 will catch almost all of the data.
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::slope_area_extraction_horizontal_intervals(double interval, double area_thin_fraction,
															string fname)
{
	int n_channels = elevations.size();
	double half_interval = interval/2;
	double target_end_interval_flowdistance;
	double target_mp_interval_flowdistance;
	double start_interval_elevations;
	double end_interval_elevations;
	double mp_interval_elevations = -9999;
	double start_interval_flowdistance;
	double end_interval_flowdistance;
	double mp_interval_flowdistance = -9999;
	int start_row, start_col;
	int mp_row = -9999;
	int mp_col = -9999;
	int end_row, end_col;
	double start_area;
	double mp_area = -9999;
	double end_area;
	double area_thin_frac_for_test;

	cout << "n channels is: " << n_channels << endl;

	ofstream SA_file;
	SA_file.open(fname.c_str());

	// loop through the channels
	for (int chan = 0; chan<n_channels; chan++)
	{
		int n_nodes_this_channel = (elevations[chan].size());
		int final_node_flag = 0;
		int end_flag;
		int mp_flag;
		int n = 0;

		//cout << "n: " << n << " final_node_flag: " << final_node_flag << " nnodes: " << n_nodes_this_channel << endl;
		// now, start at the top node (node 0) of each channel and
		// work down. The algorithm looks downstream until it hits
		// the midpoint, and then continues until it hits
		// the end point
		// if it encounters the end of the channel before it hits
		// the end point the loop exits
		while (final_node_flag == 0 && n < n_nodes_this_channel)
		{
			start_interval_elevations = elevations[chan][n];
			start_interval_flowdistance = flow_distances[chan][n];
			start_area = drainage_areas[chan][n];
			start_row = row_indices[chan][n];
			start_col = col_indices[chan][n];

			target_end_interval_flowdistance = 	start_interval_flowdistance-interval;
			target_mp_interval_flowdistance = start_interval_flowdistance-half_interval;

			int search_node = n+1;

			// reset midpoint and end flags
			mp_flag = 0;
			end_flag = 0;
			//cout << "n is: " << n << " and end_flag is: " << end_flag << endl;

			// now work downstream
			while (search_node < n_nodes_this_channel && end_flag == 0)
			{
				//cout << "search_node: " << search_node << " elev: " << elevations[chan][search_node]
				//     << " and target mp, end: " << target_mp_interval_elevations << " " << target_end_interval_elevations << endl;

				// see if search node is the midpoint node
				if ( flow_distances[chan][search_node] <= target_mp_interval_flowdistance && mp_flag  == 0)
				{
					mp_interval_elevations = elevations[chan][search_node];
					mp_interval_flowdistance = flow_distances[chan][search_node];
					mp_area = drainage_areas[chan][search_node];
					mp_row = row_indices[chan][search_node];
					mp_col = col_indices[chan][search_node];

					// set midpoint flag so it doens't collect downstream nodes
					mp_flag = 1;
				}

				// see if the search node is the end node
				if (flow_distances[chan][search_node] <= target_end_interval_flowdistance)
				{
					end_interval_elevations = elevations[chan][search_node];
					end_interval_flowdistance = flow_distances[chan][search_node];
					end_area = drainage_areas[chan][search_node];
					end_row = row_indices[chan][search_node];
					end_col = col_indices[chan][search_node];

					// set end flag so the search node is reset
					end_flag = 1;
				}

				search_node++;
			}

			// if the end flag == 0 (that means it found the end interval) then print the information
			// to file. If it didn't reach the end flag that means the previous node was the final
			// flag
			if (end_flag == 1)
			{
				double slope = (start_interval_elevations-end_interval_elevations)/
				               (start_interval_flowdistance-end_interval_flowdistance);

				// the data take a log of the slope so it is necessary to have this statement
				// the the case of a negative or zero slope
				if(slope <=0)
				{
					slope = 0.0000000001;
				}

				area_thin_frac_for_test = (end_area-start_area)/mp_area;
				if (area_thin_frac_for_test < area_thin_fraction)
				{
					SA_file << chan << " " << start_row << " " << mp_row << " " << end_row << " "
					    << start_col << " " << mp_col << " " << end_col << " "
					    << start_interval_elevations << " "
				     	<< mp_interval_elevations << " " << end_interval_elevations << " "
				     	<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
				     	<< end_interval_flowdistance << " "
				     	<< start_area << " " << mp_area << " " << end_area << " " << slope
				     	<< " " << log10(mp_area) << " " << log10(slope) << endl;
				}
			}
			else
			{
				final_node_flag = 1;
			}

			n++;
		}
	}
	SA_file.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function creates a 2D array that can be used to create a raster of channel properties
// it includes a switch that tells the function what data member to write to the array
// code for data members.
// 1 elevations
// 2  chis
// 3  chi_m_means
// 4  chi_m_standard_deviations
// 5  chi_m_standard_errors
// 6  chi_b_means
// 7  chi_b_standard_deviations
// 8  chi_b_standard_errors
// 9  chi_DW_means
// 10 chi_DW_standard_deviations
// 11 chi_DW_standard_errors
// 12 all_fitted_DW_means
// 13 all_fitted_DW_standard_deviations
// 14 all_fitted_DW_standard_errors
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<double> LSDChiNetwork::data_to_array(int data_member)
{
	// this vecvec will be used to create the array
	vector< vector<double> > data_vecvec;

	// fill the vecvec with the appropriate data based on the switch
	switch ( data_member )
	{
		case 1:
			data_vecvec = elevations;
			break;
		case 2:
			if (chis.size() == elevations.size())
			{
				data_vecvec = chis;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; you need to calculate chi! Raster will have elevations." << endl;
				data_vecvec = elevations;
			}
			break;
		case 3:
			if (chi_m_means.size() == elevations.size())
			{
				data_vecvec = chi_m_means;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 4:
			if (chi_m_standard_deviations.size() == elevations.size())
			{
				data_vecvec = chi_m_standard_deviations;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 5:
			if (chi_m_standard_errors.size() == elevations.size())
			{
				data_vecvec = chi_m_standard_errors;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 6:
			if (chi_b_means.size() == elevations.size())
			{
				data_vecvec = chi_b_means;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 7:
			if (chi_b_standard_deviations.size() == elevations.size())
			{
				data_vecvec = chi_b_standard_deviations;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 8:
			if (chi_b_standard_errors.size() == elevations.size())
			{
				data_vecvec = chi_b_standard_errors;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 9:
			if (chi_DW_means.size() == elevations.size())
			{
				data_vecvec = chi_DW_means;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 10:
			if (chi_DW_standard_deviations.size() == elevations.size())
			{
				data_vecvec = chi_DW_standard_deviations;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 11:
			if (chi_DW_standard_errors.size() == elevations.size())
			{
				data_vecvec = chi_DW_standard_errors;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;

		case 12:
			if (all_fitted_elev_means.size() == elevations.size())
			{
				data_vecvec = all_fitted_elev_means;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 13:
			if (all_fitted_elev_standard_deviations.size() == elevations.size())
			{
				data_vecvec = all_fitted_elev_standard_deviations;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		case 14:
			if (all_fitted_elev_standard_errors.size() == elevations.size())
			{
				data_vecvec = all_fitted_elev_standard_errors;
			}
			else
			{
				cout << "LSDChiNetwork::data_to_array; something is wrong, the channels are the wrong length. Printing elevations" << endl;
				data_vecvec = elevations;
			}
			break;
		default:
			data_vecvec = elevations;
	}
	// now initialize the array with nodata
	Array2D<double> DataArray(NRows,NCols,NoDataValue);
	int row, col;
	double this_data;

	// now loop through channels and nodes filling in the data
	int n_channels = int(data_vecvec.size());
	for (int chan = n_channels-1; chan >= 0; chan--)
	{
		int n_nodes_in_channel = int(data_vecvec[chan].size());
		// see if the algorithm for testing if the channel is long enough
		// has been called
		if( int(is_tributary_long_enough.size()) == n_channels)
		{
			if(is_tributary_long_enough[chan] == 1)
			{
				for (int n = 0; n<n_nodes_in_channel; n++)
				{
					row = row_indices[chan][n];
					col = col_indices[chan][n];
					this_data = data_vecvec[chan][n];

					DataArray[row][col] = this_data;
				}
			}
		}
		else
		{
			for (int n = 0; n<n_nodes_in_channel; n++)
			{
				row = row_indices[chan][n];
				col = col_indices[chan][n];
				this_data = data_vecvec[chan][n];

				DataArray[row][col] = this_data;
			}
		}
	}

	return DataArray;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function calculates the chi values for the channel network using the rectangle rule
// note: the entire network must be caluculated because the chi values of the tributaries
// depend on the chi values of the mainstem
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::calculate_chi(double A_0, double m_over_n)
{
	double dx;				// spacing between nodes

	int n_channels = elevations.size();
	for (int c = 0; c<n_channels; c++)
	{
		// extract the raw data
		vector<double> elevation = elevations[c];
		vector<double> flow_distance = flow_distances[c];
		vector<double> drainage_area = drainage_areas[c];

		// get the number of nodes in channel
		int n_nodes_in_channel = int(elevation.size());

		// initiate a chi vector
		vector<double> chi(elevation.size(),0.0);

		// get the contributing channel and downstream chi
		if (receiver_channel[c] > c)
		{
			cout << "contributing channel has not been calcualted: improper channel ordering" << endl;
			exit(EXIT_FAILURE);
		}

		// if the receiver channel is a downstream channel, get the chi value from the downstream channel
		if (receiver_channel[c] != c)
		{
			// get the chi values from the receiver channel
			vector<double> ds_chi = chis[receiver_channel[c]];

			// set the downstream chi value (which is in the last node on the channel since
			// the data is organized with the furthest upstream first
			chi[n_nodes_in_channel-1] = ds_chi[node_on_receiver_channel[c]];
		}

		// now loop up through the channel, adding chi values
		// note, the channel index are arranges with upstream element first, so you need to go through the channel
		// in reverse order
		for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
		{

			dx = flow_distance[ChIndex]-flow_distance[ChIndex+1];


			chi[ChIndex] = dx*(pow( (A_0/drainage_area[ChIndex] ),
			                    m_over_n))
		                       + chi[ChIndex+1];
			//cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
			//     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
		}	// end loop for this channel

		chis[c] = chi;
	}		// end loop for all the channels


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calucaltes the skip parameter of the main stem (the longest channel
// The maximum length of the dataset will be in the main stem so this will determine the
// target spacing of all the tributaries
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiNetwork::calculate_skip(int target_nodes)
{
	int n_nodes = chis[0].size();

	double sk = double(n_nodes)/double(target_nodes);

	int rem = n_nodes%target_nodes;
	int N;

	if (sk >= 1.5)
	{
		N = int(sk);
		//cout << "fourthirds is: " << fourthirds << " and sk is: " << sk << endl;
	}
	else
	{
		if (sk < 1)
		{
			N = 0;
		}
		else
		{
			int sk2=n_nodes/rem;
			N = -(sk2-1);
		}
	}

	return N;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calucaltes the skip parameter based on a vector of chi values
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiNetwork::calculate_skip(int target_nodes, vector<double>& sorted_chis)
{
	int n_nodes = sorted_chis.size();

	double sk = double(n_nodes)/double(target_nodes);

	int rem = n_nodes%target_nodes;
	int N;

	if (sk >= 1.5)
	{
		N = int(sk);
		//cout << "fourthirds is: " << fourthirds << " and sk is: " << sk << endl;
	}
	else
	{
		if (sk < 1)
		{
			N = 0;
		}
		else
		{
			int sk2=n_nodes/rem;
			N = -(sk2-1);
		}
	}

	return N;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calucaltes the skip parameter of the main stem (the longest channel
// The maximum length of the dataset will be in the main stem so this will determine the
// target spacing of all the tributaries
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiNetwork::calculate_skip(int target_nodes, int channel_number)
{

	if (channel_number > int(elevations.size()) )
	{
		channel_number = 0;
	}

	int n_nodes = chis[channel_number].size();

	double sk = double(n_nodes)/double(target_nodes);

	int rem = n_nodes%target_nodes;
	int N;

	if (sk >= 1.5)
	{
		N = int(sk);
		//cout << "fourthirds is: " << fourthirds << " and sk is: " << sk << endl;
	}
	else
	{
		if (sk < 1)
		{
			N = 0;
		}
		else
		{
			int sk2=n_nodes/rem;
			N = -(sk2-1);
		}
	}

	return N;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calucaltes the chi spacing of the main stem channel (the longest channel
// The maximum length of the dataset will be in the main stem so this will determine the
// target spacing of all the tributaries
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::calculate_optimal_chi_spacing(int target_nodes)
{
	double dchi;
	vector<double>::iterator viter_begin = chis[0].begin();
	vector<double>::iterator viter_end= chis[0].end();
	viter_end--;			// this is necessary since the .end() member function
							// gets the value of the vector one past the end
	//cout << "LSDChiNetwork line 245 begin: " << *viter_begin << " and end: " << *viter_end << endl;
	double chi_length = *viter_begin-*viter_end;
	//cout << "LSDChiNetwork line 247 chi length is: " << chi_length << endl;
	dchi = chi_length/(double(target_nodes));

	return dchi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function gets the most likely channel segments for a particular channel
// channel is the index into the channel
// miniumum segment length is how many nodes the mimimum segment will have
// sigma is the standard deviation of error on elevation data
// this function replaces the b, m, r2 and DW values of each segment into vectors
// it also returns the fitted elevation and the index into the original channel (since this is done
//		with thinned data)
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::find_most_likeley_segments(int channel,
						 int minimum_segment_length,
						 double sigma, int N, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc )
{
	vector<int> empty_vec;

	vector<double> reverse_Chi = chis[channel];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[channel];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());
	vector<int> node_ref;

	LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

	//int n_nodes = reverse_Chi.size();
	channel_MLE_finder.thin_data_skip(N, node_ref);

	// now create a single sigma value vector
	vector<double> sigma_values;
	sigma_values.push_back(sigma);

	// compute the best fit AIC
	channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

	channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
					     	r2_vec, DW_vec, fitted_elev,these_segment_lengths,
					       	this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

	// get the thinned chi locations, elevation of the data at these locations
	// the fitted data is also returned in the fitted_elev vector
	thinned_chi = channel_MLE_finder.get_x_data();
	thinned_elev = channel_MLE_finder.get_y_data();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=







//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function gets the most likely channel segments for a particular channel
// channel is the index into the channel
// miniumum segment length is how many nodes the mimimum segment will have
// sigma is the standard deviation of error on elevation data
// this function replaces the b, m, r2 and DW values of each segment into vectors
// it also returns the fitted elevation and the index into the original channel (since this is done
//		with thinned data)
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::find_most_likeley_segments_dchi(int channel,
						 int minimum_segment_length,
						 double sigma, double dchi, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc )
{
	vector<int> empty_vec;

	vector<double> reverse_Chi = chis[channel];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[channel];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());

	LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

	int n_nodes = reverse_Chi.size();
	channel_MLE_finder.thin_data_target_dx_preserve_data(dchi, node_reference);
	n_nodes = node_reference.size();

	// now create a single sigma value vector
	vector<double> sigma_values;
	sigma_values.push_back(sigma);

	// compute the best fit AIC
	channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

	channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
					     	r2_vec, DW_vec, fitted_elev,these_segment_lengths,
					       	this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

	// get the thinned chi locations, elevation of the data at these locations
	// the fitted data is also returned in the fitted_elev vector
	thinned_chi = channel_MLE_finder.get_x_data();
	thinned_elev = channel_MLE_finder.get_y_data();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the most likely segments but uses the monte carlo data thinning method
// the mean_dchi is the mean value of dchi, and the variation is how miuch the
// chi chan vary, such that minimum_dchi = dchi-variation_dchi.
// The expectation is that this will be used repeatedly on channels to generate statistics of the
// best fit segments by individual nodes in the channel network.
//
// SMM 01/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::find_most_likeley_segments_monte_carlo(int channel, int minimum_segment_length,
						 double sigma, int mean_skip, int skip_range,  vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc )
{
	// first create a segment finder object
        //cout << "making MLEfinder object, " << endl;
	vector<int> empty_vec;

	vector<double> reverse_Chi = chis[channel];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[channel];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());


	LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

	int n_nodes = reverse_Chi.size();

	// now thin the data, preserving the data (not interpolating)
	channel_MLE_finder.thin_data_monte_carlo_skip(mean_skip, skip_range, node_reference);
	n_nodes = node_reference.size();

	// now create a single sigma value vector
	vector<double> sigma_values;
	sigma_values.push_back(sigma);

	// compute the best fit AIC
	channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);


	channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
					     	r2_vec, DW_vec, fitted_elev,these_segment_lengths,
					       	this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

	// get the thinned chi locations, elevation of the data at these locations
	// the fitted data is also returned in the fitted_elev vector
	thinned_chi = channel_MLE_finder.get_x_data();
	thinned_elev = channel_MLE_finder.get_y_data();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the most likely segments but uses the monte carlo data thinning method
// the mean_dchi is the mean value of dchi, and the variation is how miuch the
// chi chan vary, such that minimum_dchi = dchi-variation_dchi.
// The expectation is that this will be used repeatedly on channels to generate statistics of the
// best fit segments by individual nodes in the channel network.
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::find_most_likeley_segments_monte_carlo_dchi(int channel, int minimum_segment_length,
						 double sigma, double mean_dchi, double variation_dchi, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc )
{
	// first create a segment finder object
        //cout << "making MLEfinder object, " << endl;
	vector<int> empty_vec;

	vector<double> reverse_Chi = chis[channel];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[channel];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());


	LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

	int n_nodes = reverse_Chi.size();

	// now thin the data, preserving the data (not interpolating)
	channel_MLE_finder.thin_data_monte_carlo_dchi(mean_dchi, variation_dchi, node_reference);
	n_nodes = node_reference.size();

	// now create a single sigma value vector
	vector<double> sigma_values;
	sigma_values.push_back(sigma);

	// compute the best fit AIC
	channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);


	channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
					     	r2_vec, DW_vec, fitted_elev,these_segment_lengths,
					       	this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);


	// get the thinned chi locations, elevation of the data at these locations
	// the fitted data is also returned in the fitted_elev vector
	thinned_chi = channel_MLE_finder.get_x_data();
	thinned_elev = channel_MLE_finder.get_y_data();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Monte carlo segment fitter
// This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
// standard error information
//
// the fraction_dchi_for_variation is the fration of the optimal dchi that dchi can vary over. So for example
// if this = 0.4 then the variation of dchi will be 0.4*mean_dchi and the minimum dchi will be
// min_dchi = (1-0.4)*mean_dchi
//
// Note: this is _extremely_ computationally and data intensive.
//
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit(double A_0, double m_over_n, int n_iterations,
				int mean_skip, int skip_range,
				int minimum_segment_length, double sigma)
{

  	int n_channels = chis.size();
	calculate_chi(A_0, m_over_n);
	m_over_n_for_fitted_data = m_over_n;		// store this m_over_n value
	A_0_for_fitted_data = A_0;

	// vecvecvecs for storing information. The top level
	// is the channel. The second level is the node
	// the third level is the individual data elements
  	vector< vector< vector<double> > > b_vecvecvec(n_channels);
  	vector< vector< vector<double> > > m_vecvecvec(n_channels);
  	vector< vector< vector<double> > > DW_vecvecvec(n_channels);
  	vector< vector< vector<double> > > r2_vecvecvec(n_channels);
  	vector< vector< vector<double> > > fitted_elev_vecvecvec(n_channels);
  	//vector< vector< vector<int> > > these_segment_lengths_vecvec(n_channels);

  	// iterators to navigate these data elements
  	vector< vector< vector<double> > >::iterator first_level_doub_iter;
  	vector< vector< vector<double> > >::iterator second_level_doub_iter;

	// now expand all of these vecvecvecs to be the correct size
	for (int cn = 0; cn<n_channels; cn++)
	{
		vector< vector<double> > temp_vecvec_double;
		//vector< vector<int> > temp_vecvec_int;
		vector<double> empty_double_vec;
		//vector<int> empty_int_vec;

		// expand the vecvec elements
		int nodes_in_channel = int(chis[cn].size());
		for (int n = 0; n<nodes_in_channel; n++)
		{
			temp_vecvec_double.push_back(empty_double_vec);
			//temp_vecvec_int.push_back(empty_int_vec);
		}

		// now add this to the vecvecvecs
		b_vecvecvec[cn] = temp_vecvec_double;
		m_vecvecvec[cn]= temp_vecvec_double;
		DW_vecvecvec[cn] = temp_vecvec_double;
		r2_vecvecvec[cn] =temp_vecvec_double;
		fitted_elev_vecvecvec[cn] = temp_vecvec_double;
		//these_segment_lengths_vecvec[cn] = temp_vecvec_int;
	}

	// now the vectors that will be replaced by the fitting algorithm
  	// theyare from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;
	int n_segments;

	// loop through the iterations
	for (int it = 1; it <= n_iterations; it++)
	{
		if (it%10 == 0)
		{
			cout << "LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit, iteration: " << it << endl;
		}

      	// now loop through channels
      	for (int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << endl << "LINE 501 LSDChiNetwork channel: " << chan << endl;

			// get the most liekely segments
	  	  	find_most_likeley_segments_monte_carlo(chan,minimum_segment_length, sigma,
	  	  								mean_skip, skip_range,
					    				b_vec, m_vec, r2_vec, DW_vec, chi_thinned, elev_thinned,
					   				 	elev_fitted, node_ref_thinned, these_segment_lengths,
					    				this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

			// print the segment properties for bug checking
			//cout << " channel number: " << chan << " n_segments: " << these_segment_lengths.size() << endl;
			//for (int i = 0; i< int(b_vec.size()); i++)
			//{
			//	cout << "segment: " << i << " m: " << m_vec[i] << " b: " << b_vec[i] << endl;
			//}

			// now assign the m, b, r2 and DW values for segments to all the nodes in the thinned data
			vector<double> m_per_node;
			vector<double> b_per_node;
			vector<double> r2_per_node;
			vector<double> DW_per_node;
			n_segments = int(b_vec.size());
			for(int seg = 0; seg<n_segments; seg++)
			{
				//cout << "segment length: " << these_segment_lengths[seg] << endl;
				for(int n = 0; n < these_segment_lengths[seg]; n++)
				{
					m_per_node.push_back(m_vec[seg]);
					b_per_node.push_back(b_vec[seg]);
					r2_per_node.push_back(r2_vec[seg]);
					DW_per_node.push_back(DW_vec[seg]);
				}
			}

			//cout << "size thinned: " << chi_thinned.size() << " and m_node: " << m_per_node.size() << endl;

			// now assign the values of these variable to the vecvecvecs
			int n_nodes_in_chan = int(chi_thinned.size());
			for (int n = 0; n< n_nodes_in_chan; n++)
			{
				int this_node = node_ref_thinned[n];
  				b_vecvecvec[chan][this_node].push_back(b_per_node[n]);
  				m_vecvecvec[chan][this_node].push_back(m_per_node[n]);
  				DW_vecvecvec[chan][this_node].push_back(DW_per_node[n]);
  				r2_vecvecvec[chan][this_node].push_back(r2_per_node[n]);
  				fitted_elev_vecvecvec[chan][this_node].push_back(elev_fitted[n]);
			}

			// NOTE: one could include a cumualtive AIC calculator here to compare
			// multiple instances of AICc for a further test of the best fit m over n

		}			// end channel loop
	}				// end iteration loop

	// reset the data holding the fitted network properties
	vector< vector<double> > empty_vecvec;
	vector< vector<int> > empty_int_vecvec;
	chi_m_means = empty_vecvec;
	chi_m_standard_deviations = empty_vecvec;
	chi_m_standard_errors = empty_vecvec;
	chi_b_means = empty_vecvec;
	chi_b_standard_deviations = empty_vecvec;
	chi_b_standard_errors = empty_vecvec;
	all_fitted_elev_means = empty_vecvec;
	all_fitted_elev_standard_deviations = empty_vecvec;
	all_fitted_elev_standard_errors = empty_vecvec;
	chi_DW_means = empty_vecvec;
	chi_DW_standard_deviations = empty_vecvec;
	chi_DW_standard_errors = empty_vecvec;
	n_data_points_used_in_stats = empty_int_vecvec;

	// vectors that accept the data from the vecvecvec and are passed to
	// the get_common_statistics function
	vector<double> b_datavec;
	vector<double> m_datavec;
	vector<double> DW_datavec;
	vector<double> elev_datavec;
	vector<double> common_stats;

	// now go through the channel nodes and see how many data elements there are.
	for (int chan = 0; chan<n_channels; chan++)
	{
		// initialize vectors for storing the statistics of the
		// monte carlo fitted network properties
		int n_nodes_in_chan = int(chis[chan].size());

		vector<double> m_means(n_nodes_in_chan);
		vector<double> m_standard_deviations(n_nodes_in_chan);
		vector<double> m_standard_error(n_nodes_in_chan);
		vector<double> b_means(n_nodes_in_chan);
		vector<double> b_standard_deviations(n_nodes_in_chan);
		vector<double> b_standard_error(n_nodes_in_chan);
		vector<double> DW_means(n_nodes_in_chan);
		vector<double> DW_standard_deviations(n_nodes_in_chan);
		vector<double> DW_standard_error(n_nodes_in_chan);
		vector<double> fitted_elev_means(n_nodes_in_chan);
		vector<double> fitted_elev_standard_deviations(n_nodes_in_chan);
		vector<double> fitted_elev_standard_error(n_nodes_in_chan);
		vector<int> n_data_points_in_this_channel_node(n_nodes_in_chan);

		// now loop through each node, calculating how many data points there are in each
		for (int n = 0; n< n_nodes_in_chan; n++)
		{
			// get the number of data points in this channel.
			int n_data_points_itc = int(b_vecvecvec[chan][n].size());
			n_data_points_in_this_channel_node[n] = n_data_points_itc;
			b_datavec = b_vecvecvec[chan][n];
			m_datavec = m_vecvecvec[chan][n];
			DW_datavec = DW_vecvecvec[chan][n];
			elev_datavec = fitted_elev_vecvecvec[chan][n];

			// calcualte statistics, but only if there is data
			if (n_data_points_itc > 0)
			{
				common_stats = get_common_statistics(b_datavec);
				b_means[n] = common_stats[0];
				b_standard_deviations[n] = common_stats[2];
				b_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(m_datavec);
				m_means[n] = common_stats[0];
				m_standard_deviations[n] = common_stats[2];
				m_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(DW_datavec);
				DW_means[n] = common_stats[0];
				DW_standard_deviations[n] = common_stats[2];
				DW_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(elev_datavec);
				fitted_elev_means[n] = common_stats[0];
				fitted_elev_standard_deviations[n] = common_stats[2];
				fitted_elev_standard_error[n] = common_stats[3];
			}
			else	// if there is no data, use the previous node. This works because there is always data in the 1st node
			{
				b_means[n] = b_means[n-1];
				b_standard_deviations[n] = b_standard_deviations[n-1];
				b_standard_error[n] = b_standard_error[n-1];

				m_means[n] = m_means[n-1];
				m_standard_deviations[n] = m_standard_deviations[n-1];
				m_standard_error[n] = m_standard_error[n-1];

				DW_means[n] = DW_means[n-1];
				DW_standard_deviations[n] = DW_standard_deviations[n-1];
				DW_standard_error[n] = DW_standard_error[n-1];

				fitted_elev_means[n] = fitted_elev_means[n-1];
				fitted_elev_standard_deviations[n] = fitted_elev_standard_deviations[n-1];
				fitted_elev_standard_error[n] = fitted_elev_standard_error[n-1];
			}

		}	// end looping through channel nodes

		// all of the data vectors get reversed by the data fitting algorithm so they need to
		// be reinverted before they get inserted into the vecvecs
		reverse(m_means.begin(), m_means.end());
		reverse(m_standard_deviations.begin(), m_standard_deviations.end());
		reverse(m_standard_error.begin(), m_standard_error.end());

		reverse(b_means.begin(), b_means.end());
		reverse(b_standard_deviations.begin(), b_standard_deviations.end());
		reverse(b_standard_error.begin(), b_standard_error.end());

		reverse(DW_means.begin(), DW_means.end());
		reverse(DW_standard_deviations.begin(), DW_standard_deviations.end());
		reverse(DW_standard_error.begin(), DW_standard_error.end());

		reverse(fitted_elev_means.begin(), fitted_elev_means.end());
		reverse(fitted_elev_standard_deviations.begin(), fitted_elev_standard_deviations.end());
		reverse(fitted_elev_standard_error.begin(), fitted_elev_standard_error.end());

		reverse(n_data_points_in_this_channel_node.begin(),n_data_points_in_this_channel_node.end());

		// now store all the data
		chi_m_means.push_back(m_means);
		chi_m_standard_deviations.push_back(m_standard_deviations);
		chi_m_standard_errors.push_back(m_standard_error);
		chi_b_means.push_back(b_means);
		chi_b_standard_deviations.push_back(b_standard_deviations);
		chi_b_standard_errors.push_back(b_standard_error);
		chi_DW_means.push_back(DW_means);
		chi_DW_standard_deviations.push_back(DW_standard_deviations);
		chi_DW_standard_errors.push_back(DW_standard_error);
		all_fitted_elev_means.push_back(fitted_elev_means);
		all_fitted_elev_standard_deviations.push_back(fitted_elev_standard_deviations);
		all_fitted_elev_standard_errors.push_back(fitted_elev_standard_error);
		n_data_points_used_in_stats.push_back(n_data_points_in_this_channel_node);

	}				// end channel loop for processing data
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Monte carlo segment fitter
// This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
// standard error information
//
// the break nodes vector tells the algorithm where the breaks in the channel occur
// this function is called repeatedly until the target skip equals the all of the this_skip values
//
// This function continues to split the channel into segments until the target skip is achieved
//
// SMM 01/04/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::monte_carlo_split_channel(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes,
				int minimum_segment_length, double sigma, int chan, vector<int>& break_nodes)
{
	int mean_skip;
	int skip_range;
	vector<int> node_reference;

  	int n_channels = chis.size();

  	if(chan > n_channels-1)
  	{
		cout << "You have selected a channel that doesn't exist. Switching to mainstem" << endl;
		chan = 0;
	}

	calculate_chi(A_0, m_over_n);
	m_over_n_for_fitted_data = m_over_n;		// store this m_over_n value
	A_0_for_fitted_data = A_0;

	// we need to find the maximum length of the segments for skip analysis
	// if skip > 0, then the maximum number of nodes is target_nodes*skip
	// if skip == 0 then the maximum number of nodes is target_nodes
	// if skip < 0 then the maximum number of nodes is target_nodes*( (-skip+2)/(-skip+1) )
	int max_nodes_in_section;
	if (target_skip > 0)
	{
		max_nodes_in_section = target_nodes*(target_skip+1);
	}
	else if (target_skip == 0)
	{
		max_nodes_in_section = target_nodes;
	}
	else
	{
		max_nodes_in_section = (target_nodes*(-target_skip+2))/(-target_skip+1);
	}

	// get the data from the channel
	vector<double> reverse_Chi = chis[chan];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[chan];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());

	// these data members keep track of the  breaks
	int n_nodes = int(reverse_Chi.size());
	list<int> breaks;
	list<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	breaks.push_back(n_nodes-1);
	int start_of_last_break = 0;
	int length_of_break_segment;

	//cout << "LSDCN, LINE 1730, max_nodes_in_section: " << max_nodes_in_section
	//	 << " and n nodes: " << n_nodes <<  " and break: " << *(breaks.begin()) << endl;

	// now we enter a recursive splitting loop that
	// loops through the breaks, and if there is a break it
	// returns to check if that break needs to be broken
	br_iter = breaks.begin();
	while(br_iter != breaks.end())
	{
		// get the length of this break
		length_of_break_segment = (*br_iter)-start_of_last_break+1;

		//cout << "LSDCN Line 1772, break " << (*br_iter) << " and start of last break: " << start_of_last_break << endl;
		//cout << " and length of segment: " << length_of_break_segment << " and max length: " << max_nodes_in_section << endl;

		// if this break is longer than the maximum length of the
		// break segment, start the monte carlo algorithm to break it
		if (length_of_break_segment > max_nodes_in_section)
		{
			// get the skips
			double sk = double(n_nodes)/double(target_nodes);

			int rem = n_nodes%target_nodes;
			int N;

			if (sk >= 1.5)
			{
				N = int(sk);
				//cout << "fourthirds is: " << fourthirds << " and sk is: " << sk << endl;
			}
			else
			{
				if (sk < 1)
				{
					N = 0;
				}
				else
				{
					int sk2=length_of_break_segment/rem;
					N = -(sk2-1);
				}
			}
			mean_skip = N;
			skip_range = N*2;
			if (skip_range ==0)
			{
				skip_range = 2;
			}
			if (skip_range < 0)
			{
				skip_range = -skip_range;
			}

			// now prepare the vectors for the analysis
			vec_iter_start = reverse_Chi.begin()+start_of_last_break;
			vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_chi;
			br_chi.assign(vec_iter_start,vec_iter_end);

			vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
			vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_elev;
			br_elev.assign(vec_iter_start,vec_iter_end);
			int n_br = br_elev.size();

			//cout << endl;
			//cout << "main data: " << endl;
			//for (int i = 0; i<n_nodes; i++)
			//{
			//	cout << i << " " << reverse_Chi[i] << " " << reverse_Elevation[i] << endl;
			//}

			//cout << endl;
			//cout << "break data " << endl;
			//for (int i = 0; i<n_br; i++)
			//{
			//	cout << i << " " << br_chi[i] << " " << br_elev[i] << endl;
			//}

			// initiate the vecvecs for finding the breaks
			vector< vector<double> > b_vecvec(n_br);
			vector< vector<double> > m_vecvec(n_br);
			vector< vector<double> > seg_number_vecvec(n_br);

			// initiate the vectors for holding the means
			vector<double> b_means(n_br);
			vector<double> m_means(n_br);
			vector<double> seg_number_means(n_br);

			// now run the monte carlo algotithm thorugh N iterations
			for (int iteration = 0; iteration < n_iterations; iteration++)
			{
				//cout << "LSDCN Line 1841 iteration is: " << iteration << endl;

				// now the vectors that will be replaced by the fitting algorithm
				// they are from the individual channels, which are replaced each time a new channel is analyzed
				vector<double> m_vec;
				vector<double> b_vec;
				vector<double> r2_vec;
				vector<double> DW_vec;
				int n_data_nodes;
				int this_n_segments;
				double this_MLE, this_AIC, this_AICc;
				vector<int> these_segment_lengths;
				vector<double> fitted_elev;
				vector<int> n_segements_each_iteration;
				int n_segments;

				// now assign the m, b, r2 and DW values for segments to all the nodes in the thinned data
				vector<double> m_per_node;
				vector<double> b_per_node;
				vector<double> seg_number_per_node;

				// Run the monte carlo algorithm on the break segment

				LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

				// now thin the data, preserving the data (not interpolating)
				channel_MLE_finder.thin_data_monte_carlo_skip(mean_skip, skip_range, node_reference);
				n_data_nodes = node_reference.size();

				// now create a single sigma value vector
				vector<double> sigma_values;
				sigma_values.push_back(sigma);

				// compute the best fit AIC
				channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

				// get the segments
				channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
										r2_vec, DW_vec, fitted_elev,these_segment_lengths,
										this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

				//cout << "Line 1887, n_nodes: " << n_nodes << endl;

				n_segments = int(b_vec.size());
				for(int seg = 0; seg<n_segments; seg++)
				{
					//cout << "segment length: " << these_segment_lengths[seg] << endl;
					for(int n = 0; n < these_segment_lengths[seg]; n++)
					{
						m_per_node.push_back(m_vec[seg]);
						b_per_node.push_back(b_vec[seg]);
						seg_number_per_node.push_back( double(seg) );
					}
				}
				n_segements_each_iteration.push_back(this_n_segments);

				// now assign the values of these variable to the vecvecvecs
				for (int n = 0; n< n_data_nodes; n++)
				{
					int this_node = node_reference[n];
					b_vecvec[this_node].push_back(b_per_node[n]);
					m_vecvec[this_node].push_back(m_per_node[n]);
					seg_number_vecvec[this_node].push_back(seg_number_per_node[n]);
				}
			}			// finished the monte carlo iteration

			// now get the averages and look for a break
			vector<double> b_datavec;
			vector<double> m_datavec;
			vector<double> seg_number_datavec;
			vector<double> common_stats;


			for (int br_node = 0; br_node < n_br; br_node++)
			{
				b_datavec = b_vecvec[br_node];
				m_datavec = m_vecvec[br_node];
				seg_number_datavec = seg_number_vecvec[br_node];

				int n_data_points_itc = b_datavec.size();

				// calcualte statistics, but only if there is data
				if (n_data_points_itc > 0)
				{
					common_stats = get_common_statistics(b_datavec);
					b_means[br_node] = common_stats[0];

					common_stats = get_common_statistics(m_datavec);
					m_means[br_node] = common_stats[0];

					common_stats = get_common_statistics(seg_number_datavec);
					seg_number_means[br_node] = common_stats[0];
				}
				else	// if there is no data, use the previous node. This works because there is always data in the 1st node
				{
					b_means[br_node] = b_means[br_node-1];

					m_means[br_node] = m_means[br_node-1];

					seg_number_means[br_node] = seg_number_means[br_node-1];
				}
			}

			// now show the data
			//for (int br_node = 0; br_node < n_br; br_node++)
			//{
			//	cout << "i: " << br_node << " b: " << b_means[br_node] << " m: " << m_means[br_node]
			//		 << " seg_num: " << seg_number_means[br_node] << endl;
			//}

			// now we want to split the data. We do this where the segment number is intermediate
			// between two integers
			double this_seg = 0;
			vector<int> this_segment_breaks;
			for (int br_node = 0; br_node < n_br; br_node++)
			{
				if( seg_number_means[br_node] > this_seg+0.5)
				{
					//cout << "breaking at node: " << br_node-1+start_of_last_break << endl;
					this_seg = this_seg+1;
					this_segment_breaks.push_back(br_node-1+start_of_last_break);
				}
			}

			// if there is no break, break the segment in the middle
			if(this_seg == 0)
			{
				int mid_break = n_br/2;
				this_segment_breaks.push_back(mid_break+start_of_last_break);
			}


			// now insert these breaks in the break list
			int n_new_breaks = this_segment_breaks.size();
			for (int i = 0; i<n_new_breaks; i++)
			{
				breaks.insert(br_iter,this_segment_breaks[i]);
			}

			// now the break algorithm will have to revisit these breaks, so the iterator needs to go back
			//cout << "iterator now: " << *br_iter;
			for (int i = 0; i<n_new_breaks; i++)
			{
				br_iter--;
			}
			//cout << " and after: " << *br_iter << endl;

		}		// finished attempting to break the segment

		else	// if the section is within the allowed node limit, move on to the next segment
		{
			start_of_last_break = (*br_iter)+1;
			br_iter++;
		}

	}

	// put the breaks into the break_nodes vector
	br_iter = breaks.begin();
	vector<int> br_nds;
	while(br_iter != breaks.end())
	{
		//cout << "break at: " << (*br_iter) << endl;
		br_nds.push_back( (*br_iter) );
		br_iter++;
	}

	break_nodes = br_nds;
	// now that you have the splits, it is time to accumulate the data

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Monte carlo segment fitter
// This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
// standard error information
//
// the break nodes vector tells the algorithm where the breaks in the channel occur
// this function is called repeatedly until the target skip equals the all of the this_skip values
//
// This function continues to split the channel into segments until the target skip is achieved
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::monte_carlo_split_channel_colinear(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes,
				int minimum_segment_length, double sigma,
				vector<double> reverse_Chi, vector<double> reverse_Elevation, vector<int>& break_nodes)
{
	//cout << "Line 2034, starting to break the channel" << endl;

	int mean_skip;
	int skip_range;
	vector<int> node_reference;

	m_over_n_for_fitted_data = m_over_n;		// store this m_over_n value
	A_0_for_fitted_data = A_0;

	// we need to find the maximum length of the segments for skip analysis
	// if skip > 0, then the maximum number of nodes is target_nodes*skip
	// if skip == 0 then the maximum number of nodes is target_nodes
	// if skip < 0 then the maximum number of nodes is target_nodes*( (-skip+2)/(-skip+1) )
	int max_nodes_in_section;
	if (target_skip > 0)
	{
		max_nodes_in_section = target_nodes*(target_skip+1);
	}
	else if (target_skip == 0)
	{
		max_nodes_in_section = target_nodes;
	}
	else
	{
		max_nodes_in_section = (target_nodes*(-target_skip+2))/(-target_skip+1);
	}

	//cout << "target skip: " << target_skip << " and max nodes: " << max_nodes_in_section << endl;


	// these data members keep track of the  breaks
	int n_nodes = int(reverse_Chi.size());
	list<int> breaks;
	list<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	breaks.push_back(n_nodes-1);
	int start_of_last_break = 0;
	int length_of_break_segment;

	// now we enter a recursive splitting loop that
	// loops through the breaks, and if there is a break it
	// returns to check if that break needs to be broken
	br_iter = breaks.begin();
	while(br_iter != breaks.end())
	{
		// get the length of this break
		length_of_break_segment = (*br_iter)-start_of_last_break+1;

		//cout << "Line 2084, length of break: " << length_of_break_segment << endl;

		// if this break is longer than the maximum length of the
		// break segment, start the monte carlo algorithm to break it
		if (length_of_break_segment > max_nodes_in_section)
		{
			// get the skips
			double sk = double(n_nodes)/double(target_nodes)+0.5;

			int rem = n_nodes%target_nodes;
			int N;

			if (sk >= 1.5)
			{
				N = int(sk);
				//cout << "sk is: " << sk << " and N: " << N << endl;
			}
			else
			{
				if (sk < 1)
				{
					N = 0;
				}
				else
				{
					int sk2=length_of_break_segment/rem;
					N = -(sk2-1);
				}
			}
			mean_skip = N;
			skip_range = N*2;
			if (skip_range ==0)
			{
				skip_range = 2;
			}
			if (skip_range < 0)
			{
				skip_range = -skip_range;
			}

			// now prepare the vectors for the analysis
			vec_iter_start = reverse_Chi.begin()+start_of_last_break;
			vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_chi;
			br_chi.assign(vec_iter_start,vec_iter_end);

			vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
			vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_elev;
			br_elev.assign(vec_iter_start,vec_iter_end);
			int n_br = br_elev.size();

			//cout << "LINE 2136 assigned vectors" << endl;

			// initiate the vecvecs for finding the breaks
			vector< vector<double> > b_vecvec(n_br);
			vector< vector<double> > m_vecvec(n_br);
			vector< vector<double> > seg_number_vecvec(n_br);

			// initiate the vectors for holding the means
			vector<double> b_means(n_br);
			vector<double> m_means(n_br);
			vector<double> seg_number_means(n_br);

			// now run the monte carlo algotithm thorugh N iterations
			for (int iteration = 0; iteration < n_iterations; iteration++)
			{
				//cout << "LINE 2151 iteration: " << iteration << endl;

				// now the vectors that will be replaced by the fitting algorithm
				// they are from the individual channels, which are replaced each time a new channel is analyzed
				vector<double> m_vec;
				vector<double> b_vec;
				vector<double> r2_vec;
				vector<double> DW_vec;
				int n_data_nodes;
				int this_n_segments;
				double this_MLE, this_AIC, this_AICc;
				vector<int> these_segment_lengths;
				vector<double> fitted_elev;
				vector<int> n_segements_each_iteration;
				int n_segments;

				// now assign the m, b, r2 and DW values for segments to all the nodes in the thinned data
				vector<double> m_per_node;
				vector<double> b_per_node;
				vector<double> seg_number_per_node;

				// Run the monte carlo algorithm on the break segment

				LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

				// now thin the data, preserving the data (not interpolating)
				channel_MLE_finder.thin_data_monte_carlo_skip(mean_skip, skip_range, node_reference);
				n_data_nodes = node_reference.size();

				//cout << "n data Nodes: " << n_data_nodes << endl;
				//vector<double> thinned_chi = channel_MLE_finder.get_x_data();
				//vector<double> thinned_elev = channel_MLE_finder.get_y_data();
				//for (int i = 0; i< n_data_nodes; i++)
				//{
				//	cout << "nr["<< i << "]: " << node_reference[i] << " " << thinned_chi[i] << " " << thinned_elev[i] << endl;
				//}

				// now create a single sigma value vector
				vector<double> sigma_values;
				sigma_values.push_back(sigma);

				// compute the best fit AIC
				channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

				//cout << "line 2193, got linear segments" << endl;

				// get the segments
				channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
										r2_vec, DW_vec, fitted_elev,these_segment_lengths,
										this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

				//cout << "line 2200, got data" << endl;

				n_segments = int(b_vec.size());
				for(int seg = 0; seg<n_segments; seg++)
				{
					for(int n = 0; n < these_segment_lengths[seg]; n++)
					{
						m_per_node.push_back(m_vec[seg]);
						b_per_node.push_back(b_vec[seg]);
						seg_number_per_node.push_back( double(seg) );
					}
				}
				n_segements_each_iteration.push_back(this_n_segments);




				// now assign the values of these variable to the vecvecvecs
				for (int n = 0; n< n_data_nodes; n++)
				{
					int this_node = node_reference[n];
					b_vecvec[this_node].push_back(b_per_node[n]);
					m_vecvec[this_node].push_back(m_per_node[n]);
					seg_number_vecvec[this_node].push_back(seg_number_per_node[n]);
				}
			}			// finished the monte carlo iteration

			// now get the averages and look for a break
			vector<double> b_datavec;
			vector<double> m_datavec;
			vector<double> seg_number_datavec;
			vector<double> common_stats;


			for (int br_node = 0; br_node < n_br; br_node++)
			{
				b_datavec = b_vecvec[br_node];
				m_datavec = m_vecvec[br_node];
				seg_number_datavec = seg_number_vecvec[br_node];

				int n_data_points_itc = b_datavec.size();
				//cout << "LINE 2228 br_node: " << br_node << " and n data points: " << n_data_points_itc << endl;

				// calculate statistics, but only if there is data
				if (n_data_points_itc > 0)
				{
					common_stats = get_common_statistics(b_datavec);
					b_means[br_node] = common_stats[0];

					common_stats = get_common_statistics(m_datavec);
					m_means[br_node] = common_stats[0];

					common_stats = get_common_statistics(seg_number_datavec);
					seg_number_means[br_node] = common_stats[0];
				}
				else	// if there is no data, use the previous node. This works because there is always data in the 1st node
				{
					b_means[br_node] = b_means[br_node-1];

					m_means[br_node] = m_means[br_node-1];

					seg_number_means[br_node] = seg_number_means[br_node-1];
				}
			}

			// now we want to split the data. We do this where the segment number is intermediate
			// between two integers
			double this_seg = 0;
			vector<int> this_segment_breaks;
			for (int br_node = 0; br_node < n_br; br_node++)
			{
				if( seg_number_means[br_node] > this_seg+0.5)
				{
					//cout << "breaking at node: " << br_node-1+start_of_last_break << endl;
					this_seg = this_seg+1;
					this_segment_breaks.push_back(br_node-1+start_of_last_break);
				}
			}

			// if there is no break, break the segment in the middle
			if(this_seg == 0)
			{
				int mid_break = n_br/2;
				this_segment_breaks.push_back(mid_break+start_of_last_break);
			}


			// now insert these breaks in the break list
			int n_new_breaks = this_segment_breaks.size();
			for (int i = 0; i<n_new_breaks; i++)
			{
				breaks.insert(br_iter,this_segment_breaks[i]);
			}

			// now the break algorithm will have to revisit these breaks, so the iterator needs to go back
			//cout << "iterator now: " << *br_iter;
			for (int i = 0; i<n_new_breaks; i++)
			{
				br_iter--;
			}
			//cout << " and after: " << *br_iter << endl;

		}		// finished attempting to break the segment

		else	// if the section is within the allowed node limit, move on to the next segment
		{
			start_of_last_break = (*br_iter)+1;
			br_iter++;
		}

	}

	// put the breaks into the break_nodes vector
	br_iter = breaks.begin();
	vector<int> br_nds;
	while(br_iter != breaks.end())
	{
		//cout << "break at: " << (*br_iter) << endl;
		br_nds.push_back( (*br_iter) );
		br_iter++;
	}

	break_nodes = br_nds;
	// now that you have the splits, it is time to accumulate the data

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function splits all the channels and writes the data into the break_nodes_vecvec data member
//
// SMM 01/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::split_all_channels(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes, int minimum_segment_length, double sigma)
{
	int n_channels = chis.size();
	vector<int> break_nodes;
	vector< vector<int> > this_break_vecvecvec;

	// loop through the channels, breaking into smaller bits
	for (int chan = 0; chan<n_channels; chan++)
	{
		monte_carlo_split_channel(A_0, m_over_n, n_iterations, target_skip, target_nodes,
				minimum_segment_length, sigma, chan, break_nodes);

				this_break_vecvecvec.push_back(break_nodes);
	}

	break_nodes_vecvec = this_break_vecvecvec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// this calculates the AICc of a channel after it has been broken
// it also replaces the data members n_total_segments, int& n_total_nodes, double& cumulative_MLE
// so that they can be used to get a cumulative AICc of multiple channels
//
//
// SMM 01/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::calculate_AICc_after_breaks(double A_0, double m_over_n,
				int skip, int minimum_segment_length, double sigma, int chan, vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE)

{

	calculate_chi(A_0, m_over_n);

	// get the data from the channel
	vector<double> reverse_Chi = chis[chan];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[chan];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());

	// these data members keep track of the  breaks
	//int n_nodes = int(reverse_Chi.size());
	vector<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	int start_of_last_break = 0;
	int length_of_break_segment;

	// these are vectors that hold information about each segment
	vector<double> MLE_in_this_break;
	vector<double> nodes_in_this_break;
	vector<double> segments_in_this_break;

	//int n_breaks = break_nodes.size();

	// now loop though breaks, getting the best fit.
	br_iter = break_nodes.begin();
	while(br_iter != break_nodes.end())
	{
		length_of_break_segment = (*br_iter)-start_of_last_break+1;

		// get the data of this break
		vec_iter_start = reverse_Chi.begin()+start_of_last_break;
		vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
		vector<double> br_chi;
		br_chi.assign(vec_iter_start,vec_iter_end);

		vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
		vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
		vector<double> br_elev;
		br_elev.assign(vec_iter_start,vec_iter_end);

		// now the vectors that will be replaced by the fitting algorithm
		// they are from the individual channels, which are replaced each time a new channel is analyzed
		vector<double> m_vec;
		vector<double> b_vec;
		vector<double> r2_vec;
		vector<double> DW_vec;
		int n_data_nodes;
		int this_n_segments;
		double this_MLE, this_AIC, this_AICc;
		vector<int> these_segment_lengths;
		vector<double> fitted_elev;
		vector<int> n_segements_each_iteration;
		vector<int> node_reference;

		// Run the monte carlo algorithm on the break segment
		LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

		// now thin the data, preserving the data (not interpolating)
		channel_MLE_finder.thin_data_skip(skip, node_reference);
		n_data_nodes = node_reference.size();

		// now create a single sigma value vector
		vector<double> sigma_values;
		sigma_values.push_back(sigma);

		// compute the best fit AIC
		channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

		// get the segments
		channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
								r2_vec, DW_vec, fitted_elev,these_segment_lengths,
								this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

		//cout << "this_n_segments: " <<  this_n_segments << endl;
		//for(int s = 0; s<this_n_segments; s++)
		//{
		//	cout << "segment: " << s << " m: " << m_vec[s] << " b: " << b_vec[s] << endl;
		//}

		// store the information about this segment that will be used for AICc
		MLE_in_this_break.push_back(this_MLE);
		nodes_in_this_break.push_back(n_data_nodes);
		segments_in_this_break.push_back(this_n_segments);

		// reset the starting node
		start_of_last_break = (*br_iter)+1;

		br_iter++;
	}

	//now calculate the cumulative AICc
	double thismn_AIC;
	double thismn_AICc;

	n_total_segments = 0;
	n_total_nodes = 0;
	cumulative_MLE = 1;
	double log_cum_MLE = 0;

	int n_total_breaks = MLE_in_this_break.size();

	// get the cumulative maximum likelihood estimators
	for (int br = 0; br<n_total_breaks; br++)
	{
		n_total_segments += segments_in_this_break[br];
		n_total_nodes += nodes_in_this_break[br];
		cumulative_MLE = MLE_in_this_break[br]*cumulative_MLE;

		if(MLE_in_this_break[br] <= 0)
		{
			log_cum_MLE = log_cum_MLE-1000;
		}
		else
		{
			log_cum_MLE = log(MLE_in_this_break[br])+log_cum_MLE;
		}


	}

	// these AIC and AICc values are cumulative for a given m_over_n
	//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
															   // for each segment there are 2 parameters
	thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
	thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);

	return thismn_AICc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// this calculates the AICc of a channel after it has been broken
// it also replaces the data members n_total_segments, int& n_total_nodes, double& cumulative_MLE
// so that they can be used to get a cumulative AICc of multiple channels
//
// this version uses a monte carlo approach and returns AICc values in a vector
// that has the length of the number of iterations
//
//
// SMM 01/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDChiNetwork::calculate_AICc_after_breaks_monte_carlo(double A_0, double m_over_n,
				int target_skip, int minimum_segment_length, double sigma, int chan, vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE,
				int n_iterations)

{

	calculate_chi(A_0, m_over_n);

	// get the data from the channel
	vector<double> reverse_Chi = chis[chan];
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = elevations[chan];
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());

	// these data members keep track of the  breaks
	//int n_nodes = int(reverse_Chi.size());
	vector<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	int start_of_last_break;
	int length_of_break_segment;
	int skip_range = 2*target_skip;
	if (skip_range ==0)
	{
		skip_range = 2;
	}
	if (skip_range < 0)
	{
		skip_range = -skip_range;
	}

	// these are vectors that hold information about each segment
	vector<double> MLE_in_this_break;
	vector<double> nodes_in_this_break;
	vector<double> segments_in_this_break;
	vector<double> empty_vec;
	vector<double> these_AICcs;

	cout << " looping, iteration";
	for (int iteration = 0; iteration<n_iterations; iteration++)
	{

		if (iteration %50 == 0)
		{
			cout << " " << iteration;
		}

		MLE_in_this_break = empty_vec;
		nodes_in_this_break = empty_vec;
		segments_in_this_break = empty_vec;
		start_of_last_break = 0;

		// now loop though breaks, getting the best fit.
		br_iter = break_nodes.begin();
		while(br_iter != break_nodes.end())
		{
			length_of_break_segment = (*br_iter)-start_of_last_break+1;

			// get the data of this break
			vec_iter_start = reverse_Chi.begin()+start_of_last_break;
			vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_chi;
			br_chi.assign(vec_iter_start,vec_iter_end);

			vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
			vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_elev;
			br_elev.assign(vec_iter_start,vec_iter_end);

			// now the vectors that will be replaced by the fitting algorithm
			// they are from the individual channels, which are replaced each time a new channel is analyzed
			vector<double> m_vec;
			vector<double> b_vec;
			vector<double> r2_vec;
			vector<double> DW_vec;
			int n_data_nodes;
			int this_n_segments;
			double this_MLE, this_AIC, this_AICc;
			vector<int> these_segment_lengths;
			vector<double> fitted_elev;
			vector<int> n_segements_each_iteration;
			vector<int> node_reference;

			// Run the monte carlo algorithm on the break segment
			LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

			// now thin the data, preserving the data (not interpolating)
			channel_MLE_finder.thin_data_monte_carlo_skip(target_skip, skip_range, node_reference);
			n_data_nodes = node_reference.size();

			// now create a single sigma value vector
			vector<double> sigma_values;
			sigma_values.push_back(sigma);

			// compute the best fit AIC
			channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

			// get the segments
			channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
									r2_vec, DW_vec, fitted_elev,these_segment_lengths,
									this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

			//cout << "this_n_segments: " <<  this_n_segments << endl;
			//for(int s = 0; s<this_n_segments; s++)
			//{
			//	cout << "segment: " << s << " m: " << m_vec[s] << " b: " << b_vec[s] << endl;
			//}

			// store the information about this segment that will be used for AICc
			MLE_in_this_break.push_back(this_MLE);
			nodes_in_this_break.push_back(n_data_nodes);
			segments_in_this_break.push_back(this_n_segments);

			// reset the starting node
			start_of_last_break = (*br_iter)+1;

			br_iter++;
		}

		//now calculate the cumulative AICc
		double thismn_AIC;
		double thismn_AICc;

		n_total_segments = 0;
		n_total_nodes = 0;
		cumulative_MLE = 1;
		double log_cum_MLE = 0;

		int n_total_breaks = MLE_in_this_break.size();

		// get the cumulative maximum likelihood estimators
		for (int br = 0; br<n_total_breaks; br++)
		{
			n_total_segments += segments_in_this_break[br];
			n_total_nodes += nodes_in_this_break[br];
			cumulative_MLE = MLE_in_this_break[br]*cumulative_MLE;

			if(MLE_in_this_break[br] <= 0)
			{
				log_cum_MLE = log_cum_MLE-1000;
			}
			else
			{
				log_cum_MLE = log(MLE_in_this_break[br])+log_cum_MLE;
			}

		}

		// these AIC and AICc values are cumulative for a given m_over_n
		//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
																   // for each segment there are 2 parameters
		thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
		thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);

		these_AICcs.push_back(thismn_AICc);
	}
	cout << endl;

	return these_AICcs;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// this calculates the AICc of a channel after it has been broken. Uses a combined colinear channel.
// it also replaces the data members n_total_segments, int& n_total_nodes, double& cumulative_MLE
// so that they can be used to get a cumulative AICc of multiple channels
//
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::calculate_AICc_after_breaks_colinear(double A_0, double m_over_n,
				int skip, int minimum_segment_length, double sigma,
				vector<double> reverse_Chi, vector<double> reverse_Elevation,
				vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE)
{

	// these data members keep track of the  breaks
	//int n_nodes = int(reverse_Chi.size());
	vector<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	int start_of_last_break = 0;
	int length_of_break_segment;

	// these are vectors that hold information about each segment
	vector<double> MLE_in_this_break;
	vector<double> nodes_in_this_break;
	vector<double> segments_in_this_break;

	//int n_breaks = break_nodes.size();

	// now loop though breaks, getting the best fit.
	br_iter = break_nodes.begin();
	while(br_iter != break_nodes.end())
	{
		length_of_break_segment = (*br_iter)-start_of_last_break+1;

		// get the data of this break
		vec_iter_start = reverse_Chi.begin()+start_of_last_break;
		vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
		vector<double> br_chi;
		br_chi.assign(vec_iter_start,vec_iter_end);

		vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
		vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
		vector<double> br_elev;
		br_elev.assign(vec_iter_start,vec_iter_end);

		// now the vectors that will be replaced by the fitting algorithm
		// they are from the individual channels, which are replaced each time a new channel is analyzed
		vector<double> m_vec;
		vector<double> b_vec;
		vector<double> r2_vec;
		vector<double> DW_vec;
		int n_data_nodes;
		int this_n_segments;
		double this_MLE, this_AIC, this_AICc;
		vector<int> these_segment_lengths;
		vector<double> fitted_elev;
		vector<int> n_segements_each_iteration;
		vector<int> node_reference;

		// Run the algorithm on the break segment
		LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

		// now thin the data, preserving the data (not interpolating)
		channel_MLE_finder.thin_data_skip(skip, node_reference);
		n_data_nodes = node_reference.size();

		// now create a single sigma value vector
		vector<double> sigma_values;
		sigma_values.push_back(sigma);

		// compute the best fit AIC
		channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

		// get the segments
		channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
								r2_vec, DW_vec, fitted_elev,these_segment_lengths,
								this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

		//cout << "this_n_segments: " <<  this_n_segments << endl;
		//for(int s = 0; s<this_n_segments; s++)
		//{
		//	cout << "segment: " << s << " m: " << m_vec[s] << " b: " << b_vec[s] << endl;
		//}

		// store the information about this segment that will be used for AICc
		MLE_in_this_break.push_back(this_MLE);
		nodes_in_this_break.push_back(n_data_nodes);
		segments_in_this_break.push_back(this_n_segments);

		// reset the starting node
		start_of_last_break = (*br_iter)+1;

		br_iter++;
	}

	//now calculate the cumulative AICc
	double thismn_AIC;
	double thismn_AICc;

	n_total_segments = 0;
	n_total_nodes = 0;
	cumulative_MLE = 1;
	double log_cum_MLE = 0;

	int n_total_breaks = MLE_in_this_break.size();

	// get the cumulative maximum likelihood estimators
	for (int br = 0; br<n_total_breaks; br++)
	{
		n_total_segments += segments_in_this_break[br];
		n_total_nodes += nodes_in_this_break[br];
		cumulative_MLE = MLE_in_this_break[br]*cumulative_MLE;

		if(MLE_in_this_break[br] <= 0)
		{
			log_cum_MLE = log_cum_MLE-1000;
		}
		else
		{
			log_cum_MLE = log(MLE_in_this_break[br])+log_cum_MLE;
		}
	}

	// these AIC and AICc values are cumulative for a given m_over_n
	//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
															   // for each segment there are 2 parameters
	thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
	thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);

	return thismn_AICc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// this calculates the AICc of a channel after it has been broken. Uses a combined colinear channel.
// it also replaces the data members n_total_segments, int& n_total_nodes, double& cumulative_MLE
// so that they can be used to get a cumulative AICc of multiple channels
//
// this is the montecarlo version, it returns a vecotr of AICc values from each iteration
//
//
// SMM 01/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDChiNetwork::calculate_AICc_after_breaks_colinear_monte_carlo(double A_0, double m_over_n,
				int skip, int minimum_segment_length, double sigma,
				vector<double> reverse_Chi, vector<double> reverse_Elevation,
				vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE,
				int n_iterations)
{

	// these data members keep track of the  breaks
	//int n_nodes = int(reverse_Chi.size());
	vector<int>::iterator br_iter;
	vector<double>::iterator vec_iter_start;
	vector<double>::iterator vec_iter_end;

	int start_of_last_break;
	int length_of_break_segment;

	// these are vectors that hold information about each segment
	vector<double> MLE_in_this_break;
	vector<double> nodes_in_this_break;
	vector<double> segments_in_this_break;
	vector<double> empty_vec;
	vector<double> AICc_values;


	int skip_range = 2*skip;
	if (skip_range ==0)
	{
		skip_range = 2;
	}
	if (skip_range < 0)
	{
		skip_range = -skip_range;
	}

	cout << " LSDChiNetwork Line 2615, iteration ";
	for (int iteration = 0; iteration<n_iterations; iteration++)
	{
		//cout << "iter: " << iteration << endl;
		if (iteration % 50 == 0)
		{
			cout << " " << iteration;
		}

		start_of_last_break = 0;
		MLE_in_this_break = empty_vec;
		nodes_in_this_break = empty_vec;
		segments_in_this_break = empty_vec;

		// now loop though breaks, getting the best fit.
		br_iter = break_nodes.begin();
		int br = 0;
		while(br_iter != break_nodes.end())
		{
			//cout << "break is: " << br << " and segment is: " << (*br_iter) << endl;
			br++;

			length_of_break_segment = (*br_iter)-start_of_last_break+1;
			//cout << "length of segment is: " << length_of_break_segment << endl;

			// get the data of this break
			vec_iter_start = reverse_Chi.begin()+start_of_last_break;
			vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_chi;
			br_chi.assign(vec_iter_start,vec_iter_end);

			vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
			vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
			vector<double> br_elev;
			br_elev.assign(vec_iter_start,vec_iter_end);

			// now the vectors that will be replaced by the fitting algorithm
			// they are from the individual channels, which are replaced each time a new channel is analyzed
			vector<double> m_vec;
			vector<double> b_vec;
			vector<double> r2_vec;
			vector<double> DW_vec;
			int n_data_nodes;
			int this_n_segments;
			double this_MLE, this_AIC, this_AICc;
			vector<int> these_segment_lengths;
			vector<double> fitted_elev;
			vector<int> n_segements_each_iteration;
			vector<int> node_reference;

			// Run the algorithm on the break segment
			LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

			// now thin the data, preserving the data (not interpolating)
			channel_MLE_finder.thin_data_monte_carlo_skip(skip, skip_range, node_reference);
			n_data_nodes = node_reference.size();

			// now create a single sigma value vector
			vector<double> sigma_values;
			sigma_values.push_back(sigma);

			// compute the best fit AIC
			channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

			// get the segments
			channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
									r2_vec, DW_vec, fitted_elev,these_segment_lengths,
									this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

			// store the information about this segment that will be used for AICc
			MLE_in_this_break.push_back(this_MLE);
			nodes_in_this_break.push_back(n_data_nodes);
			segments_in_this_break.push_back(this_n_segments);

			// reset the starting node
			start_of_last_break = (*br_iter)+1;

			br_iter++;
		}
		//cout << "did breaks" << endl;


		//now calculate the cumulative AICc
		double thismn_AIC;
		double thismn_AICc;

		n_total_segments = 0;
		n_total_nodes = 0;
		cumulative_MLE = 1;

		double log_cum_MLE = 0;

		int n_total_breaks = MLE_in_this_break.size();

		// get the cumulative maximum likelihood estimators
		for (int br = 0; br<n_total_breaks; br++)
		{
			n_total_segments += segments_in_this_break[br];
			n_total_nodes += nodes_in_this_break[br];
			cumulative_MLE = MLE_in_this_break[br]*cumulative_MLE;

			if(MLE_in_this_break[br] <= 0)
			{
				log_cum_MLE = log_cum_MLE-1000;
			}
			else
			{
				log_cum_MLE = log(MLE_in_this_break[br])+log_cum_MLE;
			}

			//cout << "break: " << br << " segs: " << segments_in_this_break[br]
			//     << " nodes: " << nodes_in_this_break[br] << " MLE: " << MLE_in_this_break[br] << " and log: " << log(MLE_in_this_break[br])
			//     << " and logcum: " << log_cum_MLE << endl;

		}
		//cout << "got cumulative" << endl;

		// these AIC and AICc values are cumulative for a given m_over_n
		//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
		//														   // for each segment there are 2 parameters

		thismn_AIC = 4*n_total_segments-2*log_cum_MLE;

		//cout << "TESTING AIC, test: " << test_AIC << " old_version: " << thismn_AIC << endl;

		thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);

		AICc_values.push_back(thismn_AICc);
		//cout << "this_AICc: " << thismn_AICc << endl;
	}

	return AICc_values;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=







//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Monte carlo segment fitter
// This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
// standard error information
//
// the fraction_dchi_for_variation is the fration of the optimal dchi that dchi can vary over. So for example
// if this = 0.4 then the variation of dchi will be 0.4*mean_dchi and the minimum dchi will be
// min_dchi = (1-0.4)*mean_dchi
//
// Note: this is _extremely_ computationally and data intensive.
//
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit_dchi(double A_0, double m_over_n, int n_iterations,
				double fraction_dchi_for_variation,
				int minimum_segment_length, double sigma, int target_nodes_mainstem)
{

  	int n_channels = chis.size();
	calculate_chi(A_0, m_over_n);
	m_over_n_for_fitted_data = m_over_n;		// store this m_over_n value
	A_0_for_fitted_data = A_0;

	// set up sacing and variation of chi for iterative thinning
  	double mean_dchi;			// spacing of the chi vector. This is calculated to be optimal for the
  									// main stem
	mean_dchi = calculate_optimal_chi_spacing(target_nodes_mainstem);
	double dchi_variation = mean_dchi*fraction_dchi_for_variation;

	// vecvecvecs for storing information. The top level
	// is the channel. The second level is the node
	// the third level is the individual data elements
  	vector< vector< vector<double> > > b_vecvecvec(n_channels);
  	vector< vector< vector<double> > > m_vecvecvec(n_channels);
  	vector< vector< vector<double> > > DW_vecvecvec(n_channels);
  	vector< vector< vector<double> > > r2_vecvecvec(n_channels);
  	vector< vector< vector<double> > > fitted_elev_vecvecvec(n_channels);
  	//vector< vector< vector<int> > > these_segment_lengths_vecvec(n_channels);

  	// iterators to navigate these data elements
  	vector< vector< vector<double> > >::iterator first_level_doub_iter;
  	vector< vector< vector<double> > >::iterator second_level_doub_iter;
  	//vector< vector< vector<int> > >::iterator first_level_int_iter;
  	//vector< vector< vector<int> > >::iterator second_level_int_iter;

	// now expand all of these vecvecvecs to be the correct size
	for (int cn = 0; cn<n_channels; cn++)
	{
		vector< vector<double> > temp_vecvec_double;
		//vector< vector<int> > temp_vecvec_int;
		vector<double> empty_double_vec;
		//vector<int> empty_int_vec;

		// expand the vecvec elements
		int nodes_in_channel = int(chis[cn].size());
		for (int n = 0; n<nodes_in_channel; n++)
		{
			temp_vecvec_double.push_back(empty_double_vec);
			//temp_vecvec_int.push_back(empty_int_vec);
		}

		// now add this to the vecvecvecs
		b_vecvecvec[cn] = temp_vecvec_double;
		m_vecvecvec[cn]= temp_vecvec_double;
		DW_vecvecvec[cn] = temp_vecvec_double;
		r2_vecvecvec[cn] =temp_vecvec_double;
		fitted_elev_vecvecvec[cn] = temp_vecvec_double;
		//these_segment_lengths_vecvec[cn] = temp_vecvec_int;
	}

	// now the vectors that will be replaced by the fitting algorithm
  	// theyare from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;
	int n_segments;

	// loop through the iterations
	for (int it = 1; it <= n_iterations; it++)
	{
		if (it%10 == 0)
		{
			cout << "LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit, iteration: " << it << endl;
		}

      	// now loop through channels
      	for (int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << endl << "LINE 501 LSDChiNetwork channel: " << chan << endl;

			// get the most liekely segments
	  	  	find_most_likeley_segments_monte_carlo_dchi(chan,minimum_segment_length, sigma,
	  	  								mean_dchi, dchi_variation,
					    				b_vec, m_vec, r2_vec, DW_vec, chi_thinned, elev_thinned,
					   				 	elev_fitted, node_ref_thinned, these_segment_lengths,
					    				this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

			// print the segment properties for bug checking
			//cout << " channel number: " << chan << " n_segments: " << these_segment_lengths.size() << endl;
			//for (int i = 0; i< int(b_vec.size()); i++)
			//{
			//	cout << "segment: " << i << " m: " << m_vec[i] << " b: " << b_vec[i] << endl;
			//}

			// now assign the m, b, r2 and DW values for segments to all the nodes in the thinned data
			vector<double> m_per_node;
			vector<double> b_per_node;
			vector<double> r2_per_node;
			vector<double> DW_per_node;
			n_segments = int(b_vec.size());
			for(int seg = 0; seg<n_segments; seg++)
			{
				//cout << "segment length: " << these_segment_lengths[seg] << endl;
				for(int n = 0; n < these_segment_lengths[seg]; n++)
				{
					m_per_node.push_back(m_vec[seg]);
					b_per_node.push_back(b_vec[seg]);
					r2_per_node.push_back(r2_vec[seg]);
					DW_per_node.push_back(DW_vec[seg]);
				}
			}

			//cout << "size thinned: " << chi_thinned.size() << " and m_node: " << m_per_node.size() << endl;

			// now assign the values of these variable to the vecvecvecs
			int n_nodes_in_chan = int(chi_thinned.size());
			for (int n = 0; n< n_nodes_in_chan; n++)
			{
				int this_node = node_ref_thinned[n];
  				b_vecvecvec[chan][this_node].push_back(b_per_node[n]);
  				m_vecvecvec[chan][this_node].push_back(m_per_node[n]);
  				DW_vecvecvec[chan][this_node].push_back(DW_per_node[n]);
  				r2_vecvecvec[chan][this_node].push_back(r2_per_node[n]);
  				fitted_elev_vecvecvec[chan][this_node].push_back(elev_fitted[n]);
			}

			// NOTE: one could include a cumualtive AIC calculator here to compare
			// multiple instances of AICc for a further test of the best fit m over n

		}			// end channel loop
	}				// end iteration loop

	// reset the data holding the fitted network properties
	vector< vector<double> > empty_vecvec;
	vector< vector<int> > empty_int_vecvec;
	chi_m_means = empty_vecvec;
	chi_m_standard_deviations = empty_vecvec;
	chi_m_standard_errors = empty_vecvec;
	chi_b_means = empty_vecvec;
	chi_b_standard_deviations = empty_vecvec;
	chi_b_standard_errors = empty_vecvec;
	all_fitted_elev_means = empty_vecvec;
	all_fitted_elev_standard_deviations = empty_vecvec;
	all_fitted_elev_standard_errors = empty_vecvec;
	chi_DW_means = empty_vecvec;
	chi_DW_standard_deviations = empty_vecvec;
	chi_DW_standard_errors = empty_vecvec;
	n_data_points_used_in_stats = empty_int_vecvec;

	// vectors that accept the data from the vecvecvec and are passed to
	// the get_common_statistics function
	vector<double> b_datavec;
	vector<double> m_datavec;
	vector<double> DW_datavec;
	vector<double> elev_datavec;
	vector<double> common_stats;

	// now go through the channel nodes and see how many data elements there are.
	for (int chan = 0; chan<n_channels; chan++)
	{
		// initialize vectors for storing the statistics of the
		// monte carlo fitted network properties
		int n_nodes_in_chan = int(chis[chan].size());

		vector<double> m_means(n_nodes_in_chan);
		vector<double> m_standard_deviations(n_nodes_in_chan);
		vector<double> m_standard_error(n_nodes_in_chan);
		vector<double> b_means(n_nodes_in_chan);
		vector<double> b_standard_deviations(n_nodes_in_chan);
		vector<double> b_standard_error(n_nodes_in_chan);
		vector<double> DW_means(n_nodes_in_chan);
		vector<double> DW_standard_deviations(n_nodes_in_chan);
		vector<double> DW_standard_error(n_nodes_in_chan);
		vector<double> fitted_elev_means(n_nodes_in_chan);
		vector<double> fitted_elev_standard_deviations(n_nodes_in_chan);
		vector<double> fitted_elev_standard_error(n_nodes_in_chan);
		vector<int> n_data_points_in_this_channel_node(n_nodes_in_chan);


		// now loop through each node, calculating how many data points there are in each
		for (int n = 0; n< n_nodes_in_chan; n++)
		{
			// get the number of data points in this channel.
			int n_data_points_itc = int(b_vecvecvec[chan][n].size());
			n_data_points_in_this_channel_node[n] = n_data_points_itc;


			b_datavec = b_vecvecvec[chan][n];
			m_datavec = m_vecvecvec[chan][n];
			DW_datavec = DW_vecvecvec[chan][n];
			elev_datavec = fitted_elev_vecvecvec[chan][n];

			// calcualte statistics, but only if there is data
			if (n_data_points_itc > 0)
			{
				common_stats = get_common_statistics(b_datavec);
				b_means[n] = common_stats[0];
				b_standard_deviations[n] = common_stats[2];
				b_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(m_datavec);
				m_means[n] = common_stats[0];
				m_standard_deviations[n] = common_stats[2];
				m_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(DW_datavec);
				DW_means[n] = common_stats[0];
				DW_standard_deviations[n] = common_stats[2];
				DW_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(elev_datavec);
				fitted_elev_means[n] = common_stats[0];
				fitted_elev_standard_deviations[n] = common_stats[2];
				fitted_elev_standard_error[n] = common_stats[3];
			}
			else	// if there is no data, use the previous node. This works because there is always data in the 1st node
			{
				b_means[n] = b_means[n-1];
				b_standard_deviations[n] = b_standard_deviations[n-1];
				b_standard_error[n] = b_standard_error[n-1];

				m_means[n] = m_means[n-1];
				m_standard_deviations[n] = m_standard_deviations[n-1];
				m_standard_error[n] = m_standard_error[n-1];

				DW_means[n] = DW_means[n-1];
				DW_standard_deviations[n] = DW_standard_deviations[n-1];
				DW_standard_error[n] = DW_standard_error[n-1];

				fitted_elev_means[n] = fitted_elev_means[n-1];
				fitted_elev_standard_deviations[n] = fitted_elev_standard_deviations[n-1];
				fitted_elev_standard_error[n] = fitted_elev_standard_error[n-1];
			}

		}	// end looping through channel nodes

		// all of the data vectors get reversed by the data fitting algorithm so they need to
		// be reinverted before they get inserted into the vecvecs
		reverse(m_means.begin(), m_means.end());
		reverse(m_standard_deviations.begin(), m_standard_deviations.end());
		reverse(m_standard_error.begin(), m_standard_error.end());

		reverse(b_means.begin(), b_means.end());
		reverse(b_standard_deviations.begin(), b_standard_deviations.end());
		reverse(b_standard_error.begin(), b_standard_error.end());

		reverse(DW_means.begin(), DW_means.end());
		reverse(DW_standard_deviations.begin(), DW_standard_deviations.end());
		reverse(DW_standard_error.begin(), DW_standard_error.end());

		reverse(fitted_elev_means.begin(), fitted_elev_means.end());
		reverse(fitted_elev_standard_deviations.begin(), fitted_elev_standard_deviations.end());
		reverse(fitted_elev_standard_error.begin(), fitted_elev_standard_error.end());

		reverse(n_data_points_in_this_channel_node.begin(),n_data_points_in_this_channel_node.end());

		// now store all the data
		chi_m_means.push_back(m_means);
		chi_m_standard_deviations.push_back(m_standard_deviations);
		chi_m_standard_errors.push_back(m_standard_error);
		chi_b_means.push_back(b_means);
		chi_b_standard_deviations.push_back(b_standard_deviations);
		chi_b_standard_errors.push_back(b_standard_error);
		chi_DW_means.push_back(DW_means);
		chi_DW_standard_deviations.push_back(DW_standard_deviations);
		chi_DW_standard_errors.push_back(DW_standard_error);
		all_fitted_elev_means.push_back(fitted_elev_means);
		all_fitted_elev_standard_deviations.push_back(fitted_elev_standard_deviations);
		all_fitted_elev_standard_errors.push_back(fitted_elev_standard_error);
		n_data_points_used_in_stats.push_back(n_data_points_in_this_channel_node);

	}				// end channel loop for processing data
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// this function samples the river network using monte carlo sampling but after breaking the channels
//
//
// SMM 15/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit_after_breaks(double A_0, double m_over_n,
				int n_iterations, int skip, int minimum_segment_length, double sigma)

{
	// get the contributing channel and downstream chi
	if (break_nodes_vecvec.size() == 0)
	{
				cout << "You have not run the break algorithm on all the channels" << endl;
				exit(EXIT_FAILURE);
	}


	int n_channels = chis.size();
	calculate_chi(A_0, m_over_n);
	m_over_n_for_fitted_data = m_over_n;		// store this m_over_n value
	A_0_for_fitted_data = A_0;
	int skip_range = 2*skip;
	if (skip_range ==0)
	{
		skip_range = 2;
	}
	if (skip_range < 0)
	{
		skip_range = -skip_range;
	}

	// vecvecvecs for storing information. The top level
	// is the channel. The second level is the node
	// the third level is the individual data elements
  	vector< vector< vector<double> > > b_vecvecvec(n_channels);
  	vector< vector< vector<double> > > m_vecvecvec(n_channels);
  	vector< vector< vector<double> > > DW_vecvecvec(n_channels);
  	vector< vector< vector<double> > > r2_vecvecvec(n_channels);
  	vector< vector< vector<double> > > fitted_elev_vecvecvec(n_channels);
  	//vector< vector< vector<int> > > these_segment_lengths_vecvec(n_channels);

  	// iterators to navigate these data elements
  	vector< vector< vector<double> > >::iterator first_level_doub_iter;
  	vector< vector< vector<double> > >::iterator second_level_doub_iter;

	// now expand all of these vecvecvecs to be the correct size
	for (int cn = 0; cn<n_channels; cn++)
	{
		vector< vector<double> > temp_vecvec_double;
		//vector< vector<int> > temp_vecvec_int;
		vector<double> empty_double_vec;
		//vector<int> empty_int_vec;

		// expand the vecvec elements
		int nodes_in_channel = int(chis[cn].size());
		for (int n = 0; n<nodes_in_channel; n++)
		{
			temp_vecvec_double.push_back(empty_double_vec);
			//temp_vecvec_int.push_back(empty_int_vec);
		}

		// now add this to the vecvecvecs
		b_vecvecvec[cn] = temp_vecvec_double;
		m_vecvecvec[cn]= temp_vecvec_double;
		DW_vecvecvec[cn] = temp_vecvec_double;
		r2_vecvecvec[cn] =temp_vecvec_double;
		fitted_elev_vecvecvec[cn] = temp_vecvec_double;
		//these_segment_lengths_vecvec[cn] = temp_vecvec_int;
	}

	// now the vectors that will be replaced by the fitting algorithm
  	// theyare from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;
	vector<int> node_reference;
	int n_segments;

	cout << "LSDChiNetwork::sample after breaks, iteration: ";

	// loop through the iterations
	for (int it = 1; it <= n_iterations; it++)
	{
		if (it%10 == 0)
		{
			cout << " " << it;
		}

		//cout << "LSDChiNetwork::monte_carlo_sample_river_network_for_best_fit_after_breaks, iteration: " << it << endl;

      	// now loop through channels
      	for (int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << endl << "LINE 2560 LSDChiNetwork channel: " << chan << endl;
			vector<int> break_nodes = break_nodes_vecvec[chan];

			//for (int j = 0; j< int(break_nodes.size()); j++)
			//{
			//	cout << "break["<<j<<"]: " << break_nodes[j] << endl;
			//}

			// get the data from the channel
			vector<double> reverse_Chi = chis[chan];
			reverse(reverse_Chi.begin(), reverse_Chi.end());
			vector<double> reverse_Elevation = elevations[chan];
			reverse(reverse_Elevation.begin(), reverse_Elevation.end());

			// these data members keep track of the  breaks
			//int n_nodes = int(reverse_Chi.size());
			vector<int>::iterator br_iter;
			vector<double>::iterator vec_iter_start;
			vector<double>::iterator vec_iter_end;

			int start_of_last_break = 0;
			int length_of_break_segment;

			// now loop though breaks, getting the best fit.
			br_iter = break_nodes.begin();
			while(br_iter != break_nodes.end())
			{
				//cout << "starting break: " << start_of_last_break << " and this break: " << (*br_iter) << endl;
				length_of_break_segment = (*br_iter)-start_of_last_break+1;
				//cout << "length of break segement: " << length_of_break_segment << endl;

				// get the data of this break
				vec_iter_start = reverse_Chi.begin()+start_of_last_break;
				vec_iter_end = reverse_Chi.begin()+length_of_break_segment+start_of_last_break;
				vector<double> br_chi;
				br_chi.assign(vec_iter_start,vec_iter_end);

				vec_iter_start = reverse_Elevation.begin()+start_of_last_break;
				vec_iter_end = reverse_Elevation.begin()+length_of_break_segment+start_of_last_break;
				vector<double> br_elev;
				br_elev.assign(vec_iter_start,vec_iter_end);

				// these vectors hold the thinned data
				vector<double> m_per_node;
				vector<double> b_per_node;
				vector<double> r2_per_node;
				vector<double> DW_per_node;

				// Run the monte carlo algorithm on the break segment
				LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, br_chi, br_elev);

				// now thin the data, preserving the data (not interpolating)
				channel_MLE_finder.thin_data_monte_carlo_skip(skip, skip_range, node_reference);
				n_data_nodes = node_reference.size();
				//cout << "n_data_nodes after skip" << n_data_nodes << " and before: " << br_chi.size() << endl;

				//for (int k = 0; k<n_data_nodes; k++)
				//{
				//	cout << "nr["<<k<<"]: " << node_reference[k] << endl;
				//}

				// now create a single sigma value vector
				vector<double> sigma_values;
				sigma_values.push_back(sigma);

				// compute the best fit AIC
				channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

				// get the segments
				channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
										r2_vec, DW_vec, elev_fitted,these_segment_lengths,
										this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

				n_segments = int(b_vec.size());
				for(int seg = 0; seg<n_segments; seg++)
				{
					//cout << "segment length: " << these_segment_lengths[seg] << endl;
					for(int n = 0; n < these_segment_lengths[seg]; n++)
					{
						m_per_node.push_back(m_vec[seg]);
						b_per_node.push_back(b_vec[seg]);
						r2_per_node.push_back(r2_vec[seg]);
						DW_per_node.push_back(DW_vec[seg]);
					}
				}

				//cout << "n_data_nodes " << n_data_nodes <<  endl;

				// now assign the values of these variable to the vecvecvecs
				for (int n = 0; n< n_data_nodes; n++)
				{

					int this_node = node_reference[n]+start_of_last_break;
  					b_vecvecvec[chan][this_node].push_back(b_per_node[n]);
  					m_vecvecvec[chan][this_node].push_back(m_per_node[n]);
  					DW_vecvecvec[chan][this_node].push_back(DW_per_node[n]);
  					r2_vecvecvec[chan][this_node].push_back(r2_per_node[n]);
  					fitted_elev_vecvecvec[chan][this_node].push_back(elev_fitted[n]);
				}

				// reset the starting node
				start_of_last_break = (*br_iter)+1;
				//cout << "start_of last break: " << start_of_last_break << endl;

				br_iter++;
			}	// end break loop

		}		// end channel loop
	}			// end iteration loop
	cout << endl;

	// reset the data holding the fitted network properties
	vector< vector<double> > empty_vecvec;
	vector< vector<int> > empty_int_vecvec;
	chi_m_means = empty_vecvec;
	chi_m_standard_deviations = empty_vecvec;
	chi_m_standard_errors = empty_vecvec;
	chi_b_means = empty_vecvec;
	chi_b_standard_deviations = empty_vecvec;
	chi_b_standard_errors = empty_vecvec;
	all_fitted_elev_means = empty_vecvec;
	all_fitted_elev_standard_deviations = empty_vecvec;
	all_fitted_elev_standard_errors = empty_vecvec;
	chi_DW_means = empty_vecvec;
	chi_DW_standard_deviations = empty_vecvec;
	chi_DW_standard_errors = empty_vecvec;
	n_data_points_used_in_stats = empty_int_vecvec;

	// vectors that accept the data from the vecvecvec and are passed to
	// the get_common_statistics function
	vector<double> b_datavec;
	vector<double> m_datavec;
	vector<double> DW_datavec;
	vector<double> elev_datavec;
	vector<double> common_stats;

	// now go through the channel nodes and see how many data elements there are.
	for (int chan = 0; chan<n_channels; chan++)
	{
		// initialize vectors for storing the statistics of the
		// monte carlo fitted network properties
		int n_nodes_in_chan = int(chis[chan].size());

		vector<double> m_means(n_nodes_in_chan);
		vector<double> m_standard_deviations(n_nodes_in_chan);
		vector<double> m_standard_error(n_nodes_in_chan);
		vector<double> b_means(n_nodes_in_chan);
		vector<double> b_standard_deviations(n_nodes_in_chan);
		vector<double> b_standard_error(n_nodes_in_chan);
		vector<double> DW_means(n_nodes_in_chan);
		vector<double> DW_standard_deviations(n_nodes_in_chan);
		vector<double> DW_standard_error(n_nodes_in_chan);
		vector<double> fitted_elev_means(n_nodes_in_chan);
		vector<double> fitted_elev_standard_deviations(n_nodes_in_chan);
		vector<double> fitted_elev_standard_error(n_nodes_in_chan);
		vector<int> n_data_points_in_this_channel_node(n_nodes_in_chan);

		// now loop through each node, calculating how many data points there are in each
		for (int n = 0; n< n_nodes_in_chan; n++)
		{
			// get the number of data points in this channel.
			int n_data_points_itc = int(b_vecvecvec[chan][n].size());
			n_data_points_in_this_channel_node[n] = n_data_points_itc;
			b_datavec = b_vecvecvec[chan][n];
			m_datavec = m_vecvecvec[chan][n];
			DW_datavec = DW_vecvecvec[chan][n];
			elev_datavec = fitted_elev_vecvecvec[chan][n];

			// calcualte statistics, but only if there is data
			if (n_data_points_itc > 0)
			{
				common_stats = get_common_statistics(b_datavec);
				b_means[n] = common_stats[0];
				b_standard_deviations[n] = common_stats[2];
				b_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(m_datavec);
				m_means[n] = common_stats[0];
				m_standard_deviations[n] = common_stats[2];
				m_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(DW_datavec);
				DW_means[n] = common_stats[0];
				DW_standard_deviations[n] = common_stats[2];
				DW_standard_error[n] = common_stats[3];

				common_stats = get_common_statistics(elev_datavec);
				fitted_elev_means[n] = common_stats[0];
				fitted_elev_standard_deviations[n] = common_stats[2];
				fitted_elev_standard_error[n] = common_stats[3];
			}
			else	// if there is no data, use the previous node. This works because there is always data in the 1st node
			{
				b_means[n] = b_means[n-1];
				b_standard_deviations[n] = b_standard_deviations[n-1];
				b_standard_error[n] = b_standard_error[n-1];

				m_means[n] = m_means[n-1];
				m_standard_deviations[n] = m_standard_deviations[n-1];
				m_standard_error[n] = m_standard_error[n-1];

				DW_means[n] = DW_means[n-1];
				DW_standard_deviations[n] = DW_standard_deviations[n-1];
				DW_standard_error[n] = DW_standard_error[n-1];

				fitted_elev_means[n] = fitted_elev_means[n-1];
				fitted_elev_standard_deviations[n] = fitted_elev_standard_deviations[n-1];
				fitted_elev_standard_error[n] = fitted_elev_standard_error[n-1];
			}

		}	// end looping through channel nodes

		// all of the data vectors get reversed by the data fitting algorithm so they need to
		// be reinverted before they get inserted into the vecvecs
		reverse(m_means.begin(), m_means.end());
		reverse(m_standard_deviations.begin(), m_standard_deviations.end());
		reverse(m_standard_error.begin(), m_standard_error.end());

		reverse(b_means.begin(), b_means.end());
		reverse(b_standard_deviations.begin(), b_standard_deviations.end());
		reverse(b_standard_error.begin(), b_standard_error.end());

		reverse(DW_means.begin(), DW_means.end());
		reverse(DW_standard_deviations.begin(), DW_standard_deviations.end());
		reverse(DW_standard_error.begin(), DW_standard_error.end());

		reverse(fitted_elev_means.begin(), fitted_elev_means.end());
		reverse(fitted_elev_standard_deviations.begin(), fitted_elev_standard_deviations.end());
		reverse(fitted_elev_standard_error.begin(), fitted_elev_standard_error.end());

		reverse(n_data_points_in_this_channel_node.begin(),n_data_points_in_this_channel_node.end());

		// now store all the data
		chi_m_means.push_back(m_means);
		chi_m_standard_deviations.push_back(m_standard_deviations);
		chi_m_standard_errors.push_back(m_standard_error);
		chi_b_means.push_back(b_means);
		chi_b_standard_deviations.push_back(b_standard_deviations);
		chi_b_standard_errors.push_back(b_standard_error);
		chi_DW_means.push_back(DW_means);
		chi_DW_standard_deviations.push_back(DW_standard_deviations);
		chi_DW_standard_errors.push_back(DW_standard_error);
		all_fitted_elev_means.push_back(fitted_elev_means);
		all_fitted_elev_standard_deviations.push_back(fitted_elev_standard_deviations);
		all_fitted_elev_standard_errors.push_back(fitted_elev_standard_error);
		n_data_points_used_in_stats.push_back(n_data_points_in_this_channel_node);

	}				// end channel loop for processing data
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function uses the segment fitting tool to look for the best fit values of m over n
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_dchi(double A_0, int n_movern, double d_movern, double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes_mainstem, string fname)
{
  	double m_over_n;
  	int n_channels = chis.size();
  	double dchi;					// spacing of the chi vector. This is calcuated to be optimal for the
  									// main stem

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for individual channels
  	vector< vector<double> > b_vecvec(n_channels);
  	vector< vector<double> > m_vecvec(n_channels);
  	vector< vector<double> > DW_vecvec(n_channels);
  	vector< vector<double> > r2_vecvec(n_channels);
  	vector< vector<double> > thinned_chi_vecvec(n_channels);
  	vector< vector<double> > thinned_elev_vecvec(n_channels);
  	vector< vector<double> > fitted_elev_vecvec(n_channels);
  	vector< vector<int> > node_ref_thinned_vecvec(n_channels);
  	vector< vector<int> > these_segment_lengths_vecvec(n_channels);
  	vector<double> MLE_vec(n_channels);
  	vector<int> n_segments_vec(n_channels);
  	vector<int> n_data_nodes_vec(n_channels);
  	vector<double> AICc_vec(n_channels, 9999);
  	vector<double> best_m_over_n(n_channels);

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for the cumulative channels
  	vector< vector<double> > cum_b_vecvec(n_channels);
  	vector< vector<double> > cum_m_vecvec(n_channels);
  	vector< vector<double> > cum_DW_vecvec(n_channels);
  	vector< vector<double> > cum_r2_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_chi_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_elev_vecvec(n_channels);
  	vector< vector<double> > cum_fitted_elev_vecvec(n_channels);
  	vector< vector<int> > cum_node_ref_thinned_vecvec(n_channels);
  	vector<vector<int> > cum_these_segment_lengths_vecvec(n_channels);

  	// these data are for the cumulative AICs
  	vector<double> AICc_combined_vec(n_movern);
  	vector<double> m_over_n_vec(n_movern);

  	// these are from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;

  	// loop through m_over_n values
	for(int movn = 0; movn< n_movern; movn++)
  	{
		m_over_n = double(movn)*d_movern+start_movern;
	  	cout << "m/n: " << m_over_n << endl;

      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);
  		dchi = calculate_optimal_chi_spacing(target_nodes_mainstem);

   		vector<double> MLEs_thischan(n_channels);
		vector<int> n_segs_thischan(n_channels);
      	vector<int> n_datanodes_thischan(n_channels);

     	// now loop through channels
      	for (int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << "LINE439 LSDChiNetwork channel: " << chan << endl;
		  	// get the channels for this m over n ratio
	  	  	find_most_likeley_segments_dchi(chan, minimum_segment_length, sigma, dchi,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
							      this_AIC, this_AICc );
			//cout << "LINE443 LSDChiNetwork channel: " << chan << endl;

			// check to see if the AICc value is the smallest
	  		// if so add the data to the best fit data elements
	  		if (this_AICc < AICc_vec[chan])
	    	{
	       		b_vecvec[chan] = b_vec;
	       		m_vecvec[chan] = m_vec;
	       		DW_vecvec[chan] = DW_vec;
	       		r2_vecvec[chan] = r2_vec;
	       		thinned_chi_vecvec[chan] = chi_thinned;
	       		thinned_elev_vecvec[chan] = elev_thinned;
	       		fitted_elev_vecvec[chan] = elev_fitted;
	       		node_ref_thinned_vecvec[chan] = node_ref_thinned;
	       		these_segment_lengths_vecvec[chan] = these_segment_lengths;
	       		MLE_vec[chan] = this_MLE;
	       		n_segments_vec[chan] = this_n_segments;
	       		n_data_nodes_vec[chan] = n_data_nodes;
	       		AICc_vec[chan] = this_AICc;
	       		best_m_over_n[chan] = m_over_n;
			}

	  		// add the data from this channel to the vectors that will be used to calcualte cumulative AICc
	  		MLEs_thischan[chan] = this_MLE;
	  		n_segs_thischan[chan] = this_n_segments;
	  		n_datanodes_thischan[chan] = n_data_nodes;
		}

      	//now calculate the cumulative AICc for this m over n
      	double thismn_AIC;
      	double thismn_AICc;

      	int n_total_segments = 0;
      	int n_total_nodes = 0;
      	double cumulative_MLE = 1;
      	double log_cum_MLE = 0;

		// get the cumulative maximum likelihood estimators
      	for (int chan = 0; chan<n_channels; chan++)
		{
			n_total_segments += n_segs_thischan[chan];
	  		n_total_nodes += n_datanodes_thischan[chan];
	  		cumulative_MLE = MLEs_thischan[chan]*cumulative_MLE;

			if(MLEs_thischan[chan] <= 0)
			{
				log_cum_MLE = log_cum_MLE-1000;
			}
			else
			{
				log_cum_MLE = log(MLEs_thischan[chan])+log_cum_MLE;
			}

		}


      	// these AIC and AICc values are cumulative for a given m_over_n
      	//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
			                                                       // for each segment there are 2 parameters
      	thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
      	thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);
      	AICc_combined_vec[movn] = thismn_AICc;
      	m_over_n_vec[movn] = m_over_n;

      	//cout << "m_over_n: " << m_over_n << " and combined AICc: " << thismn_AICc << endl;
      	//cout << "this cumulative MLE: " << cumulative_MLE << " n_segs: " << n_total_segments << " and n_nodes: " << n_total_nodes << endl;

    }

    cout << "and the cumulative m_over n values"<< endl;
    double min_cum_AICc = 9999;
    double bf_cum_movn = start_movern;
    for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		cout << "m over n: " << m_over_n_vec[mn] << " and AICc: " << AICc_combined_vec[mn] << endl;
		// if this is the minimum, store the m over n value
		if(AICc_combined_vec[mn] < min_cum_AICc)
	  	{
	   	 	min_cum_AICc = AICc_combined_vec[mn];
	    	bf_cum_movn = m_over_n_vec[mn];
	  	}
	}
	cout << "LSDChiNetwork Line 537 so the best fit m over n is: " << bf_cum_movn << endl;

    // now get the cumulative best fit channels
    calculate_chi(A_0, bf_cum_movn);
    dchi = calculate_optimal_chi_spacing(target_nodes_mainstem);
    cout << "LSDChiNetwork Line 542 dchi is: " << dchi << endl;
    // now loop through channels
    for (int chan = 0; chan<n_channels; chan++)
    {

		// get the channels for this m over n ratio
		find_most_likeley_segments_dchi(chan, minimum_segment_length, sigma, dchi,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
                                            this_AIC, this_AICc );
	 	cum_b_vecvec[chan] = b_vec;
	 	cum_m_vecvec[chan] = m_vec;
        cum_DW_vecvec[chan] = DW_vec;
        cum_r2_vecvec[chan] = r2_vec;
        cum_thinned_chi_vecvec[chan] = chi_thinned;
        cum_thinned_elev_vecvec[chan] = elev_thinned;
        cum_fitted_elev_vecvec[chan] = elev_fitted;
        cum_node_ref_thinned_vecvec[chan] = node_ref_thinned;
        cum_these_segment_lengths_vecvec[chan] = these_segment_lengths;
	}

    // write a file
    ofstream best_fit_info;
	best_fit_info.open(fname.c_str());

  	best_fit_info << "N_channels: " << n_channels << endl;
  	best_fit_info << "m_over_n_for_channels ";
    for (int ch = 0; ch<n_channels; ch++)
    {
		best_fit_info << "  " << best_m_over_n[ch];
    }
  	best_fit_info << endl << "m_over_n_values ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << m_over_n_vec[mn];
	}
	best_fit_info << endl << "cumulative_AICc: ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << AICc_combined_vec[mn];
	}
	best_fit_info << endl;
	cout << "The_best_fit_cumulative_m_over_n_is: " << bf_cum_movn << endl;
	for (int chan = 0; chan<n_channels; chan++)
	{
		vector<int> seglength = cum_these_segment_lengths_vecvec[chan];
		vector<double> m_val = cum_m_vecvec[chan];
		vector<double> b_val = cum_b_vecvec[chan];
		vector<double> DW_val = cum_DW_vecvec[chan];
		vector<double> r2_val = cum_r2_vecvec[chan];
		int n_segs_this_channel = seglength.size();

		best_fit_info << "Channel " << chan << " segment_length";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << seglength[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_gradient";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << m_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_intercept";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << b_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_DW_stat";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << DW_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_r2";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << r2_val[i];
		}
		best_fit_info << endl;
	}

	best_fit_info << endl << endl << "Cumulative_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = cum_thinned_chi_vecvec[chan];
		elev_thinned = cum_thinned_elev_vecvec[chan];
		elev_fitted = cum_fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }


	best_fit_info << endl << endl << "Individually_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = thinned_chi_vecvec[chan];
		elev_thinned = thinned_elev_vecvec[chan];
		elev_fitted = fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }
     return bf_cum_movn;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function uses the segment fitting tool to look for the best fit values of m over n
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n(double A_0, int n_movern, double d_movern, double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes_mainstem, string fname)
{
  	double m_over_n;
  	int n_channels = chis.size();

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for individual channels
  	vector< vector<double> > b_vecvec(n_channels);
  	vector< vector<double> > m_vecvec(n_channels);
  	vector< vector<double> > DW_vecvec(n_channels);
  	vector< vector<double> > r2_vecvec(n_channels);
  	vector< vector<double> > thinned_chi_vecvec(n_channels);
  	vector< vector<double> > thinned_elev_vecvec(n_channels);
  	vector< vector<double> > fitted_elev_vecvec(n_channels);
  	vector< vector<int> > node_ref_thinned_vecvec(n_channels);
  	vector< vector<int> > these_segment_lengths_vecvec(n_channels);
  	vector<double> MLE_vec(n_channels);
  	vector<int> n_segments_vec(n_channels);
  	vector<int> n_data_nodes_vec(n_channels);
  	vector<double> AICc_vec(n_channels, 9999);
  	vector<double> best_m_over_n(n_channels);

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for the cumulative channels
  	vector< vector<double> > cum_b_vecvec(n_channels);
  	vector< vector<double> > cum_m_vecvec(n_channels);
  	vector< vector<double> > cum_DW_vecvec(n_channels);
  	vector< vector<double> > cum_r2_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_chi_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_elev_vecvec(n_channels);
  	vector< vector<double> > cum_fitted_elev_vecvec(n_channels);
  	vector< vector<int> > cum_node_ref_thinned_vecvec(n_channels);
  	vector<vector<int> > cum_these_segment_lengths_vecvec(n_channels);

  	// these data are for the cumulative AICs
  	vector<double> AICc_combined_vec(n_movern);
  	vector<double> m_over_n_vec(n_movern);

  	// these are from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;

	int N = calculate_skip(target_nodes_mainstem);

  	// loop through m_over_n values
	for(int movn = 0; movn< n_movern; movn++)
  	{
		m_over_n = double(movn)*d_movern+start_movern;
	  	cout << "m/n: " << m_over_n << endl;

      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);

   		vector<double> MLEs_thischan(n_channels);
		vector<int> n_segs_thischan(n_channels);
      	vector<int> n_datanodes_thischan(n_channels);

     	// now loop through channels
      	for (int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << "LINE439 LSDChiNetwork channel: " << chan << endl;
		  	// get the channels for this m over n ratio
	  	  	find_most_likeley_segments(chan, minimum_segment_length, sigma, N,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
							      this_AIC, this_AICc );
			//cout << "LINE443 LSDChiNetwork channel: " << chan << endl;

			// check to see if the AICc value is the smallest
	  		// if so add the data to the best fit data elements
	  		if (this_AICc < AICc_vec[chan])
	    	{
	       		b_vecvec[chan] = b_vec;
	       		m_vecvec[chan] = m_vec;
	       		DW_vecvec[chan] = DW_vec;
	       		r2_vecvec[chan] = r2_vec;
	       		thinned_chi_vecvec[chan] = chi_thinned;
	       		thinned_elev_vecvec[chan] = elev_thinned;
	       		fitted_elev_vecvec[chan] = elev_fitted;
	       		node_ref_thinned_vecvec[chan] = node_ref_thinned;
	       		these_segment_lengths_vecvec[chan] = these_segment_lengths;
	       		MLE_vec[chan] = this_MLE;
	       		n_segments_vec[chan] = this_n_segments;
	       		n_data_nodes_vec[chan] = n_data_nodes;
	       		AICc_vec[chan] = this_AICc;
	       		best_m_over_n[chan] = m_over_n;
			}

	  		// add the data from this channel to the vectors that will be used to calcualte cumulative AICc
	  		MLEs_thischan[chan] = this_MLE;
	  		n_segs_thischan[chan] = this_n_segments;
	  		n_datanodes_thischan[chan] = n_data_nodes;
		}

      	//now calculate the cumulative AICc for this m over n
      	double thismn_AIC;
      	double thismn_AICc;

      	int n_total_segments = 0;
      	int n_total_nodes = 0;
      	double cumulative_MLE = 1;
      	double log_cum_MLE = 0;

		// get the cumulative maximum likelihood estimators
      	for (int chan = 0; chan<n_channels; chan++)
		{
			// test if the segment is too short
			if (n_datanodes_thischan[chan] > 3*minimum_segment_length)
			{
				n_total_segments += n_segs_thischan[chan];
	  			n_total_nodes += n_datanodes_thischan[chan];
	  			cumulative_MLE = MLEs_thischan[chan]*cumulative_MLE;

				if(MLEs_thischan[chan] <= 0)
				{
					log_cum_MLE = log_cum_MLE-1000;
				}
				else
				{
					log_cum_MLE = log(MLEs_thischan[chan])+log_cum_MLE;
				}

			}
			//else
			//{
			//	cout << "channel " << chan << " rejected; too few nodes" << endl;
			//}
		}


      	// these AIC and AICc values are cumulative for a given m_over_n
      	//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
			                                                       // for each segment there are 2 parameters
      	thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
      	thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);
      	AICc_combined_vec[movn] = thismn_AICc;
      	m_over_n_vec[movn] = m_over_n;

      	//cout << "m_over_n: " << m_over_n << " and combined AICc: " << thismn_AICc << endl;
      	//cout << "this cumulative MLE: " << cumulative_MLE << " n_segs: " << n_total_segments << " and n_nodes: " << n_total_nodes << endl;

    }

    cout << "and the cumulative m_over n values"<< endl;
    double min_cum_AICc = 9999;
    double bf_cum_movn = start_movern;
    for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		cout << "m over n: " << m_over_n_vec[mn] << " and AICc: " << AICc_combined_vec[mn] << endl;
		// if this is the minimum, store the m over n value
		if(AICc_combined_vec[mn] < min_cum_AICc)
	  	{
	   	 	min_cum_AICc = AICc_combined_vec[mn];
	    	bf_cum_movn = m_over_n_vec[mn];
	  	}
	}
	cout << "LSDChiNetwork Line 537 so the best fit m over n is: " << bf_cum_movn << endl;

    // now get the cumulative best fit channels
    calculate_chi(A_0, bf_cum_movn);
    // now loop through channels
    for (int chan = 0; chan<n_channels; chan++)
    {

		// get the channels for this m over n ratio
		find_most_likeley_segments(chan, minimum_segment_length, sigma, N,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
                                            this_AIC, this_AICc );
	 	cum_b_vecvec[chan] = b_vec;
	 	cum_m_vecvec[chan] = m_vec;
        cum_DW_vecvec[chan] = DW_vec;
        cum_r2_vecvec[chan] = r2_vec;
        cum_thinned_chi_vecvec[chan] = chi_thinned;
        cum_thinned_elev_vecvec[chan] = elev_thinned;
        cum_fitted_elev_vecvec[chan] = elev_fitted;
        cum_node_ref_thinned_vecvec[chan] = node_ref_thinned;
        cum_these_segment_lengths_vecvec[chan] = these_segment_lengths;
	}

    // write a file
    ofstream best_fit_info;
	best_fit_info.open(fname.c_str());

  	best_fit_info << "N_channels: " << n_channels << endl;
  	best_fit_info << "m_over_n_for_channels ";
    for (int ch = 0; ch<n_channels; ch++)
    {
		best_fit_info << "  " << best_m_over_n[ch];
    }
  	best_fit_info << endl << "m_over_n_values ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << m_over_n_vec[mn];
	}
	best_fit_info << endl << "cumulative_AICc: ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << AICc_combined_vec[mn];
	}
	best_fit_info << endl;
	cout << "The_best_fit_cumulative_m_over_n_is: " << bf_cum_movn << endl;
	for (int chan = 0; chan<n_channels; chan++)
	{
		vector<int> seglength = cum_these_segment_lengths_vecvec[chan];
		vector<double> m_val = cum_m_vecvec[chan];
		vector<double> b_val = cum_b_vecvec[chan];
		vector<double> DW_val = cum_DW_vecvec[chan];
		vector<double> r2_val = cum_r2_vecvec[chan];
		int n_segs_this_channel = seglength.size();

		best_fit_info << "Channel " << chan << " segment_length";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << seglength[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_gradient";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << m_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_intercept";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << b_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_DW_stat";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << DW_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_r2";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << r2_val[i];
		}
		best_fit_info << endl;
	}

	best_fit_info << endl << endl << "Cumulative_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = cum_thinned_chi_vecvec[chan];
		elev_thinned = cum_thinned_elev_vecvec[chan];
		elev_fitted = cum_fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }


	best_fit_info << endl << endl << "Individually_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = thinned_chi_vecvec[chan];
		elev_thinned = thinned_elev_vecvec[chan];
		elev_fitted = fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }
     return bf_cum_movn;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function looks for the best fit values of m over n by simply testing for the least variation
// in the tributaries
//
// SMM 01/05/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_colinearity_test(double A_0, int n_movern, double d_movern,
						       double start_movern, int minimum_segment_length, double sigma,
						       int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values, vector<double>& AICc_mean, vector<double>& AICc_sdtd)
{

	cout << "starting colinearity search" << endl;

  	double m_over_n;
  	int n_channels = chis.size();

	// these vectors contain all the chi and elevation values
  	vector<double> compiled_chis;
  	vector<double> compiled_elev;
  	vector<double> sorted_chis;
  	vector<double> sorted_elev;
  	vector<double> empty_vec;

  	vector<double> this_chi;
  	vector<double> this_elev;

  	// vector for storing the sorting index
  	vector<size_t> index_map;

  	// vector for holding the AICc values and the AICc variation
  	vector<double> AICc_for_mover_n(n_movern);
  	vector<double> AICc_std_dev(n_movern);
  	vector<double> movn_values(n_movern);

	// these are data elements used by the segment finder
	vector<double> b_vec;
	vector<double> m_vec;
	vector<double> r2_vec;
	vector<double> DW_vec;
	vector<double> fitted_elev;
	vector<int> these_segment_lengths;
	double this_MLE;
	int this_n_segments;
	int n_data_nodes;
	double this_AIC;
	double this_AICc;

	cout << "looping starting line 4280" << endl;


  	// loop through m over n
	for(int movn = 0; movn< n_movern; movn++)
  	{
		//cout << "i: " << movn << endl;
		//cout << "d_movn: " << d_movern << " start: " << start_movern << endl;
		m_over_n = double(movn)*d_movern+start_movern;
		//cout << "yo, : " << double(movn)*d_movern+start_movern << endl;

		movn_values[movn] = m_over_n;
		cout << "m over n: " << movn_values[movn];

		// reset the compiled vectors
		compiled_chis = empty_vec;
		compiled_elev = empty_vec;
		sorted_chis = empty_vec;
		sorted_elev = empty_vec;

      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);

	  	// now load all the chis and elevations into individual vectors
	  	for(int chan = 0; chan<n_channels; chan++)
	  	{
			this_chi = chis[chan];
			this_elev = elevations[chan];

			// add the cis and elevations to the compiled vectors
			int n_nodes = int(this_chi.size());
			for (int node = 0; node<n_nodes; node++)
			{
				compiled_chis.push_back(this_chi[node]);
				compiled_elev.push_back(this_elev[node]);
			}
		}

		// now sort the vectors
		matlab_double_sort(compiled_chis, sorted_chis, index_map);
		matlab_double_reorder(compiled_elev, index_map, sorted_elev);

		// now we use the monte carlo sampling to look for the best fit
		vector<int> empty_vec;

		// reverse the vectors
		vector<double> reverse_Chi = sorted_chis;
		reverse(reverse_Chi.begin(), reverse_Chi.end());
		vector<double> reverse_Elevation = sorted_elev;
		reverse(reverse_Elevation.begin(), reverse_Elevation.end());

		// calculate the baseline skipping for monte carlo analyis
		int mean_skip = calculate_skip(target_nodes, reverse_Chi);
		int skip_range = mean_skip*2;
		if (skip_range ==0)
		{
			skip_range = 2;
		}
		if (skip_range < 0)
		{
			skip_range = -skip_range;
		}

		vector<double> these_AICcs;
		cout << " iteration: ";
		// now enter the iterative phase
		for (int iter = 0; iter<n_iterations; iter++)
		{

			if (iter%20 == 0)
			{
				cout << iter << " ";
			}

			LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

			vector<int> node_reference;
			// now thin the data, preserving the data (not interpolating)
			channel_MLE_finder.thin_data_monte_carlo_skip(mean_skip, skip_range, node_reference);
			//cout << "The thinned number of nodes is: " << node_reference.size() << " and overall nodes: " << reverse_Chi.size() << endl;

			// now create a single sigma value vector
			vector<double> sigma_values;
			sigma_values.push_back(sigma);

			// compute the best fit AIC
			channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);

			// return the data
			channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
								r2_vec, DW_vec, fitted_elev,these_segment_lengths,
								this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

			// now get the AICc values and intert them into the these AICs vector
			these_AICcs.push_back(this_AICc);
		}

		// now get the mean and standard deviation for this m_over_n
		AICc_for_mover_n[movn] = get_mean(these_AICcs);
		AICc_std_dev[movn] = get_standard_deviation(these_AICcs, AICc_for_mover_n[movn]);

		cout << endl;

	}

	// now calucalte the best fit m over n
	double bf_movern = start_movern;
	double bf_AICc = 100000;
	for (int i = 0; i<n_movern; i++)
	{

		if(AICc_for_mover_n[i] < bf_AICc)
		{
			bf_AICc= AICc_for_mover_n[i];
			bf_movern = movn_values[i];
		}
	}

	m_over_n_values = movn_values;
	AICc_mean = AICc_for_mover_n;
	AICc_sdtd = AICc_std_dev;

	return bf_movern;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function looks for the best fit values of m over n by simply testing for the least variation
// in the tributaries
//
// If Monte_Carlo_switch = 1, then the routine iterates on the AICc and returns mean and std_deviation
// of the AICc. If not it simply returns the single AICc value and the std_dev vector just holds
// a repeat of the mean value
//
// SMM 15/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_colinearity_test_with_breaks(double A_0, int n_movern, double d_movern,
						       double start_movern, int minimum_segment_length, double sigma,
						       int target_skip, int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values, vector<double>& AICc_mean, vector<double>& AICc_sdtd,
						       int Monte_Carlo_switch)
{

	cout << "LINE 3972 starting colinearity search using breaks" << endl;

  	double m_over_n;
  	int n_channels = chis.size();

	// these vectors contain all the chi and elevation values
  	vector<double> compiled_chis;
  	vector<double> compiled_elev;
  	vector<double> sorted_chis;
  	vector<double> sorted_elev;
  	vector<double> empty_vec;

  	vector<double> this_chi;
  	vector<double> this_elev;

  	// vector for storing the sorting index
  	vector<size_t> index_map;

  	// vector for holding the AICc values and the AICc variation
  	vector<double> AICc_for_mover_n(n_movern);
  	vector<double> AICc_std_dev(n_movern);
  	vector<double> movn_values(n_movern, start_movern);

	// these are data elements used by the segment finder
	vector<double> b_vec;
	vector<double> m_vec;
	vector<double> r2_vec;
	vector<double> DW_vec;
	vector<double> fitted_elev;
	vector<int> these_segment_lengths;
	//double this_MLE;
	//int this_n_segments;
	//int n_data_nodes;
	double this_AICc;

	cout << "looping line 4456" << endl;


  	// loop through m over n
	for(int movn = 0; movn< n_movern; movn++)
  	{
		//cout << "i: " << movn << endl;
		//cout << "d_movn: " << d_movern << " start: " << start_movern << endl;
		m_over_n = double(movn)*d_movern+start_movern;
		//cout << "yo, : " << double(movn)*d_movern+start_movern << endl;

		cout << "n_movern: " << n_movern << " sz: " << movn_values.size() << " and m_over_n: " << m_over_n <<  endl;

		///for (int i = 0; i< int(movn_values.size()); i++)
		//{
		//	cout << "movn_values["<<i<<"]: " << movn_values[i] << endl;
		//}

		movn_values[movn] = m_over_n;
		//cout << "m over n: " << movn_values[movn] << endl;

		// reset the compiled vectors
		compiled_chis = empty_vec;
		compiled_elev = empty_vec;
		sorted_chis = empty_vec;
		sorted_elev = empty_vec;

		//cout << "LINE 4483 reset_vectors" << endl;

      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);

	  	//cout << "LINE 4488 calculated_chi" << endl;

	  	// now load all the chis and elevations into individual vectors
	  	for(int chan = 0; chan<n_channels; chan++)
	  	{
			//cout << "LINE 4493 chan: " << chan << endl;

			this_chi = chis[chan];
			this_elev = elevations[chan];

			// add the chis and elevations to the compiled vectors
			int n_nodes = int(this_chi.size());
			for (int node = 0; node<n_nodes; node++)
			{
				compiled_chis.push_back(this_chi[node]);
				compiled_elev.push_back(this_elev[node]);
			}
		}

		//cout << "LINE 4507 sorting vectors" << endl;

		// now sort the vectors
		matlab_double_sort(compiled_chis, sorted_chis, index_map);
		matlab_double_reorder(compiled_elev, index_map, sorted_elev);

		// now we use the monte carlo sampling to look for the best fit
		vector<int> empty_vec;

		// reverse the vectors
		vector<double> reverse_Chi = sorted_chis;
		reverse(reverse_Chi.begin(), reverse_Chi.end());
		vector<double> reverse_Elevation = sorted_elev;
		reverse(reverse_Elevation.begin(), reverse_Elevation.end());

		//cout << "LINE 4522 reversed vectors" << endl;

		vector<double> these_AICcs;
		vector<int> break_nodes;

		int n_total_segments;
		int n_total_nodes;
		double cumulative_MLE;

		// now enter the iterative phase
		monte_carlo_split_channel_colinear(A_0, m_over_n, n_iterations,
				target_skip, target_nodes, minimum_segment_length, sigma,
				reverse_Chi, reverse_Elevation, break_nodes);

		//cout << "Line 4237, calculated breaks, break nodes are" << endl;
		//for (int i = 0; i< int(break_nodes.size()); i++)
		//{
		//	cout << break_nodes[i] << endl;
		//}

		//cout << "Line 4542, doing monte carlo AICc " << endl;
		if(Monte_Carlo_switch ==1)				// iterate: this takes a little while
		{
			these_AICcs = calculate_AICc_after_breaks_colinear_monte_carlo(A_0, m_over_n,
						target_skip, minimum_segment_length, sigma,
						reverse_Chi, reverse_Elevation, break_nodes,
						n_total_segments, n_total_nodes, cumulative_MLE,
						n_iterations);


			AICc_for_mover_n[movn] = get_mean(these_AICcs);
			AICc_std_dev[movn] =  get_standard_deviation(these_AICcs, AICc_for_mover_n[movn]);
		}
		else								// no iterating on teh AICc
		{
			this_AICc = calculate_AICc_after_breaks_colinear(A_0, m_over_n,
						target_skip, minimum_segment_length, sigma,
						reverse_Chi, reverse_Elevation, break_nodes,
						n_total_segments, n_total_nodes, cumulative_MLE);
								// now get the mean and standard deviation for this m_over_n
			AICc_for_mover_n[movn] = this_AICc;
			AICc_std_dev[movn] = this_AICc;
		}

		cout << endl;
	}

	// now calculate the best fit m over n
	double bf_movern = start_movern;
	double bf_AICc = 100000;
	for (int i = 0; i<n_movern; i++)
	{

		if(AICc_for_mover_n[i] < bf_AICc)
		{
			bf_AICc= AICc_for_mover_n[i];
			bf_movern = movn_values[i];
		}
	}
	//cout << "LINE 4581, got best fit m/n, colinear" << endl;

	m_over_n_values = movn_values;
	AICc_mean = AICc_for_mover_n;
	AICc_sdtd = AICc_std_dev;

	return bf_movern;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function calculeates best fit m/n for each channel
// these channels are ones with breaks
//
// SMM 15/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_individual_channels_with_breaks(double A_0, int n_movern, double d_movern,
						       double start_movern, int minimum_segment_length, double sigma,
						       int target_skip, int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values, vector< vector<double> >& AICc_vals)
{
	// get the number of channels:
	int n_chans = elevations.size();

	vector< vector<double> > AICc_local(n_chans);
	vector<double> chan_AICs(n_movern);
  	vector<double> m_over_n_vals(n_movern);

  	// this is used to compute a cumulative MLE (which is not at the moment used
  	int n_total_segments, n_total_nodes;
  	double cumulative_MLE;

	double AICc;
	double m_over_n;

	for(int chan = 0; chan<n_chans; chan++)
	{
		for(int movn = 0; movn< n_movern; movn++)
		{
			m_over_n = double(movn)*d_movern+start_movern;
			//cout << "start m/n: " << start_movern << " d m/m: " << d_movern << " moven: " << movn << " m/n: " << m_over_n << endl;

			m_over_n_vals[movn] = m_over_n;

			vector<int> break_nodes;

			monte_carlo_split_channel(A_0, m_over_n, n_iterations,
					target_skip, target_nodes, minimum_segment_length, sigma, chan, break_nodes);

			AICc = calculate_AICc_after_breaks(A_0, m_over_n, target_skip, minimum_segment_length, sigma, chan, break_nodes,
					n_total_segments, n_total_nodes, cumulative_MLE);

			chan_AICs[movn] = AICc;
		}
		AICc_local[chan] = chan_AICs;
	}

	AICc_vals = AICc_local;
	m_over_n_values = m_over_n_vals;

	// now get the best fit on the main stem and return this value
	chan_AICs = AICc_local[0];

	double best_AICc = 1000000;
	double bf_movern = start_movern;
	for(int i = 0; i<n_movern; i++)
	{
		if(chan_AICs[i] < best_AICc)
		{
			best_AICc = chan_AICs[i];
			bf_movern = m_over_n_vals[i];
		}
	}
	return bf_movern;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function calculeates best fit m/n for each channel
// these channels are ones with breaks
//
// this function uses a monte carlo sampling approach to give some idea of the variability and uncertanty in
// the AICc values
//
// SMM 15/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_individual_channels_with_breaks_monte_carlo(double A_0, int n_movern, double d_movern,
						       double start_movern, int minimum_segment_length, double sigma,
						       int target_skip, int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values,
						       vector< vector<double> >& AICc_means, vector< vector<double> >& AICc_stddev)
{
	// get the number of channels:
	int n_chans = elevations.size();

	vector< vector<double> > AICc_local(n_chans);
	vector< vector<double> > AICc_std(n_chans);
	vector<double> chan_AICs(n_movern);
	vector<double> chan_std(n_movern);
  	vector<double> m_over_n_vals(n_movern);

  	// this is used to compute a cumulative MLE (which is not at the moment used
  	int n_total_segments, n_total_nodes;
  	double cumulative_MLE;

	vector<double> these_AICcs;
	double m_over_n;

	for(int chan = 0; chan<n_chans; chan++)
	{
		for(int movn = 0; movn< n_movern; movn++)
		{
			m_over_n = double(movn)*d_movern+start_movern;
			//cout << "start m/n: " << start_movern << " d m/m: " << d_movern << " moven: " << movn << " m/n: " << m_over_n << endl;

			m_over_n_vals[movn] = m_over_n;

			vector<int> break_nodes;

			monte_carlo_split_channel(A_0, m_over_n, n_iterations,
					target_skip, target_nodes, minimum_segment_length, sigma, chan, break_nodes);

			cout << "LSDChiNet, m/n: " << m_over_n << " chan: " << chan;
			these_AICcs = calculate_AICc_after_breaks_monte_carlo(A_0, m_over_n, target_skip,
			                               minimum_segment_length, sigma, chan, break_nodes,
					                       n_total_segments, n_total_nodes, cumulative_MLE,
					                       n_iterations);

			chan_AICs[movn] = get_mean(these_AICcs);
			chan_std[movn] = get_standard_deviation(these_AICcs,chan_AICs[movn]);
		}
		AICc_local[chan] = chan_AICs;
		AICc_std[chan] = chan_std;
	}

	AICc_means = AICc_local;
	AICc_stddev = AICc_std;
	m_over_n_values = m_over_n_vals;

	// now get the best fit on the main stem and return this value
	chan_AICs = AICc_local[0];

	double best_AICc = 1000000;
	double bf_movern = start_movern;
	for(int i = 0; i<n_movern; i++)
	{
		if(chan_AICs[i] < best_AICc)
		{
			best_AICc = chan_AICs[i];
			bf_movern = m_over_n_vals[i];
		}
	}
	return bf_movern;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function uses the segment fitting tool to look for the best fit values of m over n
//
// SMM 15/06/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDChiNetwork::search_for_best_fit_m_over_n_seperate_ms_and_tribs(double A_0, int n_movern, double d_movern,
								double start_movern, int minimum_segment_length, double sigma,
						        int target_nodes_mainstem, string fname)
{
  	double m_over_n;
  	int n_channels = chis.size();

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for individual channels
  	vector< vector<double> > b_vecvec(n_channels);
  	vector< vector<double> > m_vecvec(n_channels);
  	vector< vector<double> > DW_vecvec(n_channels);
  	vector< vector<double> > r2_vecvec(n_channels);
  	vector< vector<double> > thinned_chi_vecvec(n_channels);
  	vector< vector<double> > thinned_elev_vecvec(n_channels);
  	vector< vector<double> > fitted_elev_vecvec(n_channels);
  	vector< vector<int> > node_ref_thinned_vecvec(n_channels);
  	vector< vector<int> > these_segment_lengths_vecvec(n_channels);
  	vector<double> MLE_vec(n_channels);
  	vector<int> n_segments_vec(n_channels);
  	vector<int> n_data_nodes_vec(n_channels);
  	vector<double> AICc_vec(n_channels, 9999);
  	vector<double> best_m_over_n(n_channels);

  	// data structures to hold information about segments from all channels
  	// these are the best fit data for the cumulative channels
  	vector< vector<double> > cum_b_vecvec(n_channels);
  	vector< vector<double> > cum_m_vecvec(n_channels);
  	vector< vector<double> > cum_DW_vecvec(n_channels);
  	vector< vector<double> > cum_r2_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_chi_vecvec(n_channels);
  	vector< vector<double> > cum_thinned_elev_vecvec(n_channels);
  	vector< vector<double> > cum_fitted_elev_vecvec(n_channels);
  	vector< vector<int> > cum_node_ref_thinned_vecvec(n_channels);
  	vector<vector<int> > cum_these_segment_lengths_vecvec(n_channels);

  	// these data are for the cumulative AICs
  	vector<double> AICc_combined_vec(n_movern);
  	vector<double> AICc_mainstem_vec(n_movern);
  	vector<double> m_over_n_vec(n_movern);

  	// these are from the individual channels, which are replaced each time a new channel is analyzed
  	vector<double> m_vec;
  	vector<double> b_vec;
  	vector<double> r2_vec;
  	vector<double> DW_vec;
  	vector<double> fitted_y;
  	int n_data_nodes;
  	int this_n_segments;
  	double this_MLE, this_AIC, this_AICc;
  	vector<int> these_segment_lengths;
  	vector<double> chi_thinned;
  	vector<double> elev_thinned;
  	vector<double> elev_fitted;
	vector<int> node_ref_thinned;

  	// loop through m_over_n values for the mainstem
	int N = calculate_skip(target_nodes_mainstem);
	int ms_N = N;
	cout << "LSDCN line 2470, ms N is: " << N << endl;
	for(int movn = 0; movn< n_movern; movn++)
  	{
		m_over_n = double(movn)*d_movern+start_movern;


      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);

		int chan = 0;
		// get the channels for this m over n ratio
		find_most_likeley_segments(chan, minimum_segment_length, sigma, N,
					b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					elev_fitted, node_ref_thinned,these_segment_lengths,
					this_MLE, this_n_segments, n_data_nodes,
							  this_AIC, this_AICc );
		//cout << "LINE443 LSDChiNetwork channel: " << chan << endl;

		// check to see if the AICc value is the smallest
		// if so add the data to the best fit data elements
		if (this_AICc < AICc_vec[chan])
		{
			b_vecvec[chan] = b_vec;
			m_vecvec[chan] = m_vec;
			DW_vecvec[chan] = DW_vec;
			r2_vecvec[chan] = r2_vec;
			thinned_chi_vecvec[chan] = chi_thinned;
			thinned_elev_vecvec[chan] = elev_thinned;
			fitted_elev_vecvec[chan] = elev_fitted;
			node_ref_thinned_vecvec[chan] = node_ref_thinned;
			these_segment_lengths_vecvec[chan] = these_segment_lengths;
			MLE_vec[chan] = this_MLE;
			n_segments_vec[chan] = this_n_segments;
			n_data_nodes_vec[chan] = n_data_nodes;
			AICc_vec[chan] = this_AICc;
			best_m_over_n[chan] = m_over_n;
		}

		AICc_mainstem_vec[movn] = this_AICc;
      	m_over_n_vec[movn] = m_over_n;
      	cout << "m/n: " << m_over_n << " and ms AICc is: " << this_AICc << endl;

    }


	// now go through all the tributaries
	int trib_N = 0;
	for (int chan = 1; chan<n_channels; chan++)
	{
		N = calculate_skip(target_nodes_mainstem,chan);
		cout << "LSDCN line 2519, channel " << chan << " and N: " << N << endl;
		if (N > trib_N)
		{
			trib_N = N;
		}
	}

	int N_ratio = (ms_N+1)/(trib_N+1);


	int trib_minimum_segment_length = minimum_segment_length/N_ratio;
	cout << "N_ratio: " << N_ratio << " and tminseglength: " << trib_minimum_segment_length << endl;

  	// loop through m_over_n values
	for(int movn = 0; movn< n_movern; movn++)
  	{
		m_over_n = double(movn)*d_movern+start_movern;
	  	cout << "m/n: " << m_over_n << endl;

      	// get the transformed channel profiles for this m_over_n
	  	calculate_chi(A_0, m_over_n);

   		vector<double> MLEs_thischan(n_channels);
		vector<int> n_segs_thischan(n_channels);
      	vector<int> n_datanodes_thischan(n_channels);

     	// now loop through channels
      	for (int chan = 1; chan<n_channels; chan++)
	  	{
			//cout << "LINE439 LSDChiNetwork channel: " << chan << endl;
		  	// get the channels for this m over n ratio
	  	  	find_most_likeley_segments(chan, trib_minimum_segment_length, sigma, trib_N,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
							      this_AIC, this_AICc );
			//cout << "LINE443 LSDChiNetwork channel: " << chan << endl;

			// check to see if the AICc value is the smallest
	  		// if so add the data to the best fit data elements
	  		if (this_AICc < AICc_vec[chan])
	    	{
	       		b_vecvec[chan] = b_vec;
	       		m_vecvec[chan] = m_vec;
	       		DW_vecvec[chan] = DW_vec;
	       		r2_vecvec[chan] = r2_vec;
	       		thinned_chi_vecvec[chan] = chi_thinned;
	       		thinned_elev_vecvec[chan] = elev_thinned;
	       		fitted_elev_vecvec[chan] = elev_fitted;
	       		node_ref_thinned_vecvec[chan] = node_ref_thinned;
	       		these_segment_lengths_vecvec[chan] = these_segment_lengths;
	       		MLE_vec[chan] = this_MLE;
	       		n_segments_vec[chan] = this_n_segments;
	       		n_data_nodes_vec[chan] = n_data_nodes;
	       		AICc_vec[chan] = this_AICc;
	       		best_m_over_n[chan] = m_over_n;
			}

	  		// add the data from this channel to the vectors that will be used to calcualte cumulative AICc
	  		MLEs_thischan[chan] = this_MLE;
	  		n_segs_thischan[chan] = this_n_segments;
	  		n_datanodes_thischan[chan] = n_data_nodes;
		}

      	//now calculate the cumulative AICc for this m over n
      	double thismn_AIC;
      	double thismn_AICc;

      	int n_total_segments = 0;
      	int n_total_nodes = 0;
      	double cumulative_MLE = 1;
      	double log_cum_MLE = 0;

		// get the cumulative maximum likelihood estimators
      	for (int chan = 1; chan<n_channels; chan++)
		{
			n_total_segments += n_segs_thischan[chan];
	  		n_total_nodes += n_datanodes_thischan[chan];
	  		cumulative_MLE = MLEs_thischan[chan]*cumulative_MLE;

			if(MLEs_thischan[chan] <= 0)
			{
				log_cum_MLE = log_cum_MLE-1000;
			}
			else
			{
				log_cum_MLE = log(MLEs_thischan[chan])+log_cum_MLE;
			}
	  		//cout << "LSDCN line 2591, chan: " << chan << " MLE: " << MLEs_thischan[chan] << endl;
		}


      	// these AIC and AICc values are cumulative for a given m_over_n
      	//thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
			                                                       // for each segment there are 2 parameters
      	thismn_AIC = 4*n_total_segments-2*log_cum_MLE;
      	thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);
      	AICc_combined_vec[movn] = thismn_AICc;
      	m_over_n_vec[movn] = m_over_n;

      	//cout << "m_over_n: " << m_over_n << " and combined AICc: " << thismn_AICc << endl;
      	//cout << "this cumulative MLE: " << cumulative_MLE << " n_segs: " << n_total_segments << " and n_nodes: " << n_total_nodes << endl;

    }

    cout << "and the mainstem m_over n values"<< endl;
    double min_cum_AICc = 9999;
    double bf_cum_movn = start_movern;
    for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		cout << "m over n: " << m_over_n_vec[mn] << " and AICc: " << AICc_mainstem_vec[mn] << endl;
		// if this is the minimum, store the m over n value
		if(AICc_mainstem_vec[mn] < min_cum_AICc)
	  	{
	   	 	min_cum_AICc = AICc_mainstem_vec[mn];
	    	bf_cum_movn = m_over_n_vec[mn];
	  	}
	}
	cout << "LSDChiNetwork Line 537 so the best fit m over n is: " << bf_cum_movn << endl;

    cout << "and the tributary m_over n values"<< endl;
    min_cum_AICc = 9999;
    double trib_bf_cum_movn = start_movern;
    for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		cout << "m over n: " << m_over_n_vec[mn] << " and AICc: " << AICc_combined_vec[mn] << endl;
		// if this is the minimum, store the m over n value
		if(AICc_combined_vec[mn] < min_cum_AICc)
	  	{
	   	 	min_cum_AICc = AICc_combined_vec[mn];
	    	trib_bf_cum_movn = m_over_n_vec[mn];
	  	}
	}
	cout << "LSDChiNetwork Line 2642 so the best fit tributary m over n is: " << trib_bf_cum_movn << endl;


    // now get the cumulative best fit channels
    calculate_chi(A_0, bf_cum_movn);
    // now loop through channels
    for (int chan = 0; chan<n_channels; chan++)
    {

		// get the channels for this m over n ratio
		find_most_likeley_segments(chan, minimum_segment_length, sigma, ms_N,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
                                            this_AIC, this_AICc );
	 	cum_b_vecvec[chan] = b_vec;
	 	cum_m_vecvec[chan] = m_vec;
        cum_DW_vecvec[chan] = DW_vec;
        cum_r2_vecvec[chan] = r2_vec;
        cum_thinned_chi_vecvec[chan] = chi_thinned;
        cum_thinned_elev_vecvec[chan] = elev_thinned;
        cum_fitted_elev_vecvec[chan] = elev_fitted;
        cum_node_ref_thinned_vecvec[chan] = node_ref_thinned;
        cum_these_segment_lengths_vecvec[chan] = these_segment_lengths;
	}

    // write a file
    ofstream best_fit_info;
	best_fit_info.open(fname.c_str());

  	best_fit_info << "N_channels: " << n_channels << endl;
  	best_fit_info << "m_over_n_for_channels ";
    for (int ch = 0; ch<n_channels; ch++)
    {
		best_fit_info << "  " << best_m_over_n[ch];
    }
  	best_fit_info << endl << "m_over_n_values ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << m_over_n_vec[mn];
	}
	best_fit_info << endl << "mainstem_AICc: ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << AICc_mainstem_vec[mn];
	}
	best_fit_info << endl << "cumulative_AICc: ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << AICc_combined_vec[mn];
	}
	best_fit_info << endl;
	cout << "The_best_fit_cumulative_m_over_n_is: " << bf_cum_movn << endl;
	for (int chan = 0; chan<n_channels; chan++)
	{
		vector<int> seglength = cum_these_segment_lengths_vecvec[chan];
		vector<double> m_val = cum_m_vecvec[chan];
		vector<double> b_val = cum_b_vecvec[chan];
		vector<double> DW_val = cum_DW_vecvec[chan];
		vector<double> r2_val = cum_r2_vecvec[chan];
		int n_segs_this_channel = seglength.size();

		best_fit_info << "Channel " << chan << " segment_length";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << seglength[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_gradient";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << m_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_intercept";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << b_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_DW_stat";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << DW_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_r2";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << r2_val[i];
		}
		best_fit_info << endl;
	}

	best_fit_info << endl << endl << "Cumulative_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = cum_thinned_chi_vecvec[chan];
		elev_thinned = cum_thinned_elev_vecvec[chan];
		elev_fitted = cum_fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }


	best_fit_info << endl << endl << "Individually_fit_channels:" << endl;
    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = thinned_chi_vecvec[chan];
		elev_thinned = thinned_elev_vecvec[chan];
		elev_fitted = fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }
     return bf_cum_movn;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function tests if the channels are long enough for a reasonable segment fit
// it writes to the is_tributary_long_enough vector
// if this equals 1, the channel is long enough. If it is zero the channel is not long enough
//
// SMM 15/03/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiNetwork::is_channel_long_enough_test(int minimum_segment_length,int N)
{

	int n_channels = int(chis.size());
	vector<int> chan_is_long(n_channels,0);
	for(int chan = 0; chan< n_channels; chan++)
	{
		vector<double> reverse_Chi = chis[chan];
		reverse(reverse_Chi.begin(), reverse_Chi.end());
		vector<double> reverse_Elevation = elevations[chan];
		reverse(reverse_Elevation.begin(), reverse_Elevation.end());
		vector<int> node_ref;

		LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);

		//int n_nodes = reverse_Chi.size();
		channel_MLE_finder.thin_data_skip(N, node_ref);
		int n_nodes_in_channel = channel_MLE_finder.get_n_nodes();

		//cout << "Testing channel lenghts, channel " << chan << " has " << n_nodes_in_channel << " nodes" << endl;

		if (n_nodes_in_channel >= 3*minimum_segment_length)
		{
			chan_is_long[chan] = 1;
		}
	}


	is_tributary_long_enough = chan_is_long;
}


#endif


