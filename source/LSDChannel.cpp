//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChannel
// Land Surface Dynamics Channel
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for retaining information of channels
//
// This is a derivative class of LSDIndexChannel.
//  LSDIndexChannel alone holds pointers to data in
//  LSDFlowInfo and LSDRaster, whereas LSDChannel
//  contains actual data about the channel such as
//  elevation and drainage area.
//
// These two objects are seperated to save on memory overhead
//  during runtime.
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

// Source code for the LSDChannel object.
// This obect is linked to a LSDIndexChannel object.
// The LSDIndexChannel object holds the row, column and node
// index for each channel. The LSDChannel object contains
// additional information such as elevation, drainage area
// and chi (the transformed coordiante for integral analysis of
// channel profiles

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------


#include <vector>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDChannel.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDMostLikelyPartitionsFinder.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChannel_CPP
#define LSDChannel_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
// creates an LSDChannel by copying from an IndexChannel
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(LSDIndexChannel& InChann)
{

	vector<double> empty_vec;
	Elevation = empty_vec;
	Chi = empty_vec;
	DrainageArea = empty_vec;

	StartJunction = InChann.get_StartJunction();
	EndJunction = InChann.get_EndJunction();
	StartNode = InChann.get_StartNode();
	EndNode = InChann.get_EndNode();

	NRows = InChann.get_NRows();
	NCols = InChann.get_NCols();
	XMinimum = InChann.get_XMinimum();
	YMinimum = InChann.get_YMinimum();
	DataResolution = InChann.get_DataResolution();
	NoDataValue = InChann.get_NoDataValue();

	RowSequence =  InChann.get_RowSequence();
	ColSequence =  InChann.get_ColSequence();
	NodeSequence =  InChann.get_NodeSequence();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calculates all the channel areas, elevations and chi parameters based on
// for a starting node index and ending node index
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJN, int EJN, double downslope_chi,
                             double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
{

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;
	double dx;
	double pixel_area = DataResolution*DataResolution;



	StartJunction = -1;
	EndJunction = -1;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;
		}
		else
		{
			curr_node = receive_node;
		}
	}

	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

	// get the number of nodes in the channel
	int n_nodes_in_channel = int(NodeSequence.size());

	// the bottom node is at chi of downslope_chi
	// initiate the chi vector
	vector<double> empty_vec;


	vector<double> chi_temp(n_nodes_in_channel,downslope_chi);


	vector<double> elev_temp(n_nodes_in_channel,double(NoDataValue));
	vector<double> area_temp(n_nodes_in_channel,double(NoDataValue));

	// get the first node
	double curr_area;

	curr_node = NodeSequence[n_nodes_in_channel-1];
	curr_row = RowI[n_nodes_in_channel-1];
	curr_col = ColI[n_nodes_in_channel-1];
	curr_area = double(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
	area_temp[n_nodes_in_channel-1] = curr_area;
	elev_temp[n_nodes_in_channel-1] = Elevation_Raster.get_data_element(curr_row,curr_col);

	// now loop up through the channel, adding chi values
	// note, the channel index are arranges with upstream element first, so you need to go through the channel
	// in reverse order
	for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
	{
	 	//cout << "ChIndex is: " << ChIndex << endl;
		curr_node = NodeSequence[ChIndex];
		FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);
		if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
		{
			dx = diag_length;
		}
		else
		{
			dx = DataResolution;
		}
		//cout << "dx is: " << dx << endl;

		curr_area = double(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
		area_temp[ChIndex] = curr_area;
		elev_temp[ChIndex] = Elevation_Raster.get_data_element(curr_row,curr_col);
		chi_temp[ChIndex] = dx*(pow( (A_0/curr_area ),
			                    m_over_n))
		                       + chi_temp[ChIndex+1];
		//cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
		//     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
	}

	Chi = chi_temp;
	Elevation = elev_temp;
	DrainageArea = area_temp;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calculates all the channel areas, elevations and chi parameters based on
// for a given LSDChannelIndex
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(double downslope_chi,
                             double m_over_n, double A_0, LSDIndexChannel& InChann, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
{


	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;
	double dx;
	double pixel_area = DataResolution*DataResolution;
	//cout << "data res: " << DataResolution << endl;



	StartJunction = InChann.get_StartJunction();
	EndJunction = InChann.get_EndJunction();
	StartNode = InChann.get_StartNode();
	EndNode = InChann.get_EndNode();

	RowSequence =  InChann.get_RowSequence();
	ColSequence =  InChann.get_ColSequence();
	NodeSequence =  InChann.get_NodeSequence();

	int curr_node = StartNode;

	// push back the data vectors with the starting node
	int curr_row, curr_col;

	// get the number of nodes in the channel
	int n_nodes_in_channel = int(NodeSequence.size());

	// the bottom node is at chi of downslope_chi
	// initiate the chi vector
	vector<double> empty_vec;
	vector<double> chi_temp(n_nodes_in_channel,downslope_chi);
	vector<double> elev_temp(n_nodes_in_channel,double(NoDataValue));
	vector<double> area_temp(n_nodes_in_channel,double(NoDataValue));

	// get the first node
	double curr_area;

	//cout << "downslope_chi: " << downslope_chi << endl;

	curr_node = NodeSequence[n_nodes_in_channel-1];
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);
	curr_area = double(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
	area_temp[n_nodes_in_channel-1] = curr_area;
	elev_temp[n_nodes_in_channel-1] = Elevation_Raster.get_data_element(curr_row,curr_col);

	// now loop up through the channel, adding chi values
	// note, the channel index are arranges with upstream element first, so you need to go through the channel
	// in reverse order
	for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
	{

		curr_node = NodeSequence[ChIndex];
		FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);

		//cout << "ChIndex is: " << ChIndex << " curr_node: " << curr_node << " row: "
                //                         << curr_row << " curr_col: " << curr_col << endl;

		if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
		{
			dx = diag_length;
		}
		else
		{
			dx = DataResolution;
		}
		//cout << "dx is: " << dx << endl;

		curr_area = double(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
		area_temp[ChIndex] = curr_area;
		elev_temp[ChIndex] = Elevation_Raster.get_data_element(curr_row,curr_col);
		chi_temp[ChIndex] = dx*(pow( (A_0/curr_area ),
			                    m_over_n))
		                       + chi_temp[ChIndex+1];
		//cout << "node " << curr_node << " and chi: " << chi_temp[ChIndex]
		//     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
	}

	Chi = chi_temp;
	Elevation = elev_temp;
	DrainageArea = area_temp;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// creates an index channel with just the node index of the starting and ending nodes
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJN, int EJN, LSDFlowInfo& FlowInfo)
{

	vector<double> empty_vec;
	Elevation = empty_vec;
	Chi = empty_vec;
	DrainageArea = empty_vec;

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	StartJunction = -1;
	EndJunction = -1;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;


		}
		else
		{
			curr_node = receive_node;
		}
	}

	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// third create routine
// creates an index channel with just the node index of the starting and ending nodes
// also includes junction information
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJ, int SJN, int EJ, int EJN, LSDFlowInfo& FlowInfo)
{

	vector<double> empty_vec;
	Elevation = empty_vec;
	Chi = empty_vec;
	DrainageArea = empty_vec;

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	StartJunction = SJ;
	EndJunction = EJ;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		//cout << "receive_node: " << receive_node << " and Endnode: " << EndNode << endl;

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;

		}
		else
		{
			curr_node = receive_node;
		}
	}
	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function uses a flow info object to calculate the chi values in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::calculate_chi(double downslope_chi, double m_over_n, double A_0, LSDFlowInfo& FlowInfo )
{
	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;
	double dx;
	double pixel_area = DataResolution*DataResolution;
	int curr_node;

	// get the number of nodes in the channel
	int n_nodes_in_channel = int(NodeSequence.size());

	// the bottom node is at chi of downslope_chi
	// initiate the chi vector
	vector<double> empty_vec;
	vector<double> chi_temp(n_nodes_in_channel,downslope_chi);

	// now loop up through the channel, adding chi values
	// note, the channel index are arranges with upstream element first, so you need to go through the channel
	// in reverse order
	for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
	{
	 	//cout << "ChIndex is: " << ChIndex << endl;
		curr_node = NodeSequence[ChIndex];
		if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
		{
			dx = diag_length;
		}
		else
		{
			dx = DataResolution;
		}
		//cout << "dx is: " << dx << endl;

		chi_temp[ChIndex] = dx*(pow( (A_0/ (double(
			                    FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area) ),
			                    m_over_n))
		                       + chi_temp[ChIndex+1];
		//cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
		//     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
	}

	Chi = chi_temp;
}


// this function gets the most likely channel segments
//
// SMM 2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::find_most_likeley_segments(int minimum_segment_length, double sigma, int target_nodes,
					     vector<double>& b_vec, vector<double>&  m_vec,
					     vector<double>& 	r2_vec,vector<double>&  DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_ref_thinned,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc )
{
	// first create a segment finder object
        //cout << "making MLEfinder object, " << endl;
	vector<int> empty_vec;
	node_ref_thinned = empty_vec;

	vector<double> reverse_Chi = Chi;
	reverse(reverse_Chi.begin(), reverse_Chi.end());
	vector<double> reverse_Elevation = Elevation;
	reverse(reverse_Elevation.begin(), reverse_Elevation.end());
	vector<int> this_node_sequence = NodeSequence;
	reverse(this_node_sequence.begin(), this_node_sequence.end());

	LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);
	//cout << "got MLE finder object" << endl;
	//cout << "rc size: " << reverse_Chi.size() << " and r_elev size: " << reverse_Elevation.size() << endl;
	//cout << "ns size: " << NodeSequence.size() << " and rns sz: " << this_node_sequence.size() << endl;
	// this needs to be thinned. Get the maximum chi value and then determine dchi
	int n_nodes = reverse_Chi.size();
	double max_chi = reverse_Chi[n_nodes-1];
	double min_chi = reverse_Chi[0];
	//cout << "LSDChannel::find_most_likeley_segments, max_chi: " << max_chi << " and min: " << min_chi << endl;
        //cout << "n_nodes is: " <<  channel_MLE_finder.get_n_nodes() << endl;

	double dchi = (max_chi-min_chi)/double(target_nodes);
	cout << "LSDChannel 533, dchi is: " << dchi << endl;


	// now thin the data, preserving the data (not interpoalting)
	vector<int> node_reference;
	channel_MLE_finder.thin_data_target_dx_preserve_data(dchi, node_reference);
	n_nodes = node_reference.size();
	//cout << "number of nodes in node reference: " << n_nodes << endl;
	for (int i = 0; i< n_nodes; i++)
	  {
	    //cout << " the node reference is: " << node_reference[i] << endl;
	    //cout << " node sequence: " << this_node_sequence[ node_reference[i]] << endl;
	    node_ref_thinned.push_back(this_node_sequence[ node_reference[i] ]);
	  }

	//cout << "thinned, n_nodes is: " <<  channel_MLE_finder.get_n_nodes() << endl;

	// now create a single sigma value vector
	vector<double> sigma_values;
	sigma_values.push_back(sigma);

	// compute the best fit AIC
	channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);


	channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
					     	r2_vec, DW_vec, fitted_elev,these_segment_lengths,
					       	this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

	thinned_chi = channel_MLE_finder.get_x_data();
	thinned_elev = channel_MLE_finder.get_y_data();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function finds the best fit m/n ratio
//
// SMM 2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::find_best_fit_m_over_n_with_segments(int n_movern, double d_movern,double start_movern,
						      double downslope_chi, double A_0, LSDFlowInfo& FlowInfo,
                                          int minimum_segment_length, double sigma, double target_nodes )
{
  // now get the details of the best fit
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

  // these are used for storing the best fits in the list of m_over_n
  vector<double> best_m_vec;
  vector<double> best_b_vec;
  vector<double> best_r2_vec;
  vector<double> best_DW_vec;
  vector<double> best_fitted_y;
  int best_n_data_nodes;
  int best_this_n_segments;
  double best_this_MLE, best_this_AIC, best_this_AICc;
  vector<int> best_these_segment_lengths;
  vector<double> best_chi_thinned;
  vector<double> best_elev_thinned;
  vector<double> best_elev_fitted;
  vector<int> best_node_ref_thinned;

  double min_AICc = 9999;
  double best_movern = start_movern;
  double m_over_n;
  for(int i = 0; i<n_movern; i++)
  {

    m_over_n = double(i)*d_movern+start_movern;

    // recalculate chi
    calculate_chi(downslope_chi, m_over_n,A_0, FlowInfo );

    find_most_likeley_segments(minimum_segment_length, sigma, target_nodes,
			 b_vec, m_vec,r2_vec,DW_vec,
			       chi_thinned, elev_thinned,elev_fitted, node_ref_thinned,
			 these_segment_lengths, this_MLE, this_n_segments, n_data_nodes,
			 this_AIC, this_AICc);

    if (this_AICc < min_AICc)
      {
	best_b_vec =b_vec;
	best_m_vec = m_vec;
	best_r2_vec = r2_vec;
	best_DW_vec = DW_vec;
	best_chi_thinned = chi_thinned;
	best_elev_thinned =  elev_thinned;
	best_elev_fitted = elev_fitted;
	best_node_ref_thinned = node_ref_thinned;
	best_these_segment_lengths = these_segment_lengths;
	best_this_MLE = this_MLE;
	best_this_n_segments =  this_n_segments;
	best_n_data_nodes =  n_data_nodes;
	best_this_AIC = this_AIC;
	best_this_AICc = this_AICc;

	min_AICc = this_AICc;
	best_movern = m_over_n;

	//cout << "best AICc: " << this_AICc << " and m_over_n: " << best_movern << endl;
      }
  }

  // now print the channel profile
  //cout << endl << endl << endl << "best fit m_over_n: " << best_movern << " with AICc: " << min_AICc << endl;
  //int n_nodes = chi_thinned.size();
  //for(int i = 0; i<n_nodes; i++)
  //  {
  //    cout << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
  //  }



}




#endif
