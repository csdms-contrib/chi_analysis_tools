//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo
// Land Surface Dynamics FlowInfo
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for organizing flow routing under the Fastscape algorithm
//  (see Braun and Willett, Geomorphology 2013, v180, p 170-179)
//
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

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo.cpp
// cpp file for the LSDFlowInfo object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.2		17/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// change log --this only starts Nov 28th 2012--
// 28/11/2012: Added write drainage area, FC, SMM
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDFlowInfo_CPP
#define LSDFlowInfo_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create function, this is empty, you need to include a filename
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::create()
{
	cout << "You need to initialize with a LSDRaster!" << endl;
	exit(EXIT_FAILURE);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create function, this creates from a pickled file
// fname is the name of the pickled flow info file
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::create(string fname)
{
	unpickle(fname);

}

//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the reciever of current_node (its node, row, and column)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_receiver_information(int current_node,
                                             int& receiver_node, int& receiver_row,
                                             int& receiver_col)
{
	int rn, rr, rc;
	rn = ReceiverVector[current_node];
	rr = RowIndex[rn];
	rc = ColIndex[rn];
	receiver_node = rn;
	receiver_row = rr;
	receiver_col = rc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the row and column of the current node
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_current_row_and_col(int current_node,int& curr_row,
                                             int& curr_col)
{
	int cr, cc;
	cr = RowIndex[current_node];
	cc = ColIndex[current_node];
	curr_row = cr;
	curr_col = cc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns the base level node with the greatest drainage area
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::retrieve_largest_base_level()
{
	int n_bl = BaseLevelNodeList.size();		// get the number of baselevel nodes
	int max_bl = 0;
	for (int i = 0; i<n_bl; i++)
	{
		if(NContributingNodes[ BaseLevelNodeList[i] ] > max_bl)
		{
			max_bl = NContributingNodes[ BaseLevelNodeList[i] ];
		}
	}
	return max_bl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function calcualtes the receiver nodes
// it returns the receiver vector r_i
// it also returns a flow direction array in this ordering:
//
// 7  0 1
// 6 -1 2
// 5  4 3
//
// note this is different from ArcMap flowdirection
//	int Arc_flowdir;			// flow direction in arcmap format
//								// 32  64  128
//								// 16  --  1
//								// 8    4  2
// one can convert nthese indices using the LSDIndexRaster object
// note in arc the row index increases down (to the south)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::create(vector<string> temp_BoundaryConditions,
										  LSDRaster& TopoRaster)
{

	// initialize several data members
	BoundaryConditions = temp_BoundaryConditions;

	NRows = TopoRaster.NRows;
	NCols = TopoRaster.NCols;
	XMinimum = TopoRaster.XMinimum;
	YMinimum = TopoRaster.YMinimum;
	NoDataValue = int(TopoRaster.NoDataValue);
	DataResolution = TopoRaster.DataResolution;

	// Declare matrices for calculating flow routing
	double one_ov_root2 = 0.707106781;
	double target_elev;				// a placeholder for the elevation of the potential receiver
	double slope;
	double max_slope;				// the maximum slope away from a node
	int max_slope_index;			// index into the maximum slope

	int row, col;						// index for the rows and column
	int receive_row,receive_col;
	string::iterator string_iterator;	// used to get characters from string

	// we need logic for all of the boundaries.
	// there are 3 kinds of edge boundaries:
	// no flux
	// base level
	// periodic
	// These are denoted in a vector of strings.
	// the vector has four elements
	// North boundary, East boundary, South bondary and West boundary
	// the strings can be any length, as long as the first letter corresponds to the
	// first letter of the boundary condition. It is not case sensitive.

	// go through the boundaries
	// A NOTE ON CARDINAL DIRECTIONS
	// If one looks at the raster data, the top row in the data corresponds to the NORTH boundary
	// This is the row that is first read into the code
	// so row 0 is the NORTH boundary
	// row NRows-1 is the SOUTH boundary
	// column 0 is the WEST boundary
	// column NCols-1 is the EAST boundary
	vector<double> slopes(8,NoDataValue);
	vector<int> row_kernal(8);
	vector<int> col_kernal(8);
	int ndv = NoDataValue;
	NDataNodes = 0; 			// the number of nodes in the raster that have data
	int one_if_a_baselevel_node;	// this is a switch used to tag baseleve nodes

	// the first thing you need to do is construct a topoglogy matrix
	// the donor, receiver, etc lists are as long as the number of nodes.
	// these are made of vectors that are exactly dimension n, which is the number of nodes with
	// data. Each of these nodes has several index vectors, that point the program to where the node is
	// we construct these index vectors first
	// we need to loop through all the data before we calcualte slopes because the
	// receiver node indices must be known before the slope calculations are run
	vector<int> empty_vec;
	RowIndex = empty_vec;
	ColIndex = empty_vec;
	BaseLevelNodeList = empty_vec;
	ReceiverVector = empty_vec;
	Array2D<int> ndv_raster(NRows,NCols,ndv);


	NodeIndex = ndv_raster.copy();
	FlowDirection = ndv_raster.copy();
	FlowLengthCode = ndv_raster.copy();


	// loop through the topo data finding places where there is actually data
	for (row = 0; row<NRows; row++)
	{
		for (col = 0; col<NCols; col++)
		{
			// only do calcualtions if there is data
			if(TopoRaster.RasterData[row][col] != NoDataValue)
			{
				RowIndex.push_back(row);
				ColIndex.push_back(col);
				NodeIndex[row][col] = NDataNodes;
				NDataNodes++;
			}
		}
	}
	// now the row and col index are populated by the row and col of the node in row i
	// and the node index has the indeces into the row and col vectors
	// next up, make d, delta, and D vectors
	vector<int> ndn_vec(NDataNodes,0);
	vector<int> ndn_nodata_vec(NDataNodes,ndv);
	vector<int> ndn_plusone_vec(NDataNodes+1,0);
	vector<int> w_vector(NDataNodes,0);

	NDonorsVector = ndn_vec;
	DonorStackVector = ndn_vec;
	DeltaVector = ndn_plusone_vec;

	SVector = ndn_nodata_vec;
	BLBasinVector = ndn_nodata_vec;

	// this vector starts out empty and then base level nodes are added to it
	for (row = 0; row<NRows; row++)
	{
		for (col = 0; col<NCols; col++)
		{

			// only do calcualtions if there is data
			if(TopoRaster.RasterData[row][col] != NoDataValue)
			{

				// calcualte 8 slopes
				// no slopes mean get NoDataValue entries
				// the algorithm loops through the neighbors to the cells, collecting
				// receiver indices. The order is
				// 7 0 1
				// 6 - 2
				// 5 4 3
				// where the above directions are cardinal directions
				// do slope 0
				row_kernal[0] = row-1;
				row_kernal[1] = row-1;
				row_kernal[2] = row;
				row_kernal[3] = row+1;
				row_kernal[4] = row+1;
				row_kernal[5] = row+1;
				row_kernal[6] = row;
				row_kernal[7] = row-1;

				col_kernal[0] = col;
				col_kernal[1] = col+1;
				col_kernal[2] = col+1;
				col_kernal[3] = col+1;
				col_kernal[4] = col;
				col_kernal[5] = col-1;
				col_kernal[6] = col-1;
				col_kernal[7] = col-1;


				// check for periodic boundary conditions
				if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
				{
					if( BoundaryConditions[2].find("P") != 0 && BoundaryConditions[2].find("p") != 0 )
					{
						cout << "WARNING!!! North boundary is periodic! Changing South boundary to periodic" << endl;
						BoundaryConditions[2] = "P";
					}
				}
				if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0 )
				{
					if( BoundaryConditions[3].find("P") != 0 && BoundaryConditions[3].find("p") != 0 )
					{
						cout << "WARNING!!! East boundary is periodic! Changing West boundary to periodic" << endl;
						BoundaryConditions[3] = "P";
					}
				}
				if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0 )
				{
					if( BoundaryConditions[0].find("P") != 0 && BoundaryConditions[0].find("p") != 0 )
					{
						cout << "WARNING!!! South boundary is periodic! Changing North boundary to periodic" << endl;
						BoundaryConditions[0] = "P";
					}
				}
				if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0 )
				{
					if( BoundaryConditions[1].find("P") != 0 && BoundaryConditions[1].find("p") != 0 )
					{
						cout << "WARNING!!! West boundary is periodic! Changing East boundary to periodic" << endl;
						BoundaryConditions[1] = "P";
					}
				}

				// reset baselevel switch for boundaries
				one_if_a_baselevel_node = 0;

				// NORTH BOUNDARY
				if (row == 0)
				{
					if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
					{
						one_if_a_baselevel_node = 1;
					}
					else
					{
						// if periodic, reflect across to south boundary
						if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
						{
							row_kernal[0] = NRows-1;
							row_kernal[1] = NRows-1;
							row_kernal[7] = NRows-1;
						}
						else
						{
							row_kernal[0] = ndv;
							row_kernal[1] = ndv;
							row_kernal[7] = ndv;
						}
					}
				}
				// EAST BOUNDAY
				if (col == NCols-1)
				{
					if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
					{
						one_if_a_baselevel_node = 1;
					}
					else
					{
						if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0)
						{
							col_kernal[1] = 0;
							col_kernal[2] = 0;
							col_kernal[3] = 0;
						}
						else
						{
							col_kernal[1] = ndv;
							col_kernal[2] = ndv;
							col_kernal[3] = ndv;
						}
					}
				}
				// SOUTH BOUNDARY
				if (row == NRows-1)
				{
					if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
					{
						one_if_a_baselevel_node = 1;
					}
					else
					{
						if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0)
						{
							row_kernal[3] = 0;
							row_kernal[4] = 0;
							row_kernal[5] = 0;
						}
						else
						{
							row_kernal[3] = ndv;
							row_kernal[4] = ndv;
							row_kernal[5] = ndv;
						}
					}
				}
				// WEST BOUNDARY
				if (col == 0)
				{
					if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
					{
						one_if_a_baselevel_node = 1;
					}
					else
					{
						if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0)
						{
							col_kernal[5] = NCols-1;
							col_kernal[6] = NCols-1;
							col_kernal[7] = NCols-1;
						}
						else
						{
							col_kernal[5] = ndv;
							col_kernal[6] = ndv;
							col_kernal[7] = ndv;
						}
					}
				}

				// now loop through the surrounding nodes, calcualting the slopes
				// slopes with NoData get NoData slopes
				// reminder of ordering:
				// 7 0 1
				// 6 - 2
				// 5 4 3
				// first logic for baselevel node
				if (one_if_a_baselevel_node == 1)
				{
					// get reciever index
					FlowDirection[row][col] = -1;
					ReceiverVector.push_back(NodeIndex[row][col]);
					FlowLengthCode[row][col] = 0;
				}
				// now the rest of the nodes
				else
				{
					FlowLengthCode[row][col] = 0;		// set flow length code to 0, this gets reset
														// if there is a maximum slope
					max_slope = 0;
					max_slope_index = -1;
					receive_row = row;
					receive_col = col;
					for (int slope_iter = 0; slope_iter<8; slope_iter++)
					{
						if (row_kernal[slope_iter] == ndv || col_kernal[slope_iter] == ndv)
						{
							slopes[slope_iter] = NoDataValue;
						}
						else
						{
							target_elev = TopoRaster.RasterData[ row_kernal[slope_iter] ][ col_kernal[slope_iter] ];
							if(target_elev == NoDataValue)
							{
								slopes[slope_iter] = NoDataValue;
							}
							else
							{
								if(slope_iter%2 == 0)
								{
									//cout << "LINE 988, cardinal direction, slope iter = " << slope_iter << endl;
									slope = TopoRaster.RasterData[row][col]-target_elev;
								}
								else
								{
									slope = one_ov_root2*(TopoRaster.RasterData[row][col]-target_elev);
								}

								if (slope > max_slope)
								{
									max_slope_index = slope_iter;
									receive_row = row_kernal[slope_iter];
									receive_col = col_kernal[slope_iter];
									max_slope = slope;
									if(slope_iter%2 == 0)
									{
										FlowLengthCode[row][col] = 1;
									}
									else
									{
										FlowLengthCode[row][col] = 2;
									}
								}
							}
						}
					}
					// get reciever index
					FlowDirection[row][col] = max_slope_index;
					ReceiverVector.push_back(NodeIndex[receive_row][receive_col]);
				}		// end if baselevel boundary  conditional

				// if the node is a base level node, add it to the base level node list
				if (FlowLengthCode[row][col] == 0)
				{
					BaseLevelNodeList.push_back(NodeIndex[row][col]);
				}
			}			// end if there is data conditional
		}				// end col loop
	}					// end row loop

	//cout << "LINE 1015, NDataNodes: " << NDataNodes
	//     << " and size reciever: " << ReceiverVector.size() << endl;

	//ofstream receiver_out;
	//receiver_out.open("receiver_out.txt");
	//cout << "LINE 414 LSDFlow info, writing receiver nodes to file" << endl;
	//for (int i =0; i<int(ReceiverVector.size()); i++)
	//{
	//	receiver_out << ReceiverVector[i] << endl;
	//}
	//receiver_out.close();

	// first create the number of donors vector
	// from braun and willett eq. 5
	for(int i = 0; i<NDataNodes; i++)
	{
		NDonorsVector[ ReceiverVector[i] ]++;
	}

	// now create the delta vector
	// this starts on the last element and works its way backwards
	// from Braun and Willett eq 7 and 8
	DeltaVector[NDataNodes] = NDataNodes;
	for(int i = NDataNodes; i>0; i--)
	{
		DeltaVector[i-1] = DeltaVector[i] -  NDonorsVector[i-1];
	}

	// now the DonorStack and the r vectors. These come from Braun and Willett
	// equation 9.
	// Note that in the manscript I have there is a typo in eqaution 9
	// (Jean Braun's code is correct)
	// it should be w_{r_i} = w_{r_i}+1
	int r_index;
	int w_index;
	int delta_index;
	for (int i = 0; i<NDataNodes; i++)
	{
		r_index = ReceiverVector[i];
		delta_index = DeltaVector[ r_index ];
		w_index = w_vector[ r_index ];
		DonorStackVector[  delta_index+w_index ] = i;
		w_vector[r_index] += 1;
		//cout << "i: " << i << " r_i: " << r_index << " delta_i: " << delta_index << " w_index: " << w_index << endl;
	}

	// now go through the base level node list, building the drainage tree for each of these nodes as one goes along
	int n_base_level_nodes;
	n_base_level_nodes = BaseLevelNodeList.size();

	int k;
	int j_index;
	int begin_delta_index, end_delta_index;
	int l_index;

	j_index = 0;
	for (int i = 0; i<n_base_level_nodes; i++)
	{
		k = BaseLevelNodeList[i];			// set k to the base level node

		// This doesn't seem to be in Braun and Willet but to get the ordering correct you
		// need to make sure that the base level node appears first in the donorstack
		// of nodes contributing to the baselevel node.
		// For example, if base level node is 4, with 4 donors
		// and the donor stack has 3 4 8 9
		// the code has to put the 4 first.
		if (DonorStackVector[ DeltaVector[k] ] != k)
		{
			int this_index = DonorStackVector[ DeltaVector[k] ];
			int bs_node = k;

			for(int ds_node = 1; ds_node < NDonorsVector[k]; ds_node++)
			{
				if( DonorStackVector[ DeltaVector[k] + ds_node ] == bs_node )
				{
					DonorStackVector[ DeltaVector[k] ] = k;
					DonorStackVector[ DeltaVector[k] + ds_node ] = this_index;
				}
			}
		}

		// now run recursive algorithm
		begin_delta_index = DeltaVector[k];
		end_delta_index = DeltaVector[k+1];

		//cout << "base_level_node is: " << k << " begin_index: " << begin_delta_index << " end: " << end_delta_index << endl;

		for (int delta_index = begin_delta_index; delta_index<end_delta_index; delta_index++)
		{
			l_index = DonorStackVector[delta_index];
			add_to_stack(l_index, j_index, k);
		}
	}

	// now calcualte the indices
	calculate_upslope_reference_indices();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// recursive add_to_stack routine, from Braun and Willett eq. 12 and 13
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::add_to_stack(int lm_index, int& j_index, int bl_node)
{
	//cout << "j_index: " << j_index << " and s_vec: " << lm_index << endl;

	SVector[j_index] = lm_index;
	BLBasinVector[j_index] = bl_node;
	j_index++;


	int begin_m,end_m;
	int l_index;
	// if donating to itself, need escape hatch
	if ( lm_index == bl_node)
	{
		begin_m = 0;
		end_m = 0;
	}
	else
	{
		begin_m = DeltaVector[lm_index];
		end_m =  DeltaVector[ lm_index+1];
	}
	//cout << "lm_index: " << lm_index << " begin_m: " << begin_m << " end m: " << end_m << endl;
	for( int m_index = begin_m; m_index<end_m; m_index++)
	{
		//cout << "recursion, begin_m: " << begin_m << " and end_m: " << end_m << endl;
		l_index = DonorStackVector[m_index];
		add_to_stack(l_index, j_index, bl_node);
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function pickles the data from the flowInfo object into a binary format
// which can be read by the unpickle function later
// the filename DOES NOT include and extension: this is added by the
// function
//
// WARNING: These files are HUGE and testing indicates they don't save much time
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::pickle(string filename)
{
	string ext = ".FIpickle";
	string hdr_ext = ".FIpickle.hdr";

	string hdr_fname = filename+hdr_ext;
	string data_fname = filename+ext;

	ofstream header_out;
	header_out.open(hdr_fname.c_str());

	int contributing_nodes = int(NContributingNodes.size());
	int BLNodes = int(BaseLevelNodeList.size());

	// print the header file
	header_out <<  "ncols         		" << NCols
			   << "\nnrows         		" << NRows
			   << "\nxllcorner     		" << setprecision(14) << XMinimum
			   << "\nyllcorner     		" << setprecision(14) << YMinimum
			   << "\ncellsize      		" << DataResolution
			   << "\nNODATA_value  		" << NoDataValue
			   << "\nNDataNodes    		" << NDataNodes
			   << "\nNBaseLevelNodes    " << BLNodes
			   << "\nNContributingNodes " << contributing_nodes
			   << "\nBoundaryConditions ";
	for(int i = 0; i<4; i++)
	{
		header_out << " " << BoundaryConditions[i];
	}
	header_out << endl;
	header_out.close();


	cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
	cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
	cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
	cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
	cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
	cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;


	// now do the main data
	ofstream data_ofs(data_fname.c_str(), ios::out | ios::binary);
	int temp;
	for (int i=0; i<NRows; ++i)
	{
		for (int j=0; j<NCols; ++j)
		{
			temp = NodeIndex[i][j];
			data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
		}
	}
	for (int i=0; i<NRows; ++i)
	{
		for (int j=0; j<NCols; ++j)
		{
			temp = FlowDirection[i][j];
			data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
		}
	}
	for (int i=0; i<NRows; ++i)
	{
		for (int j=0; j<NCols; ++j)
		{
			temp = FlowLengthCode[i][j];
			data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
		}
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = RowIndex[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = ColIndex[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<BLNodes; i++)
	{
		temp = BaseLevelNodeList[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = NDonorsVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = ReceiverVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes+1; i++)
	{
		temp = DeltaVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = DonorStackVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = SVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = BLBasinVector[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<NDataNodes; i++)
	{
		temp = SVectorIndex[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
	for (int i = 0; i<contributing_nodes; i++)
	{
		temp = NContributingNodes[i];
		data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}

	data_ofs.close();


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this unpickles a pickled flow info object. It is folded into a create function
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::unpickle(string filename)
{

	string ext = ".FIpickle";
	string hdr_ext = ".FIpickle.hdr";

	string hdr_fname = filename+hdr_ext;
	string data_fname = filename+ext;

	ifstream header_in;
	header_in.open(hdr_fname.c_str());

	string temp_str;
	int contributing_nodes;
	vector<string> bc(4);
	int BLNodes;

	header_in >> temp_str >> NCols >> temp_str >> NRows >> temp_str >> XMinimum
	          >> temp_str >> YMinimum >> temp_str >> DataResolution
			  >> temp_str >> NoDataValue >> temp_str >> NDataNodes
			  >> temp_str >> BLNodes
			  >> temp_str >> contributing_nodes
			  >> temp_str >> bc[0] >> bc[1] >> bc[2] >> bc[3];
	header_in.close();
	BoundaryConditions = bc;


	// now read the data, using the binary stream option
	ifstream ifs_data(data_fname.c_str(), ios::in | ios::binary);
	if( ifs_data.fail() )
	{
		cout << "\nFATAL ERROR: the data file \"" << data_fname
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		// initialze the arrays
		Array2D<int> data_array(NRows,NCols,NoDataValue);
		NodeIndex = data_array.copy();
		FlowDirection = data_array.copy();
		FlowLengthCode = data_array.copy();

		vector<int> data_vector(NDataNodes,NoDataValue);
		vector<int> BLvector(BLNodes,NoDataValue);
		vector<int> deltaV(NDataNodes+1,NoDataValue);
		vector<int> CNvec(contributing_nodes,NoDataValue);

		int temp;
		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
				NodeIndex[i][j] =temp;
			}
		}
		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
				FlowDirection[i][j] =temp;
			}
		}
		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
				FlowLengthCode[i][j] =temp;
			}
		}
		RowIndex = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			RowIndex[i] =temp;

		}
		ColIndex = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			ColIndex[i] =temp;

		}
		BaseLevelNodeList = BLvector;
		for (int i=0; i<BLNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			BaseLevelNodeList[i] =temp;

		}
		NDonorsVector = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			NDonorsVector[i] =temp;

		}
		ReceiverVector = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			ReceiverVector[i] =temp;

		}
		DeltaVector = deltaV;
		for (int i=0; i<NDataNodes+1; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			DeltaVector[i] =temp;

		}
		DonorStackVector = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			DonorStackVector[i] =temp;

		}
		SVector = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			SVector[i] =temp;

		}
		BLBasinVector = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			BLBasinVector[i] =temp;

		}
		SVectorIndex = data_vector;
		for (int i=0; i<NDataNodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			SVectorIndex[i] =temp;

		}
		NContributingNodes = CNvec;
		for (int i=0; i<contributing_nodes; ++i)
		{
			ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
			NContributingNodes[i] =temp;

		}


	}
	ifs_data.close();

	cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
	cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
	cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
	cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
	cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
	cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the flownet information
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_flow_info_vectors(string filename)
{
	string string_filename;
	string dot = ".";
	string extension = "txt";
	string_filename = filename+dot+extension;
	cout << "The filename is " << string_filename << endl;

	// print out all the donor, reciever and stack info
	ofstream donor_info_out;
	donor_info_out.open(string_filename.c_str());
	for(int i = 0; i<NDataNodes; i++)
	{
		donor_info_out << i << " ";
	}
	donor_info_out << endl;
	for(int i = 0; i<NDataNodes; i++)
	{
		donor_info_out << ReceiverVector[i] << " ";
	}
	donor_info_out << endl;
	for(int i = 0; i<NDataNodes; i++)
	{
		donor_info_out << NDonorsVector[i] << " ";
	}
	donor_info_out << endl;
	for(int i = 0; i<NDataNodes+1; i++)
	{
		donor_info_out << DeltaVector[i] << " ";
	}
	donor_info_out << endl;
	for(int i = 0; i<NDataNodes; i++)
	{
		donor_info_out << DonorStackVector[i] << " ";
	}
	donor_info_out << endl;
	for(int i = 0; i<NDataNodes; i++)
	{
		donor_info_out << SVector[i] << " ";
	}
	donor_info_out << endl;

	if( int(SVectorIndex.size()) == NDataNodes)
	{
		for(int i = 0; i<NDataNodes; i++)
		{
			donor_info_out << SVectorIndex[i] << " ";
		}
		donor_info_out << endl;
		for(int i = 0; i<NDataNodes; i++)
		{
			donor_info_out << NContributingNodes[i] << " ";
		}
		donor_info_out << endl;
	}

	donor_info_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// these functions write the index arrays to index rasters
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndex_to_LSDIndexRaster()
{
	LSDIndexRaster temp_nodeindex(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,NodeIndex);
	return temp_nodeindex;
}

LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster()
{
	LSDIndexRaster temp_flowdir(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirection);
	return temp_flowdir;
}

LSDIndexRaster LSDFlowInfo::write_FlowLengthCode_to_LSDIndexRaster()
{
	LSDIndexRaster temp_flc(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowLengthCode);
	return temp_flc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// This function used the S vector index and NContributing nodes to calculate the area
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NContributingNodes_to_LSDIndexRaster()
{
	Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
	int row,col;

	// loop through the node vector, adding pixels to receiver nodes
	for(int node = 0; node<NDataNodes; node++)
	{
		row = RowIndex[node];
		col = ColIndex[node];
		contributing_pixels[row][col] = NContributingNodes[node];
	}

	LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels);
	return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// writes a index raster in arc format
// LSD format:
// 7  0 1
// 6 -1 2
// 5  4 3
//	int Arc_flowdir;			// flow direction in arcmap format
//								// 32  64  128
//								// 16   0  1
//								// 8    4  2
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster_Arcformat()
{

	Array2D<int> FlowDirectionArc(NRows,NCols,NoDataValue);
	for(int row = 0; row<NRows; row++)
	{
		for (int col = 0; col<NCols; col++)
		{
			if ( FlowDirection[row][col] == -1)
			{
				FlowDirectionArc[row][col] = 0;
			}
			else if ( FlowDirection[row][col] == 0)
			{
				FlowDirectionArc[row][col] = 64;
			}
			else if ( FlowDirection[row][col] == 1)
			{
				FlowDirectionArc[row][col] = 128;
			}
			else if ( FlowDirection[row][col] == 2)
			{
				FlowDirectionArc[row][col] = 1;
			}
			else if ( FlowDirection[row][col] == 3)
			{
				FlowDirectionArc[row][col] = 2;
			}
			else if ( FlowDirection[row][col] == 4)
			{
				FlowDirectionArc[row][col] = 4;
			}
			else if ( FlowDirection[row][col] == 5)
			{
				FlowDirectionArc[row][col] = 8;
			}
			else if ( FlowDirection[row][col] == 6)
			{
				FlowDirectionArc[row][col] = 16;
			}
			else if ( FlowDirection[row][col] == 7)
			{
				FlowDirectionArc[row][col] = 32;
			}
		}
	}

	LSDIndexRaster temp_fd(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirectionArc);
	return temp_fd;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes the drainage area (number of contributing nodes
// * DataResolution^2) to an LSDRaster object
// Added by FC 15/11/12
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::write_DrainageArea_to_LSDRaster()
{
  // initialise the 2D array
  int n_i;								// node index
  double ndv = double(NoDataValue);
 	Array2D<double> DrainageArea_local(NRows,NCols,ndv);

  //get the contributing nodes
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      n_i = NodeIndex[row][col];
      DrainageArea_local[row][col] = double(NContributingNodes[n_i])*DataResolution*DataResolution;
    }
  }
  // create the LSDRaster object
  LSDRaster DrainageArea(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,DrainageArea_local);
	return DrainageArea;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
// In this function a pixel that has no donors has a contributing pixel value of 0
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::calculate_n_pixels_contributing_from_upslope()
{
	Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
	int row,col;
	int receive_row, receive_col;
	int receiver_node;

	// loop through the s vector, adding pixels to receiver nodes
	for(int node = NDataNodes-1; node>=0; node--)
	{

		row = RowIndex[SVector[node]];
		col = ColIndex[SVector[node]];
		// if the pixel exists and has no contributing pixels,
		// change from nodata to zero

		if(contributing_pixels[row][col] == NoDataValue)
		{
			contributing_pixels[row][col] = 0;
		}

		receiver_node = ReceiverVector[ SVector[node] ] ;
		receive_row = RowIndex[ receiver_node ];
		receive_col = ColIndex[ receiver_node ];

		cout << "node " << node << " pixel: " << SVector[node] << " receiver: " << receiver_node << endl;
		cout << "contributing: " << contributing_pixels[row][col] << endl;

		if ( receiver_node  == SVector[node])
		{
			// do nothing
		}
		else if ( contributing_pixels[receive_row][receive_col] == NoDataValue)
		{
			contributing_pixels[receive_row][receive_col] =
				contributing_pixels[row][col]+1;
		}
		else
		{
			contributing_pixels[receive_row][receive_col] +=
				contributing_pixels[row][col]+1;
		}

		cout << "recieving: " << contributing_pixels[receive_row][receive_col] << endl;
	}

	LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels);
	return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// In this function a pixel that has no donors contributes its own flow
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::calculate_upslope_reference_indices()
{
	vector<int> vectorized_area(NDataNodes,1);
	SVectorIndex = vectorized_area;

	int receiver_node;
	int donor_node;

	// loop through the s vector, adding pixels to receiver nodes
	for(int node = NDataNodes-1; node>=0; node--)
	{
		donor_node = SVector[node];
		receiver_node = ReceiverVector[ donor_node ];

		// every node is visited once and only once so we can map the
		// unique positions of the nodes to the SVector
		SVectorIndex[donor_node] = node;

		// add the upslope area (note no action is taken
		// for base level nodes since they donate to themselves and
		// we must avoid double counting
		if (donor_node != receiver_node)
		{
			vectorized_area[ receiver_node ] +=  vectorized_area[ donor_node ];
		}
	}

	NContributingNodes = vectorized_area;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the node numbers upslope
// of of the node with number node_number_outlet
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_upslope_nodes(int node_number_outlet)
{
	vector<int> us_nodes;

	if(node_number_outlet < 0 || node_number_outlet > NDataNodes-1)
	{
		cout << "the junction number does not exist" << endl;
		exit(EXIT_FAILURE);
	}

	int start_SVector_node = SVectorIndex[node_number_outlet];
	int end_SVector_node = start_SVector_node+NContributingNodes[node_number_outlet];

	for(int node = start_SVector_node; node < end_SVector_node; node++)
	{
		us_nodes.push_back(SVector[node]);
	}

	return us_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the chi function for all the nodes upslope a given node
// it takes a node list
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> LSDFlowInfo::get_upslope_chi(int starting_node, double m_over_n, double A_0)
{
	vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);
	vector<double> chi_vec = get_upslope_chi(upslope_pixel_list, m_over_n, A_0);
	return chi_vec;
}

vector<double> LSDFlowInfo::get_upslope_chi(vector<int>& upslope_pixel_list, double m_over_n, double A_0)
{

	int receiver_node;
	int IndexOfReceiverInUplsopePList;
	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;
	double dx;
	double pixel_area = DataResolution*DataResolution;
	int node,row,col;
	// get the number of nodes upslope
	int n_nodes_upslope = upslope_pixel_list.size();
	vector<double> chi_vec(n_nodes_upslope,0.0);

	if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
	{
		cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
		exit(EXIT_FAILURE);
	}

	int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];

	for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
	{
		node = upslope_pixel_list[n_index];
		receiver_node = ReceiverVector[ node ];
		IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
		row = RowIndex[node];
		col = ColIndex[node];

		if (FlowLengthCode[row][col] == 2)
		{
			dx = diag_length;
		}
		else
		{
			dx = DataResolution;
		}


		chi_vec[n_index] = dx*(pow( (A_0/ (double(NContributingNodes[node])*pixel_area) ),m_over_n))
		                       + chi_vec[IndexOfReceiverInUplsopePList];

		cout << "node: " << upslope_pixel_list[n_index] << " receiver: " << receiver_node
		     << " SIndexReciever: " << IndexOfReceiverInUplsopePList << " and checked: " << upslope_pixel_list[IndexOfReceiverInUplsopePList]
		     << " and chi: " << chi_vec[n_index] << endl;


	}

	return chi_vec;
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// distance from outlet function
// this is overloaded.
// if it isn't provided any argument, it calculates the distance from outlet
// of all the base level nodes
// if it is given a node index number or a row and column, then
// the distance from outlet includes all the distances upstream of that node
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::distance_from_outlet()
{
	// initialize the array2d that will become the LSDRaster
	double ndv = double(NoDataValue);
	Array2D<double> flow_distance(NRows,NCols,ndv);

	// initialize the 1/root(2)
	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;

	int row,col,bl_row,bl_col,receive_row,receive_col;

	//cout << "FlowLengthCode: " << FlowLengthCode << endl;

	int start_node = 0;
	int end_node;
	int nodes_in_bl_tree;
	int baselevel_node;
	// loop through the base level node list
	int n_base_level_nodes = BaseLevelNodeList.size();
	for(int bl = 0; bl<n_base_level_nodes; bl++)
	{

		baselevel_node = BaseLevelNodeList[bl];

		bl_row = RowIndex[baselevel_node];
		bl_col = ColIndex[baselevel_node];
		// get the number of nodes upslope and including this node
		nodes_in_bl_tree = NContributingNodes[baselevel_node];
		//cout << "LINE 938, FlowInfo, base level: " << bl << " with " << nodes_in_bl_tree << " nodes upstream" << endl;

		end_node = start_node+nodes_in_bl_tree;

		// set the distance of the outlet to zero
		flow_distance[bl_row][bl_col] = 0;

		// now loop through stack
		for(int s_node = start_node; s_node < end_node; s_node++)
		{
			//cout << "Line 953 flow info, s_node is: " << s_node << endl;

			//cout << SVector.size() << " " << ReceiverVector.size() << " " << RowIndex.size() << " " << ColIndex.size() << endl;
			row = RowIndex[ SVector[ s_node]  ];
			col = ColIndex[ SVector[ s_node]  ];
			//cout << "got rows and columns " << row << " " << col << endl;
			receive_row = RowIndex[ ReceiverVector[SVector[s_node] ]];
			receive_col = ColIndex[ ReceiverVector[SVector[s_node] ]];
			//cout <<  "get receive " << receive_row << " " << receive_col << endl;

			if ( FlowLengthCode[row][col] == 1)
			{
				flow_distance[row][col] = flow_distance[receive_row][receive_col]+DataResolution;
			}
			else if ( FlowLengthCode[row][col] == 2 )
			{
				flow_distance[row][col] = flow_distance[receive_row][receive_col]
											+ diag_length;
			}
			//cout << "Flow distance: " << flow_distance << endl;
		}
		start_node = end_node;
	}
	//cout << "LINE 971 FlowInfo Flow distance complete, flow_distance is: " << endl;
	//cout << flow_distance << endl;
	LSDRaster FlowLength(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,flow_distance);
	return FlowLength;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses the flow accumulation pixels
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold)
{
	vector<int> sources;
	int row,col;
	//int n_donors;
	int donor_row,donor_col;
	int thresh_switch;
	int donor_node;

	// drop down through the stack
	// if the node is greater than or equal to the threshold
	// check all the donors
	// if there are no donors, then the node is a source
	// if none of the donors are greater than the threshold, then it also is a source
	for (int node = 0; node<NDataNodes; node++)
	{
		row = RowIndex[node];
		col = ColIndex[node];

		// see if node is greater than threshold
		if(FlowPixels.get_data_element(row,col)>=threshold)
		{
			//cout << "node " << node << " is a potential source, it has a value of "
			//     << FlowPixels.get_data_element(row,col)
			//     << "and it has " << NDonorsVector[node] <<" donors " << endl;

			// if it doesn't have donors, it is a source
			if(NDonorsVector[node] == 0)
			{
				sources.push_back(node);
			}
			else
			{
				thresh_switch = 1;
				// figure out where the donor nodes are, and if
				// the donor node is greater than the threshold
				for(int dnode = 0; dnode<NDonorsVector[node]; dnode++)
				{
					donor_node = DonorStackVector[ DeltaVector[node]+dnode];
					donor_row = RowIndex[ donor_node ];
					donor_col = ColIndex[ donor_node ];

					// we don't double count base level nodes, which donate to themselves
					if (donor_node != node)
					{
						// if the donor node is greater than the threshold,
						// then this node is not a threhold
						if(FlowPixels.get_data_element(donor_row,donor_col)>=threshold)
						{
							thresh_switch = 0;
						}
					}

					//cout << "thresh_switch is: " << thresh_switch << endl;
				}
				// if all of the donors are below the threhold, this is a source
				if (thresh_switch == 1)
				{
					sources.push_back(node);
				}
			}
		}
	}
	return sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




#endif
