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


/** @file LSDFlowInfo.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0
@brief Object to perform flow routing.
@details This is a data object which generates and then
stores information about flow routing.
It is the object that is used to generate contributing area, etc.

@date 29/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster_lw.hpp"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDFlowInfo_H
#define LSDFlowInfo_H

/// @brief Object to perform flow routing.
class LSDFlowInfo
{
	public:
  /// @brief The create function. This is default and throws an error.
  /// @author SMM
  /// @date 01/016/12
	LSDFlowInfo()										{ create(); }
	/// @brief Creates a FlowInfo object from a binary flowinfo data.
	/// @param fname String of the binary flowinfo data file to be read.
	/// @author SMM
    /// @date 01/016/12
	LSDFlowInfo(string fname)							{ create(fname); }
	/// @brief Creates a FlowInfo object from topography.
	/// @param BoundaryConditions Vector<string> of the boundary conditions at each edge of the
    /// DEM file. Boundary conditions can start with 'P' or 'p' for periodic,
    /// 'B' or 'b' for base level, or anything else for no flux.
    /// the vector shold have 4 elements, 0 is north, 1 is east, 2 is south and 3 is west
	/// @param TopoRaster LSDRaster object containing the topographic data.
	/// @author SMM
    /// @date 01/016/12
	LSDFlowInfo(vector<string> BoundaryConditions, LSDRaster& TopoRaster)
									{ create(BoundaryConditions, TopoRaster); }

	/// @brief Copy of the LSDChannelNetwork description here when written.
	friend class LSDChannelNetwork;

	// some functions for retrieving information out of the data vectors

  ///@brief Gives the reciever information for a given node.
  ///@param current_node Integer
  ///@param reveiver_node Empty integer to be assigned the index of the reciever
  ///node.
  ///@param receiver_row Empty integer to be assigned the row index of the
  ///reciever node.
  ///@param receiver_col Empty integer to be assigned the column index of the
  ///reciever node.
	/// @author SMM
    /// @date 01/016/12
	void retrieve_receiver_information(int current_node,int& reveiver_node, int& receiver_row,
                                             int& receiver_col);

  ///@brief Get the row and column indices of a given node.
  ///@param current_node Integer index of a given node.
  ///@param curr_row Empty integer to be assigned the row index of the given
  ///node.
  ///@param curr_col Empty integer to be assigned the column index of the given
  ///node.
  	/// @author SMM
    /// @date 01/016/12
	void retrieve_current_row_and_col(int current_node,int& curr_row,
                                             int& curr_col);

  ///@brief Get the number of pixels flowing into a node.
  ///@param node Integer of node index value.
  ///@return Integer of the number of contributing pixels.
  	/// @author SMM
    /// @date 01/016/12
	int retrieve_contributing_pixels_of_node(int node)
										{ return NContributingNodes[node]; }

  ///@brief Get the FlowLengthCode of a given node.
  ///@param node Integer of node index value.
  ///@return Integer of the FlowLengthCode.
  	/// @author SMM
    /// @date 01/016/12
  int retrieve_flow_length_code_of_node(int node)
										{ return FlowLengthCode[ RowIndex[node] ][ ColIndex[node] ]; }

	// get functions

	/// @return Number of rows as an integer.
	int get_NRows() const				{ return NRows; }
	/// @return Number of columns as an integer.
  int get_NCols() const				{ return NCols; }
  /// @return Minimum X coordinate as an integer.
	double get_XMinimum() const			{ return XMinimum; }
	/// @return Minimum Y coordinate as an integer.
	double get_YMinimum() const			{ return YMinimum; }
	/// @return Data resolution as an integer.
	double get_DataResolution() const	{ return DataResolution; }
	/// @return No Data Value as an integer.
	int get_NoDataValue() const			{ return NoDataValue; }
	/// @return Number of nodes with data as an integer.
	int get_NDataNodes () const			{ return NDataNodes; }
	/// @return Vector of all base level nodes.
  vector<int> get_BaseLevelNodeList ()
										{ return BaseLevelNodeList; }
  /// @return FlowDirection values as a 2D Array.
	Array2D<int> get_FlowDirection() const { return FlowDirection; }

  ///@brief Recursive add_to_stack routine, from Braun and Willett (2012)
  ///equations 12 and 13.
  ///@param lm_index Integer
  ///@param j_index Integer
  ///@param bl_node Integer
	void add_to_stack(int lm_index, int& j_index, int bl_node);

	// some functions that print out indices to rasters
	///@brief Write NodeIndex to an LSDIndexRaster.
  ///@return LSDIndexRaster of node index data.
  	/// @author SMM
    /// @date 01/016/12
	LSDIndexRaster write_NodeIndex_to_LSDIndexRaster();
	///@brief Write FlowDirection to an LSDIndexRaster.
  ///@return LSDIndexRaster of flow directions.
  	/// @author SMM
    /// @date 01/016/12
  LSDIndexRaster write_FlowDirection_to_LSDIndexRaster();
	///@brief Write FlowLengthCode to an LSDIndexRaster.
  ///@return LSDIndexRaster of flow lengths.
  	/// @author SMM
    /// @date 01/016/12
	LSDIndexRaster write_FlowLengthCode_to_LSDIndexRaster();
	///@brief Write NContributingNodes to an LSDIndexRaster.
  ///@return LSDIndexRaster of number of contributing nodes for each cell.
  	/// @author SMM
    /// @date 01/016/12
  LSDIndexRaster write_NContributingNodes_to_LSDIndexRaster();
	///@brief Writes flow directions to an LSDIndexRaster.
  ///Flow direction in arcmap format is: \n\n
  /// 32  64  128 \n
  /// 16   0  1    \n
  /// 8    4  2     \n
  ///
  ///@return LSDIndexRaster of flow directions in arcgis format.
  	/// @author SMM
    /// @date 01/016/12
	LSDIndexRaster write_FlowDirection_to_LSDIndexRaster_Arcformat();
	///@brief
  ///@return
  ///@author Fiona Clubb
  ///@date 15/11/12
  LSDRaster write_DrainageArea_to_LSDRaster();

	///@brief Prints the flow information to file.
	///@param filename String of the output file to be written.
		/// @author SMM
    /// @date 01/016/12
	void print_flow_info_vectors(string filename);

  ///@brief Unpickles flow information data from a binary file.
	///@param filename String of the binary file to be read.
		/// @author SMM
    /// @date 01/016/12
  void unpickle(string filename);
  ///@brief Pickles flow information data from a binary file.
  ///@details WARNING!!! This creates HUGE files (sometimes 10x bigger than
  /// original file). Testing indicates reading this file takes
  /// almost as long as recalculating the flowinfo object so
  /// is probably not worth doing
	///@param filename String of the binary file to be written.
  /// @author SMM
  /// @date 01/016/12
  void pickle(string filename);

	// functions for getting flow, discharge, sediment flux, etc

	///@brief This function calculates the contributing pixels.
	///It can be converted to contributing area by multiplying by the
  ///DataResolution^2. In this function a pixel that has no donors has a
  ///contributing pixel value of 0.
  ///@return LSDIndexRaster of upslope contributing pixels.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster calculate_n_pixels_contributing_from_upslope();

	///@brief This calculates area and makes an index into the s vector for
  ///efficient calculation of the basin upslope of a given node.
  /// @author SMM
  /// @date 01/016/12
  void calculate_upslope_reference_indices();

	// algorithms for basin collection
	///@brief This function returns the base level node with the greatest
  ///drainage area.
  ///@return Integer node index.
  int retrieve_largest_base_level();
  ///@brief This function returns an integer vector containing all the node
  ///indexes upslope of of the node with number node_number_outlet.
  ///@param node_number_outlet Integer of the target node.
  ///@return Integer vector of upslope node indexes.
  /// @author SMM
  /// @date 01/016/12
	vector<int> get_upslope_nodes(int node_number_outlet);

	// algorithms for stream profile analysis

	///@brief This function calculates the chi function for all the nodes upslope
  ///of a given node.
  ///@param starting_node Integer index of node to analyse upslope of.
  ///@param m_over_n
  ///@param A_0
  ///@return Vector of chi values.
	vector<double> get_upslope_chi(int starting_node, double m_over_n, double A_0);
  ///@brief This function calculates the chi function for all the nodes upslope
  ///of a given list of nodes.
  ///@param upslope_pixel_list Vector of nodes to analyse.
  ///@param m_over_n
  ///@param A_0
  ///@return Vector of chi values.
   /// @author SMM
  /// @date 01/016/12
  vector<double> get_upslope_chi(vector<int>& upslope_pixel_list, double m_over_n, double A_0);

	///@brief Calculates the distance from outlet of all the base level nodes.
	///Distance is given in spatial units, not in pixels.
  ///@return LSDRaster of the distance to the outlet for all baselevel nodes.
  /// @author SMM
  /// @date 01/016/12
	LSDRaster distance_from_outlet();

	///@brief A get sources version that uses the flow accumulation pixels.
	///@param FlowPixels LSDIndexRaster of flow accumulation in pixels.
	///@param threshold Integer flow accumulation threshold.
  ///@return Vector of source integers.
   /// @author SMM
  /// @date 01/016/12
	vector<int> get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold);

	protected:

	///Number of rows.
  int NRows;
  ///Number of columns.
	int NCols;
	///Minimum X coordinate.
  double XMinimum;
	///Minimum Y coordinate.
	double YMinimum;

	///Data resolution.
	double DataResolution;
	///No data value.
	int NoDataValue;

	/// The number of nodes in the raster that have data.
	int NDataNodes;

	/// An array that says what node number is at a given row and column.
	Array2D<int> NodeIndex;

  /// @brief A raster of flow direction information.
  ///
	/// In the format:
	///
	/// 7  0 1 \n
	/// 6 -1 2 \n
	/// 5  4 3 \n
	///
	/// Nodes with flow direction of -1 drain to themselvs and are base level/sink nodes.
  Array2D<int> FlowDirection;

  /// @brief A code to denote the flow length from the node to its reciever node.
  /// <b>Each node has one and only one receiver.</b>
  /// \n\n
	/// 0 == no receiver/self receiver (base level) \n
	/// 1 == cardinal direction, flow length = DataResolution \n
	/// 2 == diagonal, flow length = DataResolution*(1/sqrt(2)) \n
	Array2D<int> FlowLengthCode;

	/// @brief This stores the row of a node in the vectorized
	/// node index. It, combined with ColIndex, is the
	/// inverse of NodeIndex.
	vector<int> RowIndex;

  /// @brief This stores the column of a node in the vectorized
  /// node index. It, combined with RowIndex, is the
  /// inverse of NodeIndex.
	vector<int> ColIndex;

  /// A list of base level nodes.
	vector<int> BaseLevelNodeList;

  /// Stores the number of donors to each node.
	vector<int> NDonorsVector;

  /// Stores the node index of the receiving node.
	vector<int> ReceiverVector;

  /// @brief Stores the delta vector which is used to index into the donor stack
  /// and order contributing nodes. See Braun and Willett (2012).
  vector<int> DeltaVector;

  /// This is a vector that stores the donor nodes of of the nodes and is
  /// indexed by the DeltaVector.
  vector<int> DonorStackVector;

  /// @brief This vector is used to caluculate flow accumulation. For each base
  /// level node it progresses from a hilltop to a confluence and then jumps to
  /// the next hilltop so that by cascading down through the node indices in
  /// this list one can quickly calculate drainage area, discharge, sediment
  /// flux, etc.
	vector<int> SVector;

  /// This stores the base level node for all of the nodes in the DEM.
	vector<int> BLBasinVector;

  /// This points to the starting point in the S vector of each node.
	vector<int> SVectorIndex;

  /// @brief The number of contributing nodes <b>INCULDING SELF</b> to a current
	/// pixel. It is used in conjunction with the SVectorIndex to build
	/// basins upslope of any and all nodes in the node list.
	vector<int> NContributingNodes;

  /// @brief Boundary conditions stored in a vector of four strings.
	/// The conditions are North[0] East[1] South[2] West[3].
	///
	/// There are 3 kinds of edge boundaries: no flux, base level and periodic.
	///
	/// The strings can be any length, as long as the first letter corresponds to the
	/// first letter of the boundary condition. It is not case sensitive.
  vector<string> BoundaryConditions;

	private:
	void create();
	void create(string fname);
	void create(vector<string> temp_BoundaryConditions, LSDRaster& TopoRaster);
};

#endif
