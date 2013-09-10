//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRaster
// Land Surface Dynamics Raster
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
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


/** @file LSDRaster.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0.0
@brief Main analysis object to interface with other LSD objects.
@details This object contains a diverse range of geomophological
analysis routines which can be used in conjunction with the other objects in
the package.

<b>change log</b>an
MASSIVE MERGE: Starting version 1.0.0 on 15/07/2013

@date 16/07/2013
*/

/**
@mainpage
This  is the documentation for Edinburgh Topographic Analysis Package (ETAP),
incorporating LSDRaster.

These pages will help you get started using this software.

\image html ./logo.png


Tools are included to:
- Generate topographic metrics
- Perform Chi analysis
- And other important science stuff
.

@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRaster_H
#define LSDRaster_H

///@brief Main analysis object to interface with other LSD objects.
class LSDRaster
{
	public:
	// declare the LSDFlowInfo object to be a friend class
	// this gives the LSDFlowInfo object access to the data elements
	// in the LSDRaster
	/// @brief Object to perform flow routing.
	friend class LSDFlowInfo;

  /// @brief The create function. This is default and throws an error.
	LSDRaster()										{ create(); }
	/// @brief Create an LSDRaster from a file.
  /// Uses a filename and file extension
  /// @return LSDRaster
  /// @param filename A String, the file to be loaded.
  /// @param extension A String, the file extension to be loaded.
	LSDRaster(string filename, string extension)	{ create(filename, extension); }
	/// @brief Create an LSDRaster from memory.
  /// @return LSDRaster
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A double of the minimum X coordinate.
  /// @param ymin A double of the minimum Y coordinate.
  /// @param cellsize A double of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param data An Array2D of doubles in the shape nrows*ncols,
  ///containing the data to be written.
  LSDRaster(int nrows, int ncols, double xmin, double ymin,
	          double cellsize, double ndv, Array2D<double> data)
								{ create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }

	// Get functions

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
	/// @return Raster values as a 2D Array.
	Array2D<double> get_RasterData() const { return RasterData.copy(); }

	/// Assignment operator.
	LSDRaster& operator=(const LSDRaster& LSDR);


	/// @brief Read a raster into memory from a file.
  ///
  /// The supported formats are .asc and .flt which are
  /// both exported and imported by arcmap.
	///
	/// The filename is the string of characters before the '.' in the extension
	/// and the extension is the characters after the '.'.
	///
	/// If the full filename is my_dem.01.asc then:
  /// filename = "my_dem.01" and extension = "asc".
  ///
	///
	/// For float files both a data file and a header are read
	/// the header file must have the same filename, before extention, of
	/// the raster data, and the extension must be .hdr.
	///
	/// @author SMM
  /// @date 01/01/12
  void read_raster(string filename, string extension);

  /// @brief Read a raster from memory to a file.
  ///
  /// The supported formats are .asc and .flt which are
  /// both exported and imported by arcmap.
	///
	/// The filename is the string of characters before the '.' in the extension
	/// and the extension is the characters after the '.'.
	///
	/// If the full filename is my_dem.01.asc then:
  /// filename = "my_dem.01" and extension = "asc".
	///
	/// For float files both a data file and a header are written
	/// the header file must have the same filename, before extention, of
	/// the raster data, and the extension must be .hdr.
	///
	/// @param filename a string of the filename _without_ the extension.
	/// @param extension a string of the extension _without_ the leading dot
	/// @author SMM
  /// @date 01/01/12
  void write_raster(string filename, string extension);

  /// @brief Calculate the minimum bounding rectangle for an LSDRaster Object and crop out
  /// all the surrounding NoDataValues to reduce the size and load times of output rasters.
  ///
  /// @details Ideal for use with chi analysis tools which output basin and chi m value rasters
  /// which can be predominantly no data. As an example, a 253 Mb file can be reduced to
  /// ~5 Mb with no loss or resampling of data.
  ///
  /// @return A trimmed LSDRaster object.
  /// @author SWDG
  /// @date 22/08/13
  LSDRaster RasterTrimmer();

  /// @brief Make LSDRaster object using a 'template' raster and an Array2D of data.
  /// @param InputData 2DArray of doubles to be written to LSDRaster.
  /// @return LSDRaster containing the data passed in.
  /// @author SWDG
  /// @date 29/8/13
  LSDRaster LSDRasterTemplate(Array2D<double> InputData);

	// Functions relating to shading, shadowing and shielding

  /// @brief This function generates a hillshade raster.
  ///
  /// It uses the the algorithm outlined in Burrough and McDonnell Principles
  /// of GIS (1990) and in the ArcMap web help
  /// http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/
  /// spatial_analyst_tools/how_hillshade_works.htm
  ///
  /// Default values are altitude = 45, azimuth = 315, z_factor = 1
  /// @param altitude of the illumination source in degrees.
  /// @param azimuth of the illumination source in degrees
  /// @param z_factor Scaling factor between vertical and horizontal.
  /// @return Hillshaded LSDRaster object
  /// @author SWDG
  /// @date February 2013
	LSDRaster hillshade(double altitude, double azimuth, double z_factor);



	/// @brief This function generates a topographic shielding raster using the algorithm outlined in Codilean (2006).
  ///
  /// @details Creating a raster of values between 0 and 1 which can be used as a
  /// scaling factor in Cosmo analysis.
  ///
  /// Goes further than the original algorithm allowing a theoretical theta,
  /// phi pair of 1,1 to be supplied and although this will increase the
  /// computation time significantly, it is much faster than the original
  /// Avenue and VBScript implementations.
  ///
  /// Takes 2 ints, representing the theta, phi paring required.
  /// Codilean (2006) used 5,5 as the standard values, but in reality values of
  /// 10,15 are often preferred to save processing time.
  /// @param theta_step Spacing of sampled theta values.
  /// @param phi_step Spacing of sampled phi values.
  /// @pre phi_step must be a factor of 360.
  /// @author SWDG
  /// @date 11/4/13
	LSDRaster TopoShield(int theta_step, int phi_step);

  /// @brief This looks for isolated instances of no data.
  ///
  /// Does nothing else but print their location to the screen.
  /// @author MDH, DTM
  /// @date 01/01/12
  void check_isolated_nodata();

  /// @brief Get the raster data at a specified location.
  /// @param row An integer, the X coordinate of the target cell.
  /// @param column An integer, the Y coordinate of the target cell.
  /// @return The raster value at the position (row, column).
  /// @author SMM
  /// @date 01/01/12
  double get_data_element(int row, int column)	{ return RasterData[row][column]; }

  // this calculates coefficeint matrices for calculating a variety of
  // surface metrics such as slope, aspect, curvature, etc.

  /// @brief This function calculates 6 coefficient matrices that allow the user to
  /// then calcualte slope, curvature, aspect, a classification for finding saddles and peaks
  /// and other metrics.
  ///
  /// @details The coefficient matrices are overwritten during the running of this member function.
  ///
	/// Have N simultaneous linear equations, and N unknowns.
	/// => b = Ax, where x is a 1xN array containing the coefficients we need for
	/// surface fitting.
	/// A is constructed using different combinations of x and y, thus we only need
	/// to compute this once, since the window size does not change.
	/// For 2nd order surface fitting, there are 6 coefficients, therefore A is a
	/// 6x6 matrix.
  /// Updated 15/07/2013 to use a circular mask for surface fitting. DTM
  /// @param window_radius Radius of the mask.
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @param c coefficeint c.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @param f coefficeint f.
  /// @author DTM, SMM
  /// @date 01/01/12
  void calculate_polyfit_coefficient_matrices(double window_radius,
									Array2D<double>& a, Array2D<double>& b,
									Array2D<double>& c, Array2D<double>& d,
									Array2D<double>& e, Array2D<double>& f);

	// a series of functions for retrieving derived data from the polyfit calculations

  /// @brief  This function calculates the elevation based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param f coefficeint f.
  /// @return LSDRaster of elevations.
  /// @author FC
  /// @date 24/03/13
  LSDRaster calculate_polyfit_elevation(Array2D<double>& f);
  /// @brief  This function calculates the slope based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of slope.
  /// @author DTM, SMM
  /// @date 01/01/12
	LSDRaster calculate_polyfit_slope(Array2D<double>& d, Array2D<double>& e);
  /// @brief  This function calculates the aspect based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of aspect.
  /// @author DTM, SMM
  /// @date 01/01/12
	LSDRaster calculate_polyfit_aspect(Array2D<double>& d,Array2D<double>& e);
	/// @brief  This function calculates the curvature based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @return LSDRaster of curvature.
  /// @author DTM, SMM
  /// @date 01/01/12
	LSDRaster calculate_polyfit_curvature(Array2D<double>& a,Array2D<double>& b);
	/// @brief  This function calculates the planform curvature based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @param c coefficeint c.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of planform curvature.
  /// @author DTM, SMM
  /// @date 01/01/12
  LSDRaster calculate_polyfit_planform_curvature(Array2D<double>& a, Array2D<double>& b, Array2D<double>& c,
  													Array2D<double>& d, Array2D<double>& e);
	/// @brief  This function calculates the profile curvature based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @param c coefficeint c.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of profile curvature.
  /// @author DTM, SMM
  /// @date 01/01/12
  LSDRaster calculate_polyfit_profile_curvature(Array2D<double>& a, Array2D<double>& b, Array2D<double>& c,
  													Array2D<double>& d, Array2D<double>& e);
	/// @brief  This function calculates the tangential curvature based on a polynomial fit.
  ///
  /// @details the window is determined by the calculate_polyfit_coefficient_matrices
  /// this function also calculates the a,b,c,d,e and f coefficient matrices.
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @param c coefficeint c.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of tangential curvature.
  /// @author DTM, SMM
  /// @date 01/01/12
	LSDRaster calculate_polyfit_tangential_curvature(Array2D<double>& a, Array2D<double>& b, Array2D<double>& c,
  													Array2D<double>& d, Array2D<double>& e);


  /// @brief This function identifies approximate position of stationary points within
  /// discrete surface using a threshold slope.
  ///
  /// @details The nature of the stationary point is then determined to discriminate
  /// peaks, depressions and saddles. \n
  /// 0 = Non-stationary  \n
  /// 1 = Peak             \n
  /// 2 = Depression        \n
  /// 3 = Saddle             \n
  /// @param a coefficeint a.
  /// @param b coefficeint b.
  /// @param c coefficeint c.
  /// @param d coefficeint d.
  /// @param e coefficeint e.
  /// @return LSDRaster of classified elevation data.
  /// @author DTM
  /// @date	17/09/2012
	LSDIndexRaster calculate_polyfit_classification(Array2D<double>& a, Array2D<double>& b, Array2D<double>& c,
	                                                Array2D<double>& d, Array2D<double>& e);

	/// @brief this function takes the polyfit functions and requires a window radius and a vector telling the
	/// function which rasters to print to file.
  ///
  /// @details The function is data efficient since one does not need to
  /// recalculate the polyfit coefficeint matrices. It also takes a string
  /// which is the prename of the data files the file codes in the vector are:\n
	/// 0 slope         \n
	/// 1 aspect        \n
	/// 2 curvature       \n
	/// 3 planform curvature\n
	/// 4 profile curvature   \n
	/// 6 tangential curvature  \n
	/// 7 classification          \n
  /// @param window_radius Radius of the mask.
  /// @param file_prefix Output filename string.
  /// @param file_list Vector of files to be created.
  /// @author SMM
  /// @date 19-12-2012
	void calculate_and_print_polyfit_rasters(double window_radius,
											string file_prefix, vector<int> file_list);

	/// @brief Gets the hilltop curvature raster.
	///
	/// @details Modified to take an LSDRaster of hilltops - SWDG 29/8/13
  ///
  /// @param curvature LSDRaster of curvatures.
  /// @param Hilltops LSDRaster of hilltops.
  /// @return LSDRaster of hilltop curvatures.
  /// @author DTM
  /// @date 30/04/13
	LSDRaster get_hilltop_curvature(LSDRaster& curvature, LSDRaster& Hilltops);



  	// Rock exposure index
    /// @brief This function is a wrapper to get the three roughness eigenvalues
    /// s1, s2 and s3.
    /// @param window_radius
    /// @param a_plane
    /// @param b_plane
    /// @param c_plane
    /// @author DTM
    /// @date 15/7/2013
  	void calculate_plane_coefficient_matrices(double window_radius, Array2D<double>& a_plane,
										Array2D<double>& b_plane, Array2D<double>& c_plane);
		/// @brief Create the REI raster
    ///
    /// @details
    /// @param a_plane
    /// @param b_plane
    /// @param CriticalSlope
    /// @return LSDIndexRaster of rock exposure.
    /// @author DTM
  	LSDIndexRaster calculate_REI(Array2D<double>& a_plane, Array2D<double>& b_plane, double CriticalSlope);




	// hydrology tools
	///@brief This function fills pits/sinks in a DEM by incrementing elevations for cells with
  ///no downslope neighbour. The process is repeated adnausium until no cells require
  ///incrementing.
  ///
  ///Inputs required are a DEM file in ascii raster format as created by ARCMap
  ///and a file name to create a filled DEM grid.
  ///
  ///This code was built ontop of code made available by Jon D. Pelletier as part
  ///of his book:
  ///
  ///Pelletier,J.D.,'Quantitative Modelling of Landscapes' Cambridge University Press
  ///
  ///---------------------------------------------------------------------------------
  ///
  /// v1.3 reduced fill increment to 1mm  to avoid 'overfilling'
  ///
  /// Martin Hurst, October 2011
  ///
  ///---------------------------------------------------------------------------------
  ///
  /// v1.2 modified to read *.flt files
  ///
  /// Martin Hurst, November 2010
  ///
  ///---------------------------------------------------------------------------------
  ///
  /// v1.1 function incorporated to allow the tool to fill adjacent pixels immediately
  /// after filling a given pixel, should speed things up.
  ///
  /// Martin Hurst, October 2010
  ///
  ///---------------------------------------------------------------------------------
  ///
  /// v1.0 is slow as it requires many iterations through the dem
  ///
  /// Martin Hurst, June 2010
  /// @return Filled LSDRaster.
  /// @author MDH
  /// @date 01/06/10
  LSDRaster fill();
	/// @brief This is a recursive algorithm that is called by the fill function.
  /// @param fill_data
  /// @param i
  /// @param j
  /// @author MDH
  /// @date 01/06/10
	void fill_iterator(Array2D<double>& fill_data, int i, int j);


  /// @brief This function fills pits/sinks in a DEM by checking for pits from
  ///lowest to highest elevation, starting at the DEM boundary (raster edge or
  /// adjacent to NDVs).
  ///
  /// @details Utilises a priority queue to progressively populate the stack and
  /// pop out the the lowest value before checking that the neighbouring cells
  /// that are yet to be visited must be higher in a hydrologically correct DEM.
  /// This method is substantially faster on datasets with pits consisting of
  /// multiple cells since each cell only needs to be visited once.
  ///
  /// Method taken from Wang and Liu (2006), Int. J. of GIS. 20(2), 193-213
  /// @param MinSlope The minimum slope between two Nodes once filled. If set
  /// to zero will create flats.
  /// @return Filled LSDRaster object.
  /// @author Martin Hurst
  /// @date 12/3/13
  LSDRaster fill(double& MinSlope);





	/// @brief Generate a flow area raster using a multi direction algorithm.
  ///
  /// @details Computes the proportion of all downslope flows for each cell in the input
  /// DEM, and weights them using the equation from Freeman et al 1991 and routes the
  /// flow accordingly.
  ///
  /// Paper link: http://www.sciencedirect.com/science/article/pii/009830049190048I
  ///
  /// Cardinal Weighting = (elevation_drop/total_elevation_drop)^1.1  \n
  /// Diagonal Weighting = ((elevation_drop/total_elevation_drop)*(1/root(2)))^1.1
  ///
  /// Can <b>NOT</b> handle DEMs containing flats or pits -  must be filled using the new
  /// LSDRaster fill. Function built around original c++ code by Martin Hurst.
  /// @return LSDRaster of flow area.
  /// @author SWDG
  /// @date 18/4/13
  LSDRaster FreemanMDFlow();

  /// @brief Generate a flow area raster using a multi direction algorithm.
  ///
  /// @details Computes the proportion of all downslope flows for each cell in the input
  /// DEM, and weights them using the equation from Quinn et al 1991 and routes the
  /// flow accordingly.
  ///
  /// Paper link: http://onlinelibrary.wiley.com/doi/10.1002/hyp.3360050106/abstract
  ///
  /// Cardinal Weighting = (elevation_drop/total_elevation_drop)*DataResolution/2 \n
  /// Diagonal Weighting = ((elevation_drop/total_elevation_drop)*(1/root(2)))* DataResolution*0.354
  ///
  /// Can <b>NOT</b> handle DEMs containing flats or pits -  must be filled using the new
  /// LSDRaster fill. Function built around original c++ code by Martin Hurst.
  /// @return LSDRaster of flow area.
  /// @author SWDG
  /// @date 18/4/13
	LSDRaster QuinnMDFlow();

	/// @brief Generate a flow area raster using a multi 2-direction algorithm.
  ///
  /// @details Computes the proportion of all downslope flows for each cell in the input
  /// DEM. Finds the cell of the steepest descent and then checks the two
  /// neighbouring cells slopes. If either is also downslope proportion flow
  /// between the steepest cell and the steepest neighbour. If neither neighbour
  /// is downslope 100% of flow follows the steepest path.
  ///
  /// Can <b>NOT</b> handle DEMs containing flats or pits -  must be filled using the new
  /// LSDRaster fill. Function built around original c++ code by Martin Hurst.
  /// @return LSDRaster of flow area.
  /// @author SWDG
  /// @date 02/08/13
  LSDRaster M2DFlow();



	// some tools associated with ridgeline analyis

	/// @brief Module to sample LSDRaster values running along a ridgetop network.
  ///
	/// @details Ridge network is generated from LSDChannelNetwork::ExtractRidges.
  /// @param Ridges 2D Array of ridge lines.
  /// @return Sampled LSDRaster object.
  /// @author SWDG
  /// @date 04/2013
	LSDRaster RidgeSample(Array2D<double>& Ridges);
	/// @brief Pass a smoothing window over a ridge LSDRaster object to calculate
  /// an average value running along the ridgetop.
  /// @param WindowRadius optional integer smoothing radius between 1 and 6.
  /// @pre Default smoothing radius is 2 and will revert to that value if a
  /// radius outside the vaid range (1 to 6) is passed.
  /// @return Averaged LSDRaster object.
  /// @author SWDG
  /// @date 04/2013
	LSDRaster RidgeSmoother(int WindowRadius);
	/// @brief Pass a buffer over a ridge LSDRaster object to increase sampling area.
  ///
  /// @details Buffers equally in all directions, so use with care to avoid
  /// sampling areas away from the axis of the original ridge line.
  /// @param BufferRadius optional integer buffer radius between 1 and 6.
  /// @pre Default buffer radius is 2 and will revert to that value if a
  /// radius outside the vaid range (1 to 6) is passed.
  /// @return LSDRaster object containing buffered ridges.
  /// @author SWDG
  /// @date 04/2013
	LSDRaster RidgeBuffer(int BufferRadius);
	/// @brief Module assigns an average LSDRaster value to each basin.
  ///
  /// @details Works by searching for every hilltop value that falls within a
  /// basin, summing these values and writing the final average to every cell
  /// identified as the basin in question.
  ///
  /// Very inefficent at present. Module loops through every cell in LSDRaster
  /// (2 * number of basins) + 1 times. Beware!
  /// @param Basins LSDIndexRaster of Drainage basins, generated using
  /// ChannelNetwork::ExtractBasinsOrder.
  /// \n\n Bug fixed in assignment of basin IDs - SWDG 2/9/13.
  /// @return LSDRaster of average basin value for each identified basin.
  /// @author SWDG
  /// @date 04/2013
	LSDRaster BasinAverager(LSDIndexRaster& Basins);

  /// @brief Write the area(in units of area) of each basin to the basin's pixels.
  /// @param Basins LSDIndexRaster of drainage basins to measure.
  /// @return LSDRaster of basin areas.
  /// @author SWDG
  /// @date 04/2013
  LSDRaster BasinArea(LSDIndexRaster& Basins);

  /// @brief Write hilltop metrics to text file.
  ///
  /// @details This can probably be absorbed by the main hilltop flow routing as all this does is write a text file
  /// with the hilltop pixels coded by basin id.
  /// @param Hilltops
  /// @param Basins
  /// @param HilltopRelief
  /// @param HilltopAspect
  /// @param HilltopSlope
  /// @param HilltopLength
  /// @param HilltopCurvature
  /// @author SWDG
  /// @date 27/8/13
  void BasinHilltopWriter(LSDRaster& Hilltops, LSDIndexRaster& Basins, Array2D<double>& HilltopRelief, Array2D<double>& HilltopAspect,
                                   Array2D<double>& HilltopSlope, Array2D<double>& HilltopLength, Array2D<double>& HilltopCurvature);

  /// @brief Punch basins out of an LSDRaster to create DEMs of a single catchment.
  ///
  /// @details Writes files in the user supplied format (flt or asc) and returns a vector
  /// of their filenames so they can be loaded into other functions.
  /// @param basin_ids Vector of basins to punch out.
  /// @param BasinArray Basin outlines used to punch out the LSDRasters.
  /// @param output_format The output file format.
  /// @param raster_prefix A prefix to name the output files with.
  /// @return Vector of outout filenames.
  /// @author SWDG
  /// @date 27/8/13
  vector<string> BasinPuncher(vector<int> basin_ids, LSDIndexRaster BasinArray, string output_format, string raster_prefix);




  /// @brief Generate data in two text files to create a boomerang plot as in Roering et al [2007].
  /// @param Slope LSDRaster of slope.
  /// @param D_inf D-infinity Flowarea LSDRaster.
  /// @param RasterFilename Filename used to give unique name to output data.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @author SWDG
  /// @date 27/8/13
  void Boomerang(LSDRaster& Slope, LSDRaster& D_inf, string RasterFilename, double log_bin_width = 0.1);

  /// @brief Calculate drainage density of a set of input basins.
  ///
  /// @details Calculated as flow length/basin area and written to every
  /// cell of the identified basin.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param Basins LSDIndexRaster of the basins to be analysed.
  /// @param FlowDir Array2D of flowdirections generated by FlowInfo.get_FlowDirection().
  /// @return LSDRaster of basins coded with drainage density.
  /// @author SWDG
  /// @date 04/2013
  LSDRaster DrainageDensity(LSDIndexRaster& StreamNetwork, LSDIndexRaster& Basins, Array2D<int> FlowDir);

	// Smoothing tools
	//Nonlocal Means Filtering - Default values following Baudes et al [2005]
	/// @brief Perform Non-local means filtering on a DEM following Baude et al [2005].
  ///
  /// @details Smoothes non-gaussian noise. Martin Hurst, February, 2012
	/// Modified by David Milodowski, May 2012- generates grid of recording filtered noise
  ///
	/// WindowRadius has to be <= SimilarityRadius ?
  ///
	/// Adapted from a matlab script by:
	/// Author: Jose Vicente Manjon Herrera & Antoni Buades
	/// Date: 09-03-2006
  ///
	/// Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
	/// "A non-local algorithm for image denoising"
  ///
	/// **Added soft threshold optimal correction - David Milodowski, 05/2012
  /// @param WindowRadius search window radius (defualt=2).
  /// @param SimilarityRadius similarity window radius (defualt=2).
  /// @param DegreeFiltering degree of filtering (defualt=2).
  /// @param Sigma (default=1.25).
  /// @return Filtered LSDRaster object.
  /// @author MDH, DTM
  /// @date February 2012
	LSDRaster NonLocalMeansFilter(int WindowRadius=2, int SimilarityRadius=2, int DegreeFiltering=2, double Sigma=1.25);
	/// @brief Creates a buffer around an array (of size SimilarityRadius) and
  /// gives the new border mirror symmetric values of the original array reflected
  /// across the boundary.
  ///
  /// @details SimilarityRadius should be the size of the window if filtering.
	/// New array has size nrows + 2*SimilarityRadius by ncols + 2*SimilarityRadius.
  /// @param PaddedRasterData Padded LSDRaster.
  /// @param SimilarityRadius similarity window radius.
  /// @author Martin Hurst
  /// @date February 2012
	void PadRasterSymmetric(Array2D<double>& PaddedRasterData, int& SimilarityRadius);
	/// @brief Generate gaussian weighted kernel.
  ///
  /// @details kernel array must be predeclared of size SimilarityRadius and consist of zeros:
	/// Array2D<double> Kernel(SimilarityRadius,SimilarityRadius,0.0);
	///
	/// Kernel generated using: G(x,y) = (1/2*pi*sigma^2) exp ((-x^2+y^2)/(2*sigma^2))
  /// @param Kernel
  /// @param sigma
  /// @param SimilarityRadius similarity window radius.
  /// @author Martin Hurst
  /// @date February 2012
	void MakeGaussianKernel(Array2D<double>& Kernel, double sigma, int SimilarityRadius);

  //D-infinity tools

  /// @brief D-infinity flow direction algorithm after Tarboton (1997).
  ///
  /// @details Algorithm takes a filled DEM and for each cell calculates the steepest descent
  /// based on 8 triangular facets. Flow direction is assigned as an angle from 0-360
  /// degrees with -1 used to flag unresolved areas such as pits.
  ///
  /// Code is ported and optimised from a Java implementation of the algorithm
  /// supplied under the GNU GPL licence through WhiteBox GAT:
  /// http://www.uoguelph.ca/~hydrogeo/Whitebox/ and provides identical results
  /// to the whitebox tool.
  /// @return Array of Flow directions in degrees.
  /// @author SWDG
  /// @date 26/07/13
  Array2D<double> D_inf_FlowDir();

  /// @brief Main function for generating a D-infinity flow area raster after Tarboton (1997).
  ///
  /// @details Calls the recurisve D_infAccum function to get flow area for each pixel.
  /// Returns flow area in pixels.
  ///
  /// Code is ported and optimised from a Java implementation of the algorithm
  /// supplied under the GNU GPL licence through WhiteBox GAT:
  /// http://www.uoguelph.ca/~hydrogeo/Whitebox/ and provides identical results
  /// to the whitebox tool.
  /// @param FlowDir_array Array of Flowdirections generated by D_inf_FlowDir().
  /// @return LSDRaster of D-inf flow areas in pixels.
  /// @author SWDG
  /// @date 26/07/13
  LSDRaster D_inf_FlowArea(Array2D<double> FlowDir_array);

  /// @brief Recursive function to calculate accumulating area for a given pixel.
  ///
  /// @details Called by the driver for every cell which has no contributing
  /// cells - eg the highest points on the landscape. Avoids the need to flatten
  /// and sort the DEM as required in the original Tarboton (1997)
  /// implementation. For more detail on the recursive algorithm following
  /// channels see Mark (1998) "Network Models in Geomorphology".
  ///
  /// Code is ported and optimised from a Java implementation of the algorithm
  /// supplied under the GNU GPL licence through WhiteBox GAT:
  /// http://www.uoguelph.ca/~hydrogeo/Whitebox/ and provides identical results
  /// to the whitebox tool.
  ///
  /// @param i Row index of the target cell.
  /// @param j Column index of the target cell.
  /// @param CountGrid Array showing the number of cells flowing into each cell.
  /// @param Flowarea_Raster Array of the flow areas which is populated using this function.
  /// @param FlowDir_array Array of Flowdirections generated by D_inf_FlowDir().
  /// @author SWDG
  /// @date 26/07/13
  void D_infAccum(int i, int j, Array2D<double> CountGrid, Array2D<double> Flowarea_Raster, Array2D<double> FlowDir_array);

  /// @brief Wrapper Function to create a D-infinity flow area raster with one function call.
  /// @return LSDRaster of D-inf flow areas in pixels.
  /// @author SWDG
  /// @date 26/07/13
  LSDRaster D_inf();

  /// @brief Function to write the D-infinity flow directions to an LSDRaster.
  /// @param dinflow Array of Flowdirections generated by D_inf_FlowDir().
  ///@return LSDRaster of D-inf flow directions in degrees.
  /// @author SWDG
  /// @date 26/07/13
  LSDRaster write_dinf_flowdir_to_LSDRaster(Array2D<double> dinflow);

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

	/// Raster data.
	Array2D<double> RasterData;

	private:
	void create();
	void create(string filename, string extension);
	void create(int ncols, int nrows, double xmin, double ymin,
	            double cellsize, double ndv, Array2D<double> data);

};

#endif
