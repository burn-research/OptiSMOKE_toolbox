/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

// Order of the basis functions, which refers to the order of the spline that will represent it 
//which is degree + 1, or also the number of knots that will be used to describe the spline
int m = 4;

/* Degree of the basis functions */
int g = 3;

/* Orders of magnitude of difference between the smallest and the largest
possible value of the smoothing parameter lambda */
double lambdaSearchInterval = 6;

/* Number of steps in the for cycle for minimizing the smoothing parameter
lambda */
int numberOfStepsLambda = 13;

/* Fraction of the range of a spline on the y-axis for determining which
segments of the spline count as asymptotes. If the oscillations of the spline
at one of its extremities are contained within a horizontal area with size
determined by this value, the corresponding segment is identified as an
asymptote */
double fractionOfOrdinateRangeForAsymptoteIdentification = 0.005;

/* Fraction of the range of a spline on the y-axis for determining which points
count as well-defined maxima. In order to be considered a well-defined maximum,
a point in a spline must not only have first derivative equal to 0 and negative
second derivative, it must also be sufficiently distant from the two surrounding
minima. The minimum admissible distance is determined using this variable */
double fractionOfOrdinateRangeForMaximumIdentification = 0.025;

/* Fraction of the range of the experimental data on the y-axis for determining
whether a model can be considered a flat line with ordinate = 0 when compared to
the experimental data. This is useful for identifying situations that are
problematic for the calculation of the similarity indexes */
double fractionOfExpHeightForZeroLineIdentification = 0.02;

/* Minimum range on the x-axis of the models, before the addition of segments at
the extremes, required for comparison with the experimental data, expressed as a
fraction of the range of the experimental data on the x-axis */
double fractionOfExpRangeForModelExtrapolation = 0.5;

/* Fraction of the range of the experimental data on the x-axis used to
calculate the minimum shift possible */
double fractionOfExpRangeForMinShift = 0.005;

/* Fraction of the range of the experimental data on the x-axis used to shift
the models around the position of perfect alignment between their well-defined
maximum and the well-defined maximum of the experimental data */
double fractionOfExpRangeForShiftAroundMaximum = 0.05;

/* Specifies whether negative segments on the y-axis are admissible for the
splines or whether they should be replaced with straight lines with ordinate 0
*/
bool possibleNegativeOrdinates = false;

/* Number of plausible versions of the experimental data generated during the
bootstrap procedure. numberOfBootstrapVariations = 1 : no bootstrap */
int numberOfBootstrapVariations = 20;

/* Specifies whether the program should attempt to line up the maxima of
experimental data and models during the shift procedure */
bool lineUpMaxima = true;

/* Specifies whether the program should use the sum of the four dissimilarity
indexes when calculating the shift, thus finding a single shift amount, or
whether the program should consider each of the four indexes separately, thus
obtaining four different values for the shift */
bool useSumOfIndexesForAlignment = true;

/* Number of trapezoids for the numerical calculation of the indexes */
int numberOfTrapezoids = 99;

/* Default value for the relative experimental error. Used in case relative
errors are not provided along with the experimental data */
double defaultRelativeError = 0.1;

/* Specifies whether the program should use the similarity indexes to choose the
best set of nodes for the experimental data spline */
bool useIndexesToChooseExpNodes = false;

/* Specifies whether the individual indexes (d0L2, d1L2, d0Pe, d1Pe, shift)
should be saved to .csv files */
bool saveIndexesToCsv = false;

/* Specifies whether the program should save the sequences of points
representing the splines, used for plotting with external software (e.g. Excel)
*/
bool saveGraphData = false;

/* Specifies whether to create the graphs with the points and the knots used for
the fitting, and the resulting spline */
bool graphsFittingD0 = true;

/* Specifies whether to create the graphs comparing the spline for the
experimental data with the splines for the models */
bool graphsD0 = true;

/* Specifies whether to create the graphs comparing the first derivative of the
spline for the experimental data with the first derivatives of the splines for
the models */
bool graphsD1 = false;

/* Specifies whether to create the graphs comparing the spline for the
experimental data with the splines for the shifted models */
bool graphsD0Shift = true;

/* Specifies whether to create the graphs comparing the first derivative of the
spline for the experimental data with the first derivatives of the splines for
the shifted models */
bool graphsD1Shift = false;

/* Specifies whether to create the graphs with the various versions of the
experimental data created with the bootstrap procedure */
bool graphsBootstrap = false;

/* Specifies whether at least one kind of graph is to be created */
bool graphs = true;

/* Number of points to be calculated for each spline when saving the spline to a
.R file or to a .txt for future plotting */
int graphPoints = 500;

/* Graph quality. Specifies the height of the .png files generated, in number of
pixels. 1000: standard quality, 2000: high quality */
int graphQuality = 1000;

/* Pascal's triangle */
std::vector<std::vector<double>> pascalsTriangle;
