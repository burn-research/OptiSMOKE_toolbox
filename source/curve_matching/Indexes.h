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

#include "Spline.h"
#include <fstream>

class Indexes {

public:

    /* Name of the _exp.txt and _mod.txt files with the data for the splines */
    std::string fileName;

    /* Number of models to be compared with the experimental data */
    int numberOfModels;

    /* Names of the models */
    std::vector<std::string> modelNames;

    /* Scores, before the shift. scores[i] contains the data relative to the
    comparison of modelNames[i] and the experimental data. scores[i] is equal to
    -1 in case it is not possible to calculate score i, or in case the
    calculated score is nan */
    std::vector<double> scores;

    /* Scores, after the shift. scores_shift[i] contains the data relative to
    the comparison of modelNames[i] and the experimental data. scores_shift[i]
    is equal to -1 in case it is not possible to calculate score i, or in case
    the calculated score is nan */
    std::vector<double> scores_shift;

    /* Standard deviations of the scores, before the shift. stdDev[i] contains
    the data relative to the comparison of modelNames[i] and the experimental
    data. stdDev[i] is equal to -1 in case it is not possible to calculate the
    corresponding score, or it is nan */
    std::vector<double> stdDev;

    /* Standard deviations of the scores, after the shift. stdDev_shift[i]
    contains the data relative to the comparison of modelNames[i] and the
    experimental data. stdDev_shift[i] is equal to -1 in case it is not possible
    to calculate the corresponding score, or it is nan */
    std::vector<double> stdDev_shift;

    /* Individual indexes for each variation of the experimental data obtained
    with the bootstrapping procedure. allIndexes[i][j][k]: i: variations of the
    experimental data; j: models; k: d0L2, d1L2, d0Pe, d1Pe, shift (or shifts if
    all different) */
    std::vector<std::vector<std::vector<double>>> allIndexes;

    ////////////////////////////////////////////////////////////////////////////
    bool print_indexes;
    /* Executes all the operations for the comparison of the experimental data
    with the models */
    double solve(   bool calculateShift,
               		bool logScale,
               		Spline spline_of_experiments,
               		std::vector<double> valuesSim_receive,
	       		    std::vector<double> abscissae,
	    		    int a,
			    bool print_indexes,
			    bool print_splines,
                std::string quantity_of_interest );

    double FINAL_SCORE;
    ////////////////////////////////////////////////////////////////////////////

private:

    /* Equivalent of 'scores' for every bootstrap variation */
    std::vector<std::vector<double>> bootstrapScores;

    /* Equivalent of scores_shift for every bootstrap variation */
    std::vector<std::vector<double>> bootstrapScores_shift;

    /* x-coordinates and y-coordinates of the experimental data and the models.
    inputData[0] and inputData[1] contain, respectively, the x-coordinates and
    the y-coordinates of the experimental data, inputData[i] and inputData[i+1]
    contain, respectively, the x-coordinates and the y-coordinates of
    modelNames[i/2-1] */
    std::vector<std::vector<double>> inputData;

    /* Splines for the experimental data */
    //Spline splinesExp;

    /* Splines for the models. splinesMod[i] contains the spline for
    modelNames[i] */
    //Spline splinesMod;
    //std::vector<double> Sim_values;
    /* Splines for the models, before being extended */
    Spline splinesMod_original;

    /* Path to the folder containing the input data */
    std::string folderPath;

    /* Name of the folder containing the input data */
    std::string folderName;

    /* Relative errors for each experimental data point */
    std::vector<double> relativeErrors;

    /* Specifies whether to calculate the indexes for the shifted splines */
    bool calculateShift;

	/* Specifies whether to calculate the exponential (base e) of the input
	ordinates before computing the splines */
	bool logScale;

    /* Amounts to be added to the abscissae of the models to find the shifted
    abscissae. shiftAmounts[i][j] refers to splinesExp[i] and modelNames[j] */
    std::vector<std::vector<std::vector<double>>> shiftAmounts;

    /* Positions in shiftAmounts, indicating the shifts that maximise the mean
    of d0L2, d1L2, d0Pe and d1Pe for the models. maxIndexes[i][j] refers to
    shiftAmounts[i][j] */
    std::vector<std::vector<int>> maxIndexes;

    /* Index in the splinesMod std::vector of the model being calculated */
    int indexOfModelBeingCalculated;

    /* Height of the trapezoids used used during integration, divided by 2 */
    double halfHeight;

    /* Powers of the abscissae of the sides of the trapezoids used for the
    numerical calculation of the indexes, distributed evenly between the
    extremes of the experimental data */
    std::vector<std::vector<double>> xPowers;

    /* Ordinates of the experimental data spline, with yD0[i] corresponding to
    xPowers[i][1] */
    std::vector<double> yD0;

    /* Ordinates of the first derivative of the experimental data spline, with
    yD1[i] corresponding to xPowers[i][1] */
    std::vector<double> yD1;

    /* Scores without shifts. Sorted from lowest to highest. Is equal to -1 in
    case it is not possible to calculate a score */
    std::vector<double> scoresSorted;

    /* Standard deviations without shifts. Sorted from lowest to highest. Is
    equal to -1 in case it is not possible to calculate a score */
    std::vector<double> stdDevSorted;

    /* Model names, sorted as the elements of scoresSorted */
    std::vector<std::string> modelNamesSorted;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted */
    std::vector<int> scoresSortedLineTypes;

    /* Color indexes in the .R graphs, sorted as the elements of scoresSorted */
    std::vector<int> scoresSortedColorIndexes;

    /* Scores with shifts. Sorted from lowest to highest. Is equal to -1 in case
    it is not possible to calculate a score */
    std::vector<double> scoresSorted_shift;

    /* Standard deviations with shifts. Sorted from lowest to highest. Is equal
    to -1 in case it is not possible to calculate a score */
    std::vector<double> stdDevSorted_shift;

    /* Model names, sorted as the elements of scoresSorted_shift */
    std::vector<std::string> modelNamesSorted_shift;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted_shift
    */
    std::vector<int> scoresSortedLineType_shift;

    /* Color indexes in the .R graphs, sorted as the elements of
    scoresSorted_shift */
    std::vector<int> scoresSortedColorIndexes_shift;

    /* Index in the splinesExp std::vector of the spline being considered */
    int splineExpIndex;

    /* Index in the splinesExp std::vector of the best fit of the experimental data
    */
    int indexBestSplineExp;

    /* Values of the d0L2, d1L2, d0Pe and d1Pe indexes, in this order, obtained
    before any shift. indexes[i][j] refers to splinesExp[i] and contains the
    data relative to the comparison of modelNames[j] and the experimental data.
    Is equal to -1 in case it is not possible to calculate an index */
    std::vector<std::vector<std::vector<double>>> indexes;

    /* Values of the d0L2_shift, d1L2_shift, d0Pe_shift, d1Pe_shift indexes, in
    this order, corresponding to their optimal values found while shifting.
    indexes_shift[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate an index */
    std::vector<std::vector<std::vector<double>>> indexes_shift;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximize the means of d0L2, d1L2, d0Pe and d1Pe for each
    model. shifts[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate a shift */
    std::vector<std::vector<double>> shifts;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximise each of d0L2, d1L2, d0Pe and d1Pe for each model.
    shifts[i][j] refers to splinesExp[i] and contains the data relative to the
    comparison of modelNames[j] and the experimental data. Is equal to -1 in
    case it is not possible to calculate a shift */
    std::vector<std::vector<std::vector<double>>> shifts_allDifferent;

    /* Norm of the experimental data, calculated between the extremes of the
    experimental data */
    double normD0Exp;

    /* Norm of the first derivative of the experimental data, calculated between
    the extremes of the experimental data */
    double normD1Exp;

    /* Norm of a model, calculated between the extremes of the experimental data
    */
    double normD0Mod;

    /* Norm of the first derivative of a model, calculated between the extremes
    of the experimental data */
    double normD1Mod;

    /* Specifies whether for a model it is possible to calculate the indexes for
    the shifted model */
    std::vector<std::vector<bool>> possibleToCalculateShifts;

    /* Integral of the absolute value of the first derivative of the
    experimental data */
    double integralOfAbsoluteValueD1Exp;

    /* Index of the bootstrap variation currently being considered */
    int bootstrapIndex;

    /* Used by calculateShiftAmounts */
    double maxExtrapolatedLength;

    /* Used by calculateShiftAmounts */
    double shiftAroundMaximum;

    /* Used by calculateShiftAmounts */
    double distanceBetweenShiftedPoints;

    /* Used for selecting the ordinates of the graphs in 'Graphs Bootstrap' */
    Spline referenceExpSpline;

    ////////////////////////////////////////////////////////////////////////////

    /* Reads the data relative to the experiment (x, y, relative errors) and the
    models (names, x, y). Sorts the data in ascending order relative to the
    abscissae. If two or more abscissae are the same, replaces them with a
    single abscissa with the mean of the corresponding ordinates as ordinate,
    and with the mean of the corresponding relative errors as relative error.
    The data is read from folderPath/fileName_exp.txt and from
    folderPath/fileName_mod.txt */
    void readData();

    /* Calculates the splines for the experimental data */
    void calculateSplines();

    /* Calculates the indexes for each model */
    void calculateIndexes(Spline& spline_exp, Spline& spline_mod, bool print_indexes);

    /* Calculates the scores from the indexes */
    void calculateScores();

    /* Saves the data for creating the graph of each spline and the
    corresponding data points and knots to a .R file */
    //void graphsFittingR(std::vector<Spline>& splines);

    /* Saves the data for creating the graph comparing the experimental data
    with the models (either shifted or non-shifted) to a .R file. bootstrapIndex
    is the number of the variation of the experimental data in the graph among
    all the variations generated by the bootstrap procedure. bootstrapIndex
    ranges from 1 to numberOfBootstrapVariations. bootstrapIndex == 0: save
    graph to 'Scripts' and not to 'Graphs Bootstrap', using the original
    experimental data (whose bootstrapIndex is 1) */
    // void graphR(int derivativeOrder,
    //             bool shift,
    //             int bootstrapIndex,
    //             const std::string& folder);

    // /* Saves the data necessary to plot the splines to .txt files */
    // void saveGraphDataToFiles();

    /* Given a std::vector containing the sums of d0L2, d1L2, d0Pe and d1Pe at
    various shift locations, returns the position of the maximum sum,
    corresponding to the shift in the same position in shiftAmounts. If there is
    more than one absolute maximum, returns the position of the maximum
    referring to the smallest shift. If the input std::vector contains only -1 (each
    meaning the mean could not be calculated), returns -1 */
    int findIndexOfMax(const std::vector<double>& sums);

    /* Calculates the L2 index */
    double calculateL2(bool shift, int d0d1, Spline& spline_exp, Spline& spline_mod);

    /* Calculates the Pearson index */
    double calculatePearson(bool shift, int d0d1, Spline& spline_exp, Spline& spline_mod);

    /* Calculates either the norm of a spline or of its first derivative, with
    or without a shift along the x-axis */
    double calculateNorm(Spline& spline, bool shift, int d0d1);

    /* Sorts the scoresIn std::vector from highest to lowest, and sorts stdDevIn
    accordingly. Saves the results to scoresOut and stdDevOut. scoresOutNames,
    scoresOutLineTypes and scoresOutColorIndexes are also sorted as the elements
    of scoresOut */
    void sortScores(std::vector<double> scoresIn,
                    std::vector<double>& scoresOut,
                    std::vector<double> stdDevIn,
                    std::vector<double>& stdDevOut,
                    std::vector<std::string>& scoresOutNames,
                    std::vector<int>& scoresOutLineTypes,
                    std::vector<int>& scoresOutColorIndexes);

    /* Calculates the shift amounts for each model */
    void calculateShiftAmounts(Spline& spline_exp, Spline& spline_mod);

    /* Calculates the mean and the standard deviation of the mean of the
    bootstrap scores of a model, from either bootstrapScores or
    bootstrapScores_shift */
    void statisticalAnalysis(const std::vector<std::vector<double>>& bootstrap,
                             int i, /*index of the model in modelNames*/
                             double& mean,
                             double& stdDeviation);

    /* Saves the individual indexes for all the bootstrap variations to a .csv
    file */
    void saveIndexesToCsvFile();

    /* Rounds 'number' to a number of decimal places equal to decimalPlaces and
    converts it to a std::string. |Number| must be <10 */
    std::string roundToString(double number, int decimalPlaces);

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//solve(false, false, splinesExp[i][a], Sim_values[i], ObjectInput2.Exp_data[i][0], i, ObjectInput2.print_indexes, ObjectInput2.print_splines);

double Indexes::solve(bool CalculateShift,
				    bool LogScale,
                    Spline splinesExp_receive,
                    std::vector<double> valuesSim_receive,
		    std::vector<double> experimental_abscissae,
		    int a,
		    bool print_indexes,
		    bool print_splines,
            	    std::string quantity_of_interest
                    //const std::string& FolderPath,
                    //const std::string& FolderName,
                    //const std::string& FileName,
                    ) {

    Spline splinesExp;
    print_indexes = print_indexes;
    /* Splines for the models. splinesMod[i] contains the spline for
     *     modelNames[i] */
    Spline splinesMod;
    std::vector<double> Sim_values;
    
    splinesExp = splinesExp_receive;
    Sim_values = valuesSim_receive;


    if(quantity_of_interest == "IDT"){

        std::vector<double> temporary_vector;
        temporary_vector.resize(Sim_values.size());

        for (int z; z<Sim_values.size(); ++z)
        {
                temporary_vector[z] = std::log(Sim_values[z]);
        }

        splinesMod.solve(experimental_abscissae,temporary_vector,1,0,print_splines);
    } else {
        splinesMod.solve(experimental_abscissae,Sim_values,1,0,print_splines);
    }


    //for(int i=0; i<experimental_abscissae.size(); ++i)
    //	{
    //		std::cout<< "experimental_abscissae is " << experimental_abscissae[i]<<std::endl;
    //	}
    //std::cout<<"It enters the solve function in indexes class "<<std::endl;
    //folderPath = FolderPath;
    //folderName = FolderName;
    //fileName = FileName;
    calculateShift = CalculateShift;
	logScale = LogScale;

    // Processes the data in the input files
    //this->readData();

	//// Converts the ordinates to a linear scale if they were in logarithmic
	//// form
	// if (logScale == true)
	// 	for (int a=1; a<inputData.size(); ++a)
	// 		if (a%2 == 1)
	// 			for (int b=0; b<inputData[a].size(); ++b)
	// 				inputData[a][b] = exp(inputData[a][b]);

    // Creates a matrix with each line equal to a plausible version of the
    // ordinates of the experimental data, generated randomly by taking into
    // account the measurement error. The first line contains the ordinates
    // provided in the input files, and not randomly generated ordinates

    // auto bootstrapExp = std::vector<std::vector<double>>(numberOfBootstrapVariations);
    // for (int a=0; a<numberOfBootstrapVariations; ++a)
    //     // indeed bootstrap variation are initially set to be equal to the experiments
    //     bootstrapExp[a] = inputData[1];

    // create a logScale matrix with the values of the Ordinates of the experiments
	// std::vector<double> inputDataLogScale;
	// if (logScale == true) 
    // {
    //     // it saves the value of 
	// 	inputDataLogScale = inputData[1];
	// 	for (int a=0; a<inputData[1].size(); ++a)
	// 		inputData[1][a] = exp(inputData[1][a]);
	// }

    // for (int a=0; a<inputData[1].size(); ++a) {

    //     // Calculates the standard deviation for the single point
	// 	double stdDeviation = inputData[1][a]*relativeErrors[a];
    //     // create an object of class std::default_random_engine, which generates pseudo-random numbers
    //     std::default_random_engine generator;
    //     // initializes the seed of the random_engine_generator
    //     generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    //     // initializes a normal distribution for this experimental point, with mean inputData[1][a] and standard deviation stdDeviation
    //     std::normal_distribution<double> distribution(inputData[1][a],stdDeviation);
    //     // loop over the number of Bootstrap variations
    //     for (int b=1; b<numberOfBootstrapVariations; ++b) 
    //     {
    //         // generate a random number belonging to the distribution
    //         double number = distribution(generator);
    //         // if negative ordinates are not possible, then the minimum value is 0 for the random generated value.
    //         // this should only be a safety belt
    //         if (number < 0 && possibleNegativeOrdinates == false) number = 0;
    //         // if logScale is false, returns the number
    //         // otherwise, return the log of the number to bootstrap.
	// 		logScale == false ? bootstrapExp[b][a] = number :
	// 						    bootstrapExp[b][a] = log(number);

    //     }

    // }

    // if you have chosen to use the log scale you will recover the logscale
	// if (logScale == true) inputData[1] = inputDataLogScale;

    // initializes bootstrapScores and bootstrapScores_shift
    //std::cout<<"It initializes some variables we are going to need "<<std::endl;
    bootstrapScores = std::vector<std::vector<double>>(1); // there was: numberOfBootstrapVariations
    bootstrapScores_shift = std::vector<std::vector<double>>(1); // numberOfBootstrapVariations

    // Calculates the splines for the models
    // computes the number of total models, from inputData.size()
    // -1 is needed as the first two "columns" are referred to experimental data.
    numberOfModels = 1; // there was: inputData.size()/2 - 1
    // initializes an object of Splines: it will contain 1 Spline for each model
    // <Spline> splinesMod ;
    // for (int a=0; a<numberOfModels; ++a)
        // create all the models splines [TRUST THE SPLINES]
    //splinesMod[a].solve(inputData[2*(a+1)],inputData[2*(a+1)+1],1,0);

    //splinesMod_original = splinesMod;

    // create the reference experiments spline
    
    // I HAVE TO CALL THE EXP SPLINE LIKE THIS
    // I DON'T NEED TO SOLVE THE SPLINE EVERYTIME BUT I LEAVE IT HERE TO REMEMBER
    //referenceExpSpline.solve(inputData[0],inputData[1],0,0);
    // SAME HERE
    //splinesMod[a].solve(inputData[2*(a+1)],inputData[2*(a+1)+1],1,0);

    //for all BootstrapVariations
    //for (int a=numberOfBootstrapVariations-1; a>-1; --a) {
        // stores the index
        //bootstrapIndex = a;
        //assign the value of the bootstrapExp for the new evaluations
        //inputData[1] = bootstrapExp[bootstrapIndex];
        //Calculates the splines only for the experimental data
        //this->calculateSplines();
        //Initializes several elements used in the subsequent calculations
        possibleToCalculateShifts = std::vector<std::vector<bool>>(1, std::vector<bool>(numberOfModels,false));
        
        indexes = std::vector<std::vector<std::vector<double>>>(1,
                         std::vector<std::vector<double>>(numberOfModels,
                                std::vector<double>(4,-1)));

        indexes_shift = std::vector<std::vector<std::vector<double>>>(1,
                               std::vector<std::vector<double>>(numberOfModels,
                                      std::vector<double>(4,-1)));

        shifts = std::vector<std::vector<double>>(1,
                        std::vector<double>(numberOfModels,-1));

        shifts_allDifferent = std::vector<std::vector<std::vector<double>>>(1,
                                     std::vector<std::vector<double>>(numberOfModels,
                                            std::vector<double>(4,-1)));

        shiftAmounts = std::vector<std::vector<std::vector<double>>>(1,
                              std::vector<std::vector<double>>(numberOfModels));

        maxIndexes = std::vector<std::vector<int>>(1,
                            std::vector<int>(numberOfModels,-1));

        // For each knot arrangement in the experimental data, calculates the
        // indexes for each model
        //std::cout<<"It's going to calculate the indexes "<<std::endl;
        
        
        this->calculateIndexes(splinesExp, splinesMod, print_indexes);

        //std::cout<<"It's going to calculate the scores"<<std::endl;
        // Calculates the scores from the indexes
        this->calculateScores();

        // Creates the graph(s) for this versions of the experimental data
        // if (graphsBootstrap == true) {
        //     this->graphR(0,false,a+1,"Graphs Bootstrap");
        //     if (useSumOfIndexesForAlignment == true)
        //         this->graphR(0,true,a+1,"Graphs Bootstrap");
        // }

    //}

    // scores = std::vector<double>(numberOfModels,-1);
    // scores_shift = std::vector<double>(numberOfModels,-1);
    // stdDev = std::vector<double>(numberOfModels,-1);
    // stdDev_shift = std::vector<double>(numberOfModels,-1);

    // for (int a=0; a<numberOfModels; ++a) {

    //     // No shift:
    //     bool isNA = false;
    //     for (int b=0; b<numberOfBootstrapVariations; ++b)
    //         if (bootstrapScores[b][a] == -1) {
    //             isNA = true;
    //             break;
    //         }
    //     if (isNA == false)
    //         statisticalAnalysis(bootstrapScores, a, scores[a], stdDev[a]);

    //     // Shift:
    //     isNA = false;
    //     for (int b=0; b<numberOfBootstrapVariations; ++b)
    //         if (bootstrapScores_shift[b][a] == -1) {
    //             isNA = true;
    //             break;
    //         }
    //     if (isNA == false)
    //         statisticalAnalysis(bootstrapScores_shift, a,
    //                             scores_shift[a], stdDev_shift[a]);

    // }

    // sortScores(scores,
    //            scoresSorted,
    //            stdDev,
    //            stdDevSorted,
    //            modelNamesSorted,
    //            scoresSortedLineTypes,
    //            scoresSortedColorIndexes);

    // sortScores(scores_shift,
    //            scoresSorted_shift,
    //            stdDev_shift,
    //            stdDevSorted_shift,
    //            modelNamesSorted_shift,
    //            scoresSortedLineType_shift,
    //            scoresSortedColorIndexes_shift);

    // // Saves the data necessary to create the graphs with R

    // if (graphsFittingD0 == true) {
    //     this->graphsFittingR(splinesExp);
    //     this->graphsFittingR(splinesMod);
    // }

    // if (graphsD0 == true) this->graphR(0,false,1,"Graphs");
    // if (graphsD1 == true) this->graphR(1,false,1,"Graphs");

    // if (calculateShift == true && useSumOfIndexesForAlignment == true) {
    //     if (graphsD0Shift == true) this->graphR(0,true,1,"Graphs");
    //     if (graphsD1Shift == true) this->graphR(1,true,1,"Graphs");
    // }

    // // Saves the individual indexes, computed for all the bootstrap variations,
    // // in a .csv file
    // if (saveIndexesToCsv == true) this->saveIndexesToCsvFile();

    // // Saves the data necessary to plot the splines to .txt files.
    // // Used to interface with Alice's work
    // if (saveGraphData == true) this->saveGraphDataToFiles();

    return FINAL_SCORE;
}



void Indexes::readData() {

    // Copies the data from folderPath/name_exp.txt to datastringExp

    std::vector<std::vector<std::string>> datastringExp;

    std::string fileNameWithPathExp = folderPath + fileName + "_exp.txt";
    const char* fileNameWithPathCharExp = fileNameWithPathExp.c_str();
    std::string lineExp;
    std::string elementExp;
    std::vector<std::string> emptyvectorExp;
    int rowIndexExp = 0;
    // opens a buffer where to store that file containing experiments
    std::ifstream myfileExp(fileNameWithPathCharExp);
    // while it is still possible to get a line from the buffer and assign it to lineExp
    while (getline(myfileExp,lineExp))
    {
        // pushing back an empty std::vector to first empty element 
        datastringExp.push_back(emptyvectorExp);
        // assign the c-like std::string lineExp to lineStream of type std::std::stringstream, which wrap the std::string with a stream interface, so we can use << and >> on it
        std::stringstream lineStream(lineExp);
        // this while stores all the file within a matrix, called datastringExp
        while (lineStream >> elementExp)
            datastringExp[rowIndexExp].push_back(elementExp);
        ++rowIndexExp;
    }
    //close the buffer
    myfileExp.close();

    // Converts the input data from std::string to double and saves it to
    // dataDoubleExp, with the data saved in reverse order relative to the order
    // in the input file (to speed up processing)

    auto dataDoubleExp = std::vector<std::vector<double>>(3);

    for (int a=datastringExp.size()-1; a>0; --a)
        if (datastringExp[a][0] != "NA") {
            dataDoubleExp[0].push_back(stod(datastringExp[a][0]));
            dataDoubleExp[1].push_back(stod(datastringExp[a][1]));
            // it is interesting that you can even not specify the error, there is a default one
            if (datastringExp[0].size() == 3)
                dataDoubleExp[2].push_back(stod(datastringExp[a][2]));
            else
                dataDoubleExp[2].push_back(defaultRelativeError);
        }

    // Fills dataSortedExp, corresponding to dataDoubleExp with the data sorted
    // in ascending order with respect to the abscissae
    auto dataSortedExp = std::vector<std::vector<double>>(3);

    // until there are data within the matrix of double
    while (dataDoubleExp[0].size() > 0)
        // so it goes from the final position up, one by one. 
        for (int a=dataDoubleExp[0].size()-1; a>-1; --a) 
        {
            // assign value of abscissa
            double value = dataDoubleExp[0][a];
            // assign index
            int index = a;
            // compare value with all the others above, if any abscissa is lower than value, this becomes the new value
            // and its index becomes the new index
            for (int b=dataDoubleExp[0].size()-1; b>-1; --b)
                if (dataDoubleExp[0][b] < value) 
                {
                    value = dataDoubleExp[0][b];
                    index = b;
                }
            // the minimum value become the first in sorted matrix, together with exp and error
            dataSortedExp[0].push_back(dataDoubleExp[0][index]);
            dataSortedExp[1].push_back(dataDoubleExp[1][index]);
            dataSortedExp[2].push_back(dataDoubleExp[2][index]);
            // the same value is erased from Double Exp
            dataDoubleExp[0].erase(dataDoubleExp[0].begin()+index);
            dataDoubleExp[1].erase(dataDoubleExp[1].begin()+index);
            dataDoubleExp[2].erase(dataDoubleExp[2].begin()+index);
        }

    // Copies the data from folderPath/name_mod.txt to datastd::stringMod
    std::vector<std::vector<std::string>> datastringMod;

    //aquire path +name
    std::string fileNameWithPathMod = folderPath + fileName + "_mod.txt";
    // convert to c-like std::string
    const char* fileNameWithPathCharMod = fileNameWithPathMod.c_str();
    // initializes to element std::vectors of std::strings
    std::string lineMod;
    std::string elementMod;
    // initialize an empty std::vector to be pushed back into the datastd::stringMod
    std::vector<std::string> emptyvectorMod;
    int rowIndexMod = 0;
    // open a myfileMod buffer were to store the content of the file
    std::ifstream myfileMod(fileNameWithPathCharMod);
    // until it is possible to get lines from the buffer
    while (getline(myfileMod,lineMod)) 
    {
        // push back an empty std::vector of std::strings into the structure
        datastringMod.push_back(emptyvectorMod);
        // 
        std::stringstream lineStream(lineMod);
        // assign the c-like std::string lineExp to lineStream of type std::std::stringstream, which wrap the std::string with a stream interface, so we can use << and >> on it
        // and actually populate each line of datastd::stringMod with std::strings
        while (lineStream >> elementMod)
            datastringMod[rowIndexMod].push_back(elementMod);
        ++rowIndexMod;
    }
    myfileMod.close();

    // Saves the model names to modelNames
    for (int a=0; a<datastringMod[0].size(); ++a)
        // (a%2 != 0) this is true only for 1,3,5,7,9,11 etc, so basically for 38_1, 38_2, 38_3 (look at paper)
        if (a%2 != 0) 
        {
            // take the name of the model into modelName variable
            std::string modelName = datastringMod[0][a];
            // erase the part 38_ from name and only the number remains
            modelName.erase(0,fileName.size()+1);
            // add it into something else
            modelNames.push_back(modelName);
        }
    
    // Converts the input data from std::string to double and saves it to
    // dataDoubleMod, with the data saved in reverse order relative to the order
    // in the input file (to speed up processing)

    auto dataDoubleMod = std::vector<std::vector<double>>(datastringMod[0].size());
    // looping over the size of the first row of the file
    for (int a=0; a<datastringMod[0].size(); ++a)
        // this time only for the values of the abscissa 
        if (a%2 == 0)
            //looping over the number of std::strings in ascending way
            for (int b=datastringMod.size()-1; b>0; --b)
                // only if the std::string is a number and not a NA (not a number)
                if (datastringMod[b][a] != "NA")
                {
                    // push back the double related to abscissa and sim value
                    dataDoubleMod[a].push_back(stod(datastringMod[b][a]));
                    dataDoubleMod[a+1].push_back(stod(datastringMod[b][a+1]));
                }

    // Fills dataSortedExp, corresponding to dataDoubleExp with the data sorted
    // in ascending order with respect to the abscissae

    // initialize an automatic type assigned matrix of double, which has the size
    // of the row of the _mod.txt file
    auto dataSortedMod = std::vector<std::vector<double>>(dataDoubleMod.size());
    
    for (int a=0; a<dataDoubleMod.size(); ++a)
        // only for abscissa columns
        if (a%2 == 0)
            // until the columns of the original matrix is full, keep going
            while (dataDoubleMod[a].size() > 0)
                // looping over the columns from bottom to top
                for (int b=dataDoubleMod[a].size()-1; b>-1; --b) 
                {
                    // takes the value from the matrix of doubles
                    double value = dataDoubleMod[a][b];
                    // create an index
                    int index = b;
                    // loop over the entire column searching for values, which are lower than
                    // the one assigned to value variable.
                    for (int c=dataDoubleMod[a].size()-1; c>-1; --c)
                        // in the case it finds a value which is lower than the bottom one, 
                        // value is replaced, and compared with the ones above. This way he
                        // detect the minimum value and send it to dataSortedMod.
                        if (dataDoubleMod[a][c] < value) 
                        {
                            value = dataDoubleMod[a][c];
                            index = c;
                        }
                    dataSortedMod[a].push_back(dataDoubleMod[a][index]);
                    dataSortedMod[a+1].push_back(dataDoubleMod[a+1][index]);
                    // contemporarily the added values are deleted from the matrix of non-ordered
                    // doubles
                    dataDoubleMod[a].erase(dataDoubleMod[a].begin()+index);
                    dataDoubleMod[a+1].erase(dataDoubleMod[a+1].begin()+index);
                }

    // If two or more abscissae are the same, replaces them with a single
    // abscissa with the mean of the corresponding ordinates as ordinate, and
    // (for experimental data) with the mean of the corresponding relative
    // errors as relative error. Copies the resulting abscissae and ordinates to
    // allInputData and the resulting relative errors to allRelativeErrors

    // initializes allInputData with the size of the _mod.txt matrix + 2 probably for the experiments
    auto allInputData = std::vector<std::vector<double>>(2+dataSortedMod.size());
    std::vector<double> allRelativeErrors;
    // push back the first line of exp abscissa, experiments and relative uncertainty
    allInputData[0].push_back(dataSortedExp[0][0]);
    allInputData[1].push_back(dataSortedExp[1][0]);
    allRelativeErrors.push_back(dataSortedExp[2][0]);

    int nExp = 1;
    // for the number of rows in the experimental data file
    for (int a=1; a<dataSortedExp[0].size(); ++a) 
    {
        // if the value in dataSortedExp [0][a] is not the same of the one
        // already present in the end of allInputData[0], then push back.
        if (dataSortedExp[0][a] != allInputData[0].back()) {
            allInputData[0].push_back(dataSortedExp[0][a]);
            allInputData[1].push_back(dataSortedExp[1][a]);
            allRelativeErrors.push_back(dataSortedExp[2][a]);
            nExp = 1;
        }
        // if they are the same value then do this routine
        else {
            // the two experiments are summed up and the errors as well
            allInputData[1].back() =
                allInputData[1].back()*(double)nExp + dataSortedExp[1][a];
            allRelativeErrors.back() =
                allRelativeErrors.back()*(double)nExp + dataSortedExp[2][a];
            ++nExp;
            allInputData[1].back() /= (double)nExp;
            allRelativeErrors.back() /= (double)nExp;
        }
    }

    // it does the same it did for experiments, also for simulations
    for (int a=0; a<dataSortedMod.size(); ++a)
        
        if (a%2 == 0)
            if (dataSortedMod[a].size() > 1) {
                allInputData[a+2].push_back(dataSortedMod[a][0]);
                allInputData[a+3].push_back(dataSortedMod[a+1][0]);
            }

    int nMod = 1;

    for (int a=0; a<dataSortedMod.size(); ++a)
        if (a%2 == 0)
            if (dataSortedMod[a].size() > 1)
                for (int b=1; b<dataSortedMod[a].size(); ++b) {
                    if (dataSortedMod[a][b] != allInputData[a+2].back()) {
                        allInputData[a+2].push_back(dataSortedMod[a][b]);
                        allInputData[a+3].push_back(dataSortedMod[a+1][b]);
                        nMod = 1;
                    }
                    else {
                        allInputData[a+3].back() =
                            allInputData[a+3].back()*(double)nMod +
                            dataSortedMod[a+1][b];
                        ++nMod;
                        allInputData[a+3].back() /= (double)nMod;
                    }
                }

    // Obtains inputData and relativeErrors by removing excess points with
    // ordinate 0 on the left and on the right of experimental data and models

    // INITIALIZE InputData, with the first direction of the same size of allInputData
    inputData = std::vector<std::vector<double>>(allInputData.size());

    for (int a=0; a<allInputData.size(); ++a)
        if (allInputData[a].size() > 1) {

            if (a == 0) {
                
                // initialize and gives the max/mix Ordinate of the experiments randomly
                // the same value as the double stored in the first position
                double maxOrdinate = allInputData[a+1][0];
                double minOrdinate = allInputData[a+1][0];

                // search for actual min/max Ordinate of the experiments 
                for (int b=1; b<allInputData[a].size(); ++b) {
                    if (allInputData[a+1][b] > maxOrdinate)
                        maxOrdinate = allInputData[a+1][b];
                    if (allInputData[a+1][b] < minOrdinate)
                        minOrdinate = allInputData[a+1][b];
                }
                // defines the height as the difference between the max/min ordinates 
                double height = maxOrdinate - minOrdinate;
                // compute the minHeight thorugh the use of a variable defined and assigned with a value in GlobalVariables.h
                double minHeight =
                    fractionOfOrdinateRangeForAsymptoteIdentification * height;

                // finds the first point that differs of more than minHeight from the previous, in terms of Ordinates
                int indexOfFirstPoint = 0;
                for (int b=1; b<allInputData[a].size(); ++b)
                    if (fabs(allInputData[a+1][b]-allInputData[a+1][b-1]) >
                        minHeight) {
                        indexOfFirstPoint = b-1;
                        break;
                    }

                // it applies the same concept to find last point!
                int indexOfLastPoint = allInputData[a].size()-1;
                for (int b=allInputData[a].size()-2; b>-1; --b)
                    if (fabs(allInputData[a+1][b]-allInputData[a+1][b+1]) >
                        minHeight) {
                        indexOfLastPoint = b+1;
                        break;
                    }
                // once it has initialized the max/min points indexes and values
                // he initializes std::vectors in the positions 0 and 1 to have the size of the difference between the these two indexes
                inputData[a] =
                    std::vector<double>(indexOfLastPoint-indexOfFirstPoint+1);
                inputData[a+1] =
                    std::vector<double>(indexOfLastPoint-indexOfFirstPoint+1);

                // populates these first to positions of inputData
                for (int b=0; b<inputData[a].size(); ++b)
                    inputData[a][b] = allInputData[a][indexOfFirstPoint+b];

                for (int b=0; b<inputData[a].size(); ++b)
                    inputData[a+1][b] = allInputData[a+1][indexOfFirstPoint+b];
                
                // populates also relative errors std::vector with the same size
                for (int b=indexOfFirstPoint; b<=indexOfLastPoint; ++b)
                    relativeErrors.push_back(allRelativeErrors[b]);

            }
            // for the values of the simulations instead, it takes all of them here. And it makes a lot of sense actually.
            if (a > 0 && a%2 == 0) {

                inputData[a]   = allInputData[a];
                inputData[a+1] = allInputData[a+1];

            }

        }

}



// void Indexes::calculateSplines() {

    
//     if (useIndexesToChooseExpNodes == false)
//         splinesExp = std::vector<Spline>(1);
//     else {
//         if (inputData[0].size() < 3)
//             splinesExp = std::vector<Spline>(1);
//         else if (inputData[0].size() < 5)
//             splinesExp = std::vector<Spline>(2);
//         else
//             splinesExp = std::vector<Spline>(3);
//     }
    
//     splinesExp[0].solve(inputData[0],inputData[1],0,0);
//     splinesExp[0].removeAsymptotes();

//     if (splinesExp.size() > 1) {
//         splinesExp[1].solve(inputData[0],inputData[1],0,2);
//         splinesExp[1].removeAsymptotes();
//     }

//     if (splinesExp.size() > 2) {
//         splinesExp[2].solve(inputData[0],inputData[1],0,5);
//         splinesExp[2].removeAsymptotes();
//     }

// }



void Indexes::calculateIndexes(Spline& splinesExp, Spline& splinesMod, bool print_indexes) {
    
    //std::cout<< "- enters the calculated indexes -  " << std::endl;
    // it will always be 
    // 
    //splinesMod = splinesMod_original;

    // Checks whether any model can be considered a flat line with ordinate = 0
    // when compared to the experimental data
    if (possibleNegativeOrdinates == false) {
        double minHeight =
            splinesExp.yD0Max * fractionOfExpHeightForZeroLineIdentification;
        

    if (splinesMod.possibleToCalculateSpline == true)
        if (splinesMod.yD0Max < minHeight)
                splinesMod.possibleToCalculateSpline = false;
    }
    //std::cout<<"It is possible to calculate this spline?  " <<  splinesMod.possibleToCalculateSpline<<std::endl;
    //std::cout<<"The fractionOfExpRangeForModelExtrapolation is  " <<  fractionOfExpRangeForModelExtrapolation <<std::endl;
    //std::cout<<"The splinesExp.xRange is " <<  splinesExp.xRange  <<std::endl;
    // selects just t
    double minLengthInCommon =
        splinesExp.xRange * fractionOfExpRangeForModelExtrapolation;

    //std::cout<<"minLength in common is "<< minLengthInCommon <<std::endl;

    maxExtrapolatedLength = splinesExp.xRange - minLengthInCommon;

    //std::cout<<"maxExtrapolated length is "<< maxExtrapolatedLength <<std::endl;

    // Adds segments to the left and to the right of the models so that their
    // ordinates span the entirety of the experimental data range, regardless of
    // the shift value

    if (splinesMod.possibleToCalculateSpline == true)
        splinesMod.extendSpline(maxExtrapolatedLength,
                                maxExtrapolatedLength);

    //std::cout<<"Spline has been extended"<<std::endl;
    // Normalizes the spline coefficients of the experimental data and of the
    // first derivative of the experimental data
    
    splinesExp.normalizeCoefficients(
                                    -splinesExp.yD0Min,
                                    splinesExp.yD0Max-splinesExp.yD0Min,
                                    splinesExp.yD1MaxAbs);

    //std::cout<<"coefficients have been normalized"<<std::endl;
    
    // Calculates the values of shiftAmounts
    //std::cout<<"Calculate shift? "<< calculateShift << "should be false"<<std::endl;

    if (calculateShift == true) {

        shiftAroundMaximum =
            splinesExp.xRange * fractionOfExpRangeForShiftAroundMaximum;

        distanceBetweenShiftedPoints =
            splinesExp.xRange * fractionOfExpRangeForMinShift;

        
        if (splinesMod.possibleToCalculateSpline == true)
            calculateShiftAmounts(splinesExp, splinesMod);

    }

    
    // Initializes several elements necessary for the calculation of the indexes
    double height = (splinesExp.knots.back() - splinesExp.knots[0]) /
                    (double)numberOfTrapezoids;

    halfHeight = height/2.;
    //std::cout<<"Height is "<< height <<std::endl;
    // Initialises the abscissa with 100 positions
    std::vector<double> x;
    for (int a=0; a<numberOfTrapezoids; ++a)
        x.push_back(splinesExp.knots[0]+(double)a*height);
    x.push_back(splinesExp.knots.back());

    // creates this 100 positions vector where we have 1 vector for each position
    // 1 vector of 4 elements, each element is Abscissa to the power of 0,1,2,3.
    xPowers = std::vector<std::vector<double>>(numberOfTrapezoids+1,std::vector<double>(m,1));
    for (int a=0; a<numberOfTrapezoids+1; ++a)
        for (int b=1; b<m; ++b)
            xPowers[a][b] = xPowers[a][b-1] * x[a];
    // Here it computes the corresponding 100 ordinates to x, computed with normalised coefficients
    yD0 = std::vector<double>(numberOfTrapezoids+1);
    for (int a=0; a<numberOfTrapezoids+1; ++a)
        yD0[a] = splinesExp.D0(xPowers[a]);
    //for (int a=0; a<numberOfTrapezoids+1; ++a)
    //	std::cout<<"yD0["<<a<<"]"<<" is "<< yD0[a] << std::endl;
    // Here it computes the corresponding 100 first derivative corresponding to x, computed with normalised coefficients
    yD1 = std::vector<double>(numberOfTrapezoids+1);
    for (int a=0; a<numberOfTrapezoids+1; ++a)
        yD1[a] = splinesExp.D1(xPowers[a]);
    //for (int a=0; a<numberOfTrapezoids+1; ++a)
    //    std::cout<<"yD1["<<a<<"]"<<" is "<< yD1[a] << std::endl;
    
    //std::cout<<"It has prepared x, xPowers, yD0, yD1" <<std::endl;
    // It stores in xL2Means[i] the average between x[i] and x[i+1] !!!
    std::vector<double> xL2Means;
    for (int a=0; a<numberOfTrapezoids; ++a)
        xL2Means.push_back((x[a]+x[a+1])/2.);

    //std::cout<<"It has prepared xL2Means" <<std::endl;
    // it creates powers vectors for these average absissas vector!!!
    auto xL2MeansPowers =
        std::vector<std::vector<double>>(numberOfTrapezoids,std::vector<double>(m,1));
    for (int a=0; a<numberOfTrapezoids; ++a)
        for (int	 b=1; b<m; ++b)
            xL2MeansPowers[a][b] = xL2MeansPowers[a][b-1] * xL2Means[a];
    //std::cout<<"It has prepared xL2MeansPowers" <<std::endl;
    // Calculates the norms of the experimental data
    normD0Exp = calculateNorm(splinesExp,false,0);
    normD1Exp = calculateNorm(splinesExp,false,1);
    
    //std::cout<<"normD0Exp is "<< normD0Exp <<std::endl;
    //std::cout<<"normD1Exp is "<< normD1Exp <<std::endl;
    // Calculates the indexes for each model, if possible
    //std::cout<<"splinesMod.possibleToCalculateSpline is "<< splinesMod.possibleToCalculateSpline <<std::endl;
    

        //indexOfModelBeingCalculated = a;

        // Normalizes the coefficients of the spline and of the first
        // derivative of the spline
        splinesMod.normalizeCoefficients(
                                            -splinesExp.yD0Min,
                                            splinesExp.yD0Max-splinesExp.yD0Min,
                                            splinesExp.yD1MaxAbs);

            // Checks whether the experimental data and the model have enough
            // non-extrapolated lenght in common
            	bool possibleToCalculateIndexes = true;
            	double limitLeft = splinesMod.originalAbscissae[0];
	    	if (limitLeft < splinesExp.xMinNoAsymptotes)
		{
                	limitLeft = splinesExp.xMinNoAsymptotes;
		}
            	double limitRight = splinesMod.originalAbscissae.back();
	    
        	if (limitRight > splinesExp.xMaxNoAsymptotes)
		{
                	limitRight = splinesExp.xMaxNoAsymptotes;
		}
        	if (limitRight - limitLeft < minLengthInCommon)
		{
                // here i think i have to put a penalty function if it's not possible to
                // compute the indexes
                possibleToCalculateIndexes = false;
		}
	    //std::cout << "is it possible to calculate the indexes? " << possibleToCalculateIndexes << std::endl;
            if (possibleToCalculateIndexes == true) 
            {
                // Calculates the norms for the model
                normD0Mod = calculateNorm(splinesMod,false,0);
                normD1Mod = calculateNorm(splinesMod,false,1);
                // Calculates d0L2, d1L2, d0Pe, d1Pe (no shift)
                indexes[0][0][0] = calculateL2(false,0, splinesExp, splinesMod);
                indexes[0][0][1] = calculateL2(false,1, splinesExp, splinesMod);
                indexes[0][0][2] = calculatePearson(false,0,splinesExp, splinesMod);
                indexes[0][0][3] = calculatePearson(false,1,splinesExp, splinesMod);
            	
		if(print_indexes)
		{
			std::cout<<"indexes[0][0][0] "<< indexes[0][0][0]  <<std::endl;
			std::cout<<"indexes[0][0][1] "<< indexes[0][0][1]  <<std::endl;
			std::cout<<"indexes[0][0][2] "<< indexes[0][0][2]  <<std::endl;
			std::cout<<"indexes[0][0][3] "<< indexes[0][0][3]  <<std::endl;
	    	}	
	     }

            //std::cout<<"It has calculated L2 norms and Pearsons" <<std::endl;

            if (calculateShift == true &&
                possibleToCalculateShifts[0][0] == true) {

                // indexValues will contain the values for an index calculated
                // for each value of shiftAmounts[i][a].
                // indexValues[i][0] = d0L2, indexValues[i][1] = d1L2,
                // indexValues[i][2] = d0Pe, indexValues[i][3] = d1Pe
                auto indexValues =
                    std::vector<std::vector<double>>(4,
                           std::vector<double>(shiftAmounts[0][0].size(),-1));

                auto sumsOfIndexValues =
                    std::vector<double>(shiftAmounts[0][0].size(),-1);

                // Calculates d0L2, d1L2, d0Pe and d1Pe for the various shift
                // amounts
                for (int b=0; b<shiftAmounts[0][0].size(); ++b) {

                    splinesMod.calculateShift(shiftAmounts[0][0][b]);

                    // Calculates the norms for the shifted model
                    normD0Mod = calculateNorm(splinesMod,true,0);
                    normD1Mod = calculateNorm(splinesMod,true,1);

                    // Calculates d0L2, d1L2, d0Pe, d1Pe (with shift)
                    indexValues[0][b] = calculateL2(true,0, splinesExp, splinesMod);
                    indexValues[1][b] = calculateL2(true,1, splinesExp, splinesMod);
                    indexValues[2][b] = calculatePearson(true,0, splinesExp, splinesMod);
                    indexValues[3][b] = calculatePearson(true,1, splinesExp, splinesMod);

                } // End of the for cycle for the shifts

                for (int b=0; b<shiftAmounts[0][0].size(); ++b) {

                    bool possibleToCalculateSum = true;

                    for (int c=0; c<4; ++c)
                        if (indexValues[c][b] == -1)
                            possibleToCalculateSum = false;

                    if (possibleToCalculateSum == true)
                        sumsOfIndexValues[b] =
                            indexValues[0][b] + indexValues[1][b] +
                            indexValues[2][b] + indexValues[3][b];

                }

                int indexOfMax = findIndexOfMax(sumsOfIndexValues);

                maxIndexes[0][0] = indexOfMax;

                if (useSumOfIndexesForAlignment == true)
                    if (indexOfMax != -1) {

                        indexes_shift[0][0][0] = indexValues[0][indexOfMax];
                        indexes_shift[0][0][1] = indexValues[1][indexOfMax];
                        indexes_shift[0][0][2] = indexValues[2][indexOfMax];
                        indexes_shift[0][0][3] = indexValues[3][indexOfMax];

                        double ratio = abs(shiftAmounts[0][0][indexOfMax])/
                                       splinesExp.xRange;

                        if (ratio > 1.) ratio = 1.;

                        shifts[0][0] = 1.-ratio;

                    }

                if (useSumOfIndexesForAlignment == false)
                    for (int b=0; b<4; ++b) {
                        int c = findIndexOfMax(indexValues[b]);
                        if (c != -1) {
                            indexes_shift[0][0][b] = indexValues[b][c];

                            double ratio = abs(shiftAmounts[0][0][c])/
                                           splinesExp.xRange;

                            if (ratio > 1.) ratio = 1.;

                            shifts_allDifferent[0][0][b] = 1.-ratio;
                        }
                    }

            } // End of if(calculateShift == true &&
              //           possibleToCalculateShifts[i][a] == true)

         // End of if (splinesMod[a].possibleToCalculateSpline == true)
     // End of the for cycle for each model

}



void Indexes::calculateScores() {

    auto matrix_scores = std::vector<std::vector<double>>(1,
                                std::vector<double>(numberOfModels,-1));
    auto matrix_scores_shift = std::vector<std::vector<double>>(1,
                                      std::vector<double>(numberOfModels,-1));

    // Calculates matrix_scores
    for (int i=0; i<1; ++i)
        for (int a=0; a<numberOfModels; ++a) {
            
            // initializes a boolean to check if the indexes d0, d1, p0, p1 can be used for scores estimatation.
            bool scoreIsNA = false;
            // check all of the indexes are -1, meaning it was not possible to get the index for some reason.
            for (int b=0; b<4; ++b)
                if (indexes[i][a][b] == -1)
                    scoreIsNA = true;
            // check all of the indexes are not numbers, in this case, it won't be possible to compute the CM index.
            for (int b=0; b<4; ++b)
                if (isnan(indexes[i][a][b]) == true)
                    scoreIsNA = true;

            if (scoreIsNA == false)
                matrix_scores[i][a] = (indexes[i][a][0]+indexes[i][a][1]+
                                       indexes[i][a][2]+indexes[i][a][3]) / 4.;

        }

    FINAL_SCORE = matrix_scores[0][0];
    //std::cout<<"The CM index is equal to "<< FINAL_SCORE << std::endl;
    // AB // I COULD EVEN STOP HERE // 

    // Calculates matrix_scores_shift, without considering the shift, in order
    // to select the best spline from splinesExp
    if (calculateShift == true)
        for (int i=0; i<1; ++i)
            for (int a=0; a<numberOfModels; ++a) {

                bool scoreIsNA = false;

                for (int b=0; b<4; ++b)
                    if (indexes_shift[i][a][b] == -1)
                        scoreIsNA = true;

                for (int b=0; b<4; ++b)
                    if (isnan(indexes_shift[i][a][b]) == true)
                        scoreIsNA = true;

                if (scoreIsNA == false)
                    matrix_scores_shift[i][a] =
                        (indexes_shift[i][a][0]+indexes_shift[i][a][1]+
                         indexes_shift[i][a][2]+indexes_shift[i][a][3]) / 4.;

            }

    // Selects the best spline of splinesExp
    auto meansOfScores = std::vector<double>(1,-1);
    // I don't think i need this to be honest, what do you need the mean of the scores between different models for?
    // if (calculateShift == false)
    // {
    //     for (int i=0; i<splinesExp.size(); ++i) 
    //     {
            
    //         double sumOfScores = 0;
    //         int counter = 0;
    //         for (int a=0; a<numberOfModels; ++a)
    //             if (matrix_scores[i][a] != -1) {
    //                 sumOfScores += matrix_scores[i][a];
    //                 ++counter;
    //             }
    //         if (counter != 0)
    //             meansOfScores[i] = sumOfScores / (double)counter;
    //     }
    // }
    // if (calculateShift == true)
    //     for (int i=0; i<splinesExp.size(); ++i) {
    //         double sumOfScores = 0;
    //         int counter = 0;
    //         for (int a=0; a<numberOfModels; ++a)
    //             if (matrix_scores_shift[i][a] != -1) {
    //                 sumOfScores += matrix_scores_shift[i][a];
    //                 ++counter;
    //             }
    //         if (counter != 0)
    //             meansOfScores[i] = sumOfScores / (double)counter;
    //     }
    // Even here, why the hell do we even have this in our code, jesus.
    // indexBestSplineExp = 0;
    // double maxMeanOfScores = -1.;
    // for (int i=0; i<splinesExp.size(); ++i)
    //     if (meansOfScores[i] != -1)
    //         if (meansOfScores[i] > maxMeanOfScores) {
    //             maxMeanOfScores = meansOfScores[i];
    //             indexBestSplineExp = i;
    //         }
    // rather, i would need this, 
    bootstrapScores[bootstrapIndex] = matrix_scores[indexBestSplineExp];

    bootstrapScores_shift[bootstrapIndex] =
        matrix_scores_shift[indexBestSplineExp];
    if (calculateShift == true)
        for (int a = 0; a<numberOfModels; ++a)
            if (matrix_scores_shift[indexBestSplineExp][a] != -1)
                useSumOfIndexesForAlignment == true ?

                    bootstrapScores_shift[bootstrapIndex][a] =
                        (matrix_scores_shift[indexBestSplineExp][a]*4.+
                        (shifts[indexBestSplineExp][a]*2.))/6. :

                    bootstrapScores_shift[bootstrapIndex][a] =
                        (matrix_scores_shift[indexBestSplineExp][a]*4.+
                         (shifts_allDifferent[indexBestSplineExp][a][0]+
                          shifts_allDifferent[indexBestSplineExp][a][1]+
                          shifts_allDifferent[indexBestSplineExp][a][2]+
                          shifts_allDifferent[indexBestSplineExp][a][3])/2.)/6.;

    // Rows: models, Columns: d0L2, d1L2, d0Pe, d1Pe, shift (or shifts if all
    // different)
    auto idx = std::vector<std::vector<double>>(numberOfModels);

    for (int a=0; a<numberOfModels; ++a) {

        if (calculateShift == true)
            idx[a] = indexes_shift[0][a];
        else
            idx[a] = indexes[0][a];

        // if (useSumOfIndexesForAlignment == true)
        //     idx[a].push_back(shifts[indexBestSplineExp][a]);
        // else
        //     for (int b=0; b<4; ++b)
        //         idx[a].push_back(shifts_allDifferent[indexBestSplineExp][a][b]);
    }

    allIndexes.push_back(idx);

}



// void Indexes::graphsFittingR(std::vector<Spline>& splines) {

//     for (int a=0; a<splines.size(); ++a)
//         if (splines[a].possibleToCalculateSpline == true) {

//             bool skipSpline = false;
//             if (splines[0].splineType == 0 /*Experimental data*/)
//                 if (a != indexBestSplineExp)
//                     skipSpline = true;

//             if (skipSpline == false) {

//                 std::string whichSpline;
//                 if (splines[0].splineType == 0 /*Experimental data*/)
//                     whichSpline = "Exp";
//                 if (splines[0].splineType == 1 /*Model*/)
//                     whichSpline = modelNames[a];

//                 std::string scriptNameAndPath =
//                     "./Scripts/"+folderName+"_"+fileName+".R";
//                 const char* scriptNameAndPathChar = scriptNameAndPath.c_str();
//                 std::ofstream script;
//                 script.open(scriptNameAndPathChar,std::ios::app);

//                 double xMin = splines[a].knots[0];
//                 double xMax = splines[a].knots.back();

//                 double yMin = splines[a].D0(splines[a].knots[0]);
//                 double yMax = splines[a].D0(splines[a].knots[0]);

//                 double distance = (xMax - xMin) / (double)graphPoints;

//                 script << "x = c(";
//                 for (int b=0; b<graphPoints; ++b)
//                     script << xMin+(double)b*distance << ", ";
//                 script << xMax << ")\n";

//                 double y;
//                 script << "y = c(";
//                 for (int b=0; b<graphPoints; ++b) {
//                     double t = xMin+(double)b*distance;
//                     y = splines[a].D0(t);
//                     if (y < yMin) yMin = y;
//                     if (y > yMax) yMax = y;
//                     script << y << ", ";
//                 }
//                 y = splines[a].D0(xMax);
//                 if (y < yMin) yMin = y;
//                 if (y > yMax) yMax = y;
//                 script << y << ")\n";

//                 for (int b=0; b<splines[a].originalOrdinates.size(); ++b) {
//                     if (splines[a].originalOrdinates[b] < yMin)
//                         yMin = splines[a].originalOrdinates[b];
//                     if (splines[a].originalOrdinates[b] > yMax)
//                         yMax = splines[a].originalOrdinates[b];
//                 }

//                 script << "xExp = c(";
//                 for (int b=0; b<splines[a].abscissae.size()-1; ++b)
//                     script << splines[a].abscissae[b] << ", ";
//                 script << splines[a].abscissae.back() << ")\n";

//                 script << "yExp = c(";
//                 for (int b=0; b<splines[a].ordinates.size()-1; ++b)
//                     script << splines[a].ordinates[b] << ", ";
//                 script << splines[a].ordinates.back() << ")\n";

//                 double quality = (double)graphQuality / 1000.;

//                 script << "png('./Graphs Fitting/" << folderName << "_"
//                     << fileName << "_" << whichSpline << ".png', width="
//                     << graphQuality << ", height=" << graphQuality
//                     << ", type='cairo')\n"
//                     << "par(mar=c(" << 4.*quality << "," << 4.*quality << ","
//                     << 4.*quality << "," << 2.*quality << "))\n"
//                     << "plot(xExp, yExp, xlim=c("
//                     << splines[a].originalAbscissae[0] << ","
//                     << splines[a].originalAbscissae.back() << "), ylim=c("
//                     << yMin << "," << yMax << "), main='" << folderName << "  "
//                     << fileName << "  " << whichSpline
//                     << "', pch='o', col='dimgray', xlab='', ylab='', cex="
//                     << quality << ", cex.axis=" << 1.5*quality << ", cex.lab="
//                     << 1.5*quality << ", cex.main=" << 1.6*quality << ", mgp=c("
//                     << 0 << "," << 2.*quality << "," << 0 << "))\n"
//                     << "box(lwd=" << quality << ")\n"
//                     << "axis(side=1, labels=FALSE, lwd.ticks=" << quality
//                     << ", tcl=" << -0.75*quality << ")\n"
//                     << "axis(side=2, labels=FALSE, lwd.ticks=" << quality
//                     << ", tcl=" << -0.75*quality << ")\n"
//                     << "lines(x, y, col='blue', type='l', lwd=" << 3.*quality
//                     << ")\n"
//                     << "abline(v=c(";
//                 for (int b=0; b<splines[a].numberOfKnots; ++b) {
//                     script << std::to_string(splines[a].knots[b]);
//                     if (b != splines[a].numberOfPolynomials)
//                         script << ", ";
//                     else
//                         script << ")";
//                 }
//                 script << ", col=\'dimgray\', lty=\'longdash\', lwd="
//                        << 1.*quality << ")\ndev.off()\n\n";

//                 script.close();

//             } // end of if (skipSpline == false)

//         }

// }



// void Indexes::graphR(int derivativeOrder,
//                      bool shift,
//                      int bootstrapIndex,
//                      const std::string& folder) {

//     // Saves the std::vectors with abscissae and ordinates to a .R file

//     int j = indexBestSplineExp;

//     bool bootstrap = false;
//     if (folder == "Graphs Bootstrap") bootstrap = true;

//     double xMin = splinesExp[j].knots[0];
//     for (int a=0; a<numberOfModels; ++a)
//         if (splinesMod[a].possibleToCalculateSpline == true)
//             if (splinesMod[a].originalAbscissae[0] < xMin)
//                 xMin = splinesMod[a].originalAbscissae[0];

//     double xMax = splinesExp[j].knots.back();
//     for (int a=0; a<numberOfModels; ++a)
//         if (splinesMod[a].possibleToCalculateSpline == true)
//             if (splinesMod[a].originalAbscissae.back() > xMax)
//                 xMax = splinesMod[a].originalAbscissae.back();

//     std::string scriptNameAndPath = "./Scripts/"+folderName+"_"+fileName+".R";
//     const char* scriptNameAndPathChar = scriptNameAndPath.c_str();
//     std::ofstream script;
//     script.open(scriptNameAndPathChar,std::ios::app);

//     double yMin;
//     double yMax;

//     if (derivativeOrder == 0) {

//         if (bootstrap == false) {
//             yMin = splinesExp[j].yD0Min;
//             yMax = splinesExp[j].yD0Max;
//         }
//         else {
//             double max = 0;
//             for (int a=0; a<relativeErrors.size(); ++a)
//                 if (relativeErrors[a] > max) max = relativeErrors[a];
//             yMin = referenceExpSpline.yD0Min * (1.-2*max);
//             yMax = referenceExpSpline.yD0Max * (1.+2*max);
//         }

//         for (int a=0; a<numberOfModels; ++a)
//             if (splinesMod[a].possibleToCalculateSpline == true)
//                 if (splinesMod[a].yD0MinOriginal < yMin)
//                     yMin = splinesMod[a].yD0MinOriginal;

//         for (int a=0; a<numberOfModels; ++a)
//             if (splinesMod[a].possibleToCalculateSpline == true)
//                 if (splinesMod[a].yD0MaxOriginal > yMax)
//                     yMax = splinesMod[a].yD0MaxOriginal;

//     }

//     if (derivativeOrder == 1) {

//         yMin = splinesExp[j].yD1Min;
//         yMax = splinesExp[j].yD1Max;

//         for (int a=0; a<numberOfModels; ++a)
//             if (splinesMod[a].possibleToCalculateSpline == true)
//                 if (splinesMod[a].yD1MinOriginal < yMin)
//                     yMin = splinesMod[a].yD1MinOriginal;

//         for (int a=0; a<numberOfModels; ++a)
//             if (splinesMod[a].possibleToCalculateSpline == true)
//                 if (splinesMod[a].yD1MaxOriginal > yMax)
//                     yMax = splinesMod[a].yD1MaxOriginal;

//     }

//     double xFirst = splinesExp[j].knots[0];
//     double xLast = splinesExp[j].knots.back();

//     double distance = (xLast-xFirst) / (double)(graphPoints-1);

//     auto xGraph = std::vector<double>(graphPoints);

//     script << "xExp = c(";
//     for (int a=0; a<splinesExp[j].abscissae.size()-1; ++a)
//         script << splinesExp[j].abscissae[a] << ", ";
//     script << splinesExp[j].abscissae.back() << ")\n";

//     script << "yExp = c(";
//     for (int a=0; a<splinesExp[j].ordinates.size()-1; ++a)
//         script << splinesExp[j].ordinates[a] << ", ";
//     script << splinesExp[j].ordinates.back() << ")\n";

//     script << "x0 = c(";
//     for (int b=0; b<graphPoints-1; ++b) {
//         xGraph[b] = xFirst+(double)b*distance;
//         script << xGraph[b] << ", ";
//     }
//     xGraph.back() = xLast;
//     script << xGraph.back() << ")\n";

//     double y;
//     script << "y0 = c(";
//     for (int b=0; b<graphPoints-1; ++b) {
//         if (derivativeOrder == 0)
//             y = splinesExp[j].D0(xGraph[b]);
//         else
//             y = splinesExp[j].D1(xGraph[b]);
//         script << y << ", ";
//     }
//     if (derivativeOrder == 0)
//         y = splinesExp[j].D0(xLast);
//     else
//         y = splinesExp[j].D1(xLast);
//     script << y << ")\n";

//     for (int a=0; a<numberOfModels; ++a)
//         if (splinesMod[a].possibleToCalculateSpline == true) {

//             if (shift == false) {

//                 double xFirst = splinesMod[a].knots[0];
//                 double xLast = splinesMod[a].knots.back();

//                 double distance = (xLast - xFirst) / (double)(graphPoints-1);

//                 auto xGraph = std::vector<double>(graphPoints);

//                 script << "x" << a+1 << " = c(";
//                 for (int b=0; b<graphPoints-1; ++b) {
//                     xGraph[b] = xFirst+(double)b*distance;
//                     script << xGraph[b] << ", ";
//                 }
//                 xGraph.back() = xLast;
//                 script << xGraph.back() << ")\n";

//                 double y;
//                 script << "y" << a+1 << " = c(";
//                 for (int b=0; b<graphPoints-1; ++b) {
//                     if (derivativeOrder == 0)
//                         y = splinesMod[a].D0(xGraph[b]);
//                     else
//                         y = splinesMod[a].D1(xGraph[b]);
//                     script << y << ", ";
//                 }
//                 if (derivativeOrder == 0)
//                     y = splinesMod[a].D0(xLast);
//                 else
//                     y = splinesMod[a].D1(xLast);
//                 script << y << ")\n";

//             }

//             if (shift == true && possibleToCalculateShifts[j][a] == true) {

//                 double shiftAmount = shiftAmounts[j][a][maxIndexes[j][a]];

//                 double xFirst = splinesMod[a].knots[0] + shiftAmount;
//                 double xLast = splinesMod[a].knots.back() + shiftAmount;

//                 double distance = (xLast - xFirst) / (double)(graphPoints-1);

//                 auto xGraph = std::vector<double>(graphPoints);

//                 script << "x" << a+1 << " = c(";
//                 for (int b=0; b<graphPoints-1; ++b) {
//                     xGraph[b] = xFirst+(double)b*distance;
//                     script << xGraph[b] << ", ";
//                 }
//                 xGraph.back() = xLast;
//                 script << xGraph.back() << ")\n";

//                 auto xGraphPowers =
//                     std::vector<std::vector<double>>(graphPoints,std::vector<double>(m,1));
//                 for (int a=0; a<graphPoints; ++a)
//                     for (int b=1; b<m; ++b)
//                         xGraphPowers[a][b] = xGraphPowers[a][b-1] * xGraph[a];

//                 splinesMod[a].normalizeCoefficients(0,1,1);
//                 splinesMod[a].calculateShift(shiftAmount);

//                 double y;
//                 script << "y" << a+1 << " = c(";
//                 for (int b=0; b<graphPoints-1; ++b) {
//                     if (derivativeOrder == 0)
//                         y = splinesMod[a].D0Shift(xGraphPowers[b]);
//                     else
//                         y = splinesMod[a].D1Shift(xGraphPowers[b]);
//                     script << y << ", ";
//                 }
//                 if (derivativeOrder == 0)
//                     y = splinesMod[a].D0Shift(xGraphPowers.back());
//                 else
//                     y = splinesMod[a].D1Shift(xGraphPowers.back());
//                 script << y << ")\n";

//             }

//         }

//     // Adds the instructions to create the graph to the .R file

//     double quality = (double)graphQuality / 1000.;

//     std::string bootstrapName = "";
//     std::string bootstrapTitle = "";
//     if (bootstrap == true) {
//         bootstrapName = "_Bootstrap-"+std::to_string(bootstrapIndex);
//         bootstrapTitle = "  Bootstrap: "+std::to_string(bootstrapIndex);
//     }

//     std::string titleEndUnderscore = "";
//     std::string titleEndSpace = "";
//     if (shift == true) {
//         titleEndUnderscore = "_shift";
//         titleEndSpace = "  shift";
//     }

//     double graphWidth;
//     double rightMargin;
//     if (bootstrap == false) {
//         graphWidth = graphQuality*1.2;
//         rightMargin = 18.*quality;
//     }
//     else {
//         if (useSumOfIndexesForAlignment == true) {
//             graphWidth = graphQuality*1.4;
//             rightMargin = 30.*quality;
//         }
//         else {
//             graphWidth = graphQuality*1.6;
//             rightMargin = 44.*quality;
//         }
//     }

//     script << "colors = rainbow(" << numberOfModels << ")\n";

//     // File name and size in pixels
//     script << "png('./" << folder << "/" << folderName << "_" << fileName
//            << bootstrapName << "_D" << derivativeOrder << titleEndUnderscore
//            << ".png', width=" << graphWidth << ", height=" << graphQuality
//            << ", type='cairo')\n";
//     // Margins
//     script << "par(mar=c(" << 4.*quality << "," << 4.*quality << ","
//            << 4.*quality << "," << rightMargin << "))\n";
//     // Title, plot details and experimental data spline
//     script << "plot(x0, y0, xlim=c(" << xMin << "," << xMax << "), ylim=c("
//            << yMin << "," << yMax << "), main='" << folderName << "  "
//            << fileName << "  D" << derivativeOrder << titleEndSpace
//            << bootstrapTitle << "', "
//            << "col='Black', type='l', xlab='', ylab='', lwd=" << 3.*quality
//            << ", " << "cex=" << quality << ", cex.axis=" << 1.5*quality
//            << ", cex.lab=" << 1.5*quality << ", cex.main=" << 1.6*quality
//            << ", mgp=c(" << 0 << "," << 2.*quality << "," << 0 << "))\n"
//            << "box(lwd=" << quality << ")\n";
//     // Axes
//     script << "axis(side=1, labels=FALSE, lwd.ticks=" << quality << ", tcl="
//            << -0.75*quality << ")\n"
//            << "axis(side=2, labels=FALSE, lwd.ticks=" << quality << ", tcl="
//            << -0.75*quality << ")\n";
//     // Experimental data points
//     if (bootstrap == true) script << "points(xExp,yExp)\n";
//     // Splines of the models
//     for (int a=0; a<numberOfModels; ++a)
//         if (splinesMod[a].possibleToCalculateSpline == true)
//             script << "lines(x" << a+1 << ", y" << a+1 << ", col=colors[" << a+1
//                    << "], lty=" << a%3+1 << ", type='l', lwd=" << 3.*quality
//                    << ")\n";
//         else
//             script << "lines(x0, y0, col='Black', type='l', lwd=" << 3.*quality
//                    << ")\n";
//     // Summary of main results
//     script << "par(xpd=TRUE)\n";
//     if (bootstrap == false) {
//         script << "legend(" << xMax+(xMax-xMin)*0.05 << ", " << (yMax+yMin)/2.
//                << ", legend=c('";
//         for (int a=0; a<numberOfModels; ++a) {
//             if (a > 0) script << ",'";
//             if (calculateShift == false) {
//                 if (scoresSorted[a] == -1)
//                     script << "NA                       ";
//                 else
//                     script << roundToString(scoresSorted[a],4) << " ± "
//                            << roundToString(stdDevSorted[a],4) << "  ";
//                 script << modelNamesSorted[a] << "'";
//             }
//             else {
//                 if (scoresSorted_shift[a] == -1)
//                     script << "NA                       ";
//                 else
//                     script << roundToString(scoresSorted_shift[a],4) << " ± "
//                            << roundToString(stdDevSorted_shift[a],4) << "  ";
//                 script << modelNamesSorted_shift[a] << "'";
//             }
//         }
//         script << "), col=c(";
//         for (int a=0; a<numberOfModels; ++a) {
//             if (a > 0) script << ",";
//             if (calculateShift == false)
//                 script << "colors[" << scoresSortedColorIndexes[a] << "]";
//             else
//                 script << "colors[" << scoresSortedColorIndexes_shift[a] << "]";
//         }
//         script << "), lty=c(";
//         for (int a=0; a<numberOfModels; ++a) {
//             if (a > 0) script << ",";
//             if (calculateShift == false)
//                 script << scoresSortedLineTypes[a];
//             else
//                 script << scoresSortedLineType_shift[a];
//         }
//         script << "), box.lwd=" << quality << ", cex=" << 1.5*quality
//                << ", lwd=" << 3.*quality << ")\n";
//     }
//     if (bootstrap == true) {
//         script << "legend(" << xMax+(xMax-xMin)*0.05 << ", " << (yMax+yMin)/2.
//                << ", legend=c('d0L2     d1L2     d0Pe     d1Pe    ";
//         useSumOfIndexesForAlignment == true ?
//             script << "shift'" :
//             script << "Sd0L2   Sd1L2  Sd0Pe   Sd1Pe'";
//         for (int a=0; a<numberOfModels; ++a) {
//             script << ",'";
//             for (int b=0; b<allIndexes[0][0].size(); ++b)
//                 script << roundToString(allIndexes.back()[a][b],4) << "  ";
//             script << modelNames[a] << "'";
//         }
//         script << "), col=c('White'";
//         for (int a=0; a<numberOfModels; ++a)
//             script << ",colors[" << a+1 << "]";
//         script << "), lty=c(1";
//         for (int a=0; a<numberOfModels; ++a)
//             script << "," << a%3+1;
//         script << "), box.lwd=" << quality << ", cex=" << 1.5*quality
//                << ", lwd=" << 3.*quality << ")\n";
//     }
//     // Legend: line types and corresponding splines
//     script << "legend(" << xMax+(xMax-xMin)*0.05 << ", "
//            << yMax+(yMax-yMin)*0.04 << ", legend=c('Exp'";
//     for (int a=0; a<numberOfModels; ++a)
//         script << ",'" << modelNames[a] << "'";
//     script << "), col=c('Black'";
//     for (int a=0; a<numberOfModels; ++a)
//         script << ",colors[" << a+1 << "]";
//     script << "), lty=c(1";
//     for (int a=0; a<numberOfModels; ++a)
//         script << "," << a%3+1;
//     script << "), box.lwd=" << quality << ", cex=" << quality << ", lwd="
//            << 4*quality << ")\ndev.off()\n\n";

//     script.close();

// }



// void Indexes::saveGraphDataToFiles() {

//     if (graphs == false) return;

//     int j = indexBestSplineExp;

//     auto x = std::vector<double>(graphPoints);
//     std::vector<double> xShift;
//     std::vector<std::vector<double>> xShiftPowers;
//     if (graphsD0Shift == true || graphsD1Shift == true) {
//         xShift = std::vector<double>(graphPoints);
//         xShiftPowers = std::vector<std::vector<double>>(graphPoints,std::vector<double>(m,1));
//     }

//     std::vector<double> y_D0, y_D1, y_D0_Shift, y_D1_Shift;
//     if (graphsD0 == true) y_D0 = std::vector<double>(graphPoints);
//     if (graphsD1 == true) y_D1 = std::vector<double>(graphPoints);
//     if (graphsD0Shift == true) y_D0_Shift = std::vector<double>(graphPoints);
//     if (graphsD1Shift == true) y_D1_Shift = std::vector<double>(graphPoints);

//     if (splinesExp[j].possibleToCalculateSpline == true) {

//         // Selects the points on the x-axis

//         double distance =
//             (splinesExp[j].knots.back()-splinesExp[j].knots[0]) /
//             (double)(graphPoints-1);

//         for (int b=0; b<graphPoints-1; ++b)
//             x[b] = splinesExp[j].knots[0]+(double)b*distance;
//         x.back() = splinesExp[j].knots.back();

//         // Calculates the ordinates

//         if (graphsD0 == true)
//             for (int b=0; b<graphPoints; ++b)
//                 y_D0[b] = splinesExp[j].D0(x[b]);

//         if (graphsD1 == true)
//             for (int b=0; b<graphPoints; ++b)
//                 y_D1[b] = splinesExp[j].D1(x[b]);

//         if (graphsD0Shift == true) {
//             if (graphsD0 == true)
//                 y_D0_Shift = y_D0;
//             else
//                 for (int b=0; b<graphPoints; ++b)
//                     y_D0_Shift[b] = splinesExp[j].D0(x[b]);
//         }

//         if (graphsD1Shift == true) {
//             if (graphsD1 == true)
//                 y_D1_Shift = y_D1;
//             else
//                 for (int b=0; b<graphPoints; ++b)
//                     y_D1_Shift[b] = splinesExp[j].D1(x[b]);
//         }

//         // Saves the results to .txt files

//         std::string pathAndPartialName =
//             "./Graph data/" + folderName + "_" + fileName + "_";

//         if (graphsD0 == true) {

//             std::string pathAndNamestring = pathAndPartialName + "Exp_D0.txt";

//             const char* pathAndName = pathAndNamestring.c_str();
//             std::ofstream script;
//             script.open(pathAndName,std::ios::app);

//             script << "x_Exp\t" << fileName << "_Exp\n";
//             for (int b=0; b<graphPoints; ++b)
//                 script << x[b] << "\t" << y_D0[b] << "\n";

//             script.close();

//         }

//         if (graphsD1 == true) {

//             std::string pathAndNamestring = pathAndPartialName + "Exp_D1.txt";

//             const char* pathAndName = pathAndNamestring.c_str();
//             std::ofstream script;
//             script.open(pathAndName,std::ios::app);

//             script << "x_Exp\t" << fileName << "_Exp\n";
//             for (int b=0; b<graphPoints; ++b)
//                 script << x[b] << "\t" << y_D1[b] << "\n";

//             script.close();

//         }

//         if (graphsD0Shift == true) {

//             std::string pathAndNamestring = pathAndPartialName + "Exp_D0_shift.txt";

//             const char* pathAndName = pathAndNamestring.c_str();
//             std::ofstream script;
//             script.open(pathAndName,std::ios::app);

//             script << "x_Exp\t" << fileName << "_Exp\n";
//             for (int b=0; b<graphPoints; ++b)
//                 script << x[b] << "\t" << y_D0_Shift[b] << "\n";

//             script.close();

//         }

//         if (graphsD1Shift == true) {

//             std::string pathAndNamestring = pathAndPartialName + "Exp_D1_shift.txt";

//             const char* pathAndName = pathAndNamestring.c_str();
//             std::ofstream script;
//             script.open(pathAndName,std::ios::app);

//             script << "x_Exp\t" << fileName << "_Exp\n";
//             for (int b=0; b<graphPoints; ++b)
//                 script << x[b] << "\t" << y_D1_Shift[b] << "\n";

//             script.close();

//         }

//     }

//     for (int a=0; a<numberOfModels; ++a)
//         if (splinesMod[a].possibleToCalculateSpline == true) {

//             // Selects the points on the x-axis

//             double distance =
//                 (splinesMod[a].knots.back()-splinesMod[a].knots[0]) /
//                 (double)(graphPoints-1);

//             for (int b=0; b<graphPoints-1; ++b)
//                 x[b] = splinesMod[a].knots[0]+(double)b*distance;
//             x.back() = splinesMod[a].knots.back();

//             if (possibleToCalculateShifts[j][a] == true)
//                 if (graphsD0Shift == true || graphsD1Shift == true) {

//                     double shift = shiftAmounts[j][a][maxIndexes[j][a]];

//                     for (int b=0; b<graphPoints; ++b)
//                         xShift[b] = x[b] + shift;

//                     for (int b=0; b<graphPoints; ++b)
//                         for (int c=1; c<m; ++c)
//                             xShiftPowers[b][c] = xShiftPowers[b][c-1]*xShift[b];

//                 }

//             // Calculates the ordinates

//             if (graphsD0 == true)
//                 for (int b=0; b<graphPoints; ++b)
//                     y_D0[b] = splinesMod[a].D0(x[b]);

//             if (graphsD1 == true)
//                 for (int b=0; b<graphPoints; ++b)
//                     y_D1[b] = splinesMod[a].D1(x[b]);

//             if (possibleToCalculateShifts[j][a] == true) {

//                 if (graphsD0Shift == true || graphsD1Shift == true) {
//                     splinesMod[a].normalizeCoefficients(0,1,1);
//                     splinesMod[a].calculateShift(
//                         shiftAmounts[j][a][maxIndexes[j][a]]);
//                 }

//                 if (graphsD0Shift == true)
//                     for (int b=0; b<graphPoints; ++b)
//                         y_D0_Shift[b] = splinesMod[a].D0Shift(xShiftPowers[b]);

//                 if (graphsD1Shift == true)
//                     for (int b=0; b<graphPoints; ++b)
//                         y_D1_Shift[b] = splinesMod[a].D1Shift(xShiftPowers[b]);

//             }

//             // Saves the results to .txt files

//             std::string pathAndPartialName =
//                 "./Graph data/" + folderName + "_" + fileName + "_";

//             if (graphsD0 == true) {

//                 std::string pathAndNamestring =
//                     pathAndPartialName + modelNames[a] + "_D0.txt";

//                 const char* pathAndName = pathAndNamestring.c_str();
//                 std::ofstream script;
//                 script.open(pathAndName,std::ios::app);

//                 script << "x_" << modelNames[a] << "\t" << fileName << "_"
//                     << modelNames[a] << "\n";
//                 for (int b=0; b<graphPoints; ++b)
//                     script << x[b] << "\t" << y_D0[b] << "\n";

//                 script.close();

//             }

//             if (graphsD1 == true) {

//                 std::string pathAndNamestring =
//                     pathAndPartialName + modelNames[a] + "_D1.txt";

//                 const char* pathAndName = pathAndNamestring.c_str();
//                 std::ofstream script;
//                 script.open(pathAndName,std::ios::app);

//                 script << "x_" << modelNames[a] << "\t" << fileName << "_"
//                     << modelNames[a] << "\n";
//                 for (int b=0; b<graphPoints; ++b)
//                     script << x[b] << "\t" << y_D1[b] << "\n";

//                 script.close();

//             }

//             if (graphsD0Shift == true &&
//                 possibleToCalculateShifts[j][a] == true) {

//                 std::string pathAndNamestring =
//                     pathAndPartialName +  modelNames[a] + "_D0_shift.txt";

//                 const char* pathAndName = pathAndNamestring.c_str();
//                 std::ofstream script;
//                 script.open(pathAndName,std::ios::app);

//                 script << "x_" << modelNames[a] << "\t" << fileName << "_"
//                     << modelNames[a] << "\n";
//                 for (int b=0; b<graphPoints; ++b)
//                     script << xShift[b] << "\t" << y_D0_Shift[b] << "\n";

//                 script.close();

//             }

//             if (graphsD1Shift == true &&
//                 possibleToCalculateShifts[j][a] == true) {

//                 std::string pathAndNamestring =
//                     pathAndPartialName + modelNames[a] + "_D1_shift.txt";

//                 const char* pathAndName = pathAndNamestring.c_str();
//                 std::ofstream script;
//                 script.open(pathAndName,std::ios::app);

//                 script << "x_" << modelNames[a] << "\t" << fileName << "_"
//                     << modelNames[a] << "\n";
//                 for (int b=0; b<graphPoints; ++b)
//                     script << xShift[b] << "\t" << y_D1_Shift[b] << "\n";

//                 script.close();

//             }

//         }

// }



int Indexes::findIndexOfMax(const std::vector<double>& sums) {

    int positionMax = -1;
    double max = sums[0];
    int k = 0;
    if (sums.size()%2 != 0) ++k;

    for (int a=0; a<(sums.size()+k)/2; ++a)
        if (sums[a] != -1)
            if (sums[a] >= max) {
                positionMax = a;
                max = sums[a];
            }
    for (int a=(sums.size()+k)/2; a<sums.size(); ++a)
        if (sums[a] != -1)
            if (sums[a] > max) {
                positionMax = a;
                max = sums[a];
            }

    // Returns -1 if 'sums' contains only -1
    if (positionMax == -1) return -1;

    return positionMax;

}



double Indexes::calculateL2(bool shift, int d0d1, Spline& splinesExp, Spline& splinesMod) {

    

    std::vector<double> y;

    if (shift == false) {

        if (d0d1 == 0)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y.push_back(yD0[a]-splinesMod.D0(xPowers[a]));
        if (d0d1 == 1)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y.push_back(yD1[a]-splinesMod.D1(xPowers[a]));

    }

    if (shift == true) {

        if (d0d1 == 0)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y.push_back(yD0[a]-splinesMod.D0Shift(xPowers[a]));
        if (d0d1 == 1)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y.push_back(yD1[a]-splinesMod.D1Shift(xPowers[a]));

    }

    double baseAlpha;
    double baseBeta;
    double integral = 0;

    baseAlpha = y[0]*y[0];
    for (int c=1; c<numberOfTrapezoids+1; ++c) {
        baseBeta = y[c]*y[c];
        integral += (baseAlpha+baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return 1./(1.+sqrt(integral)/splinesExp.xRange);

}



double Indexes::calculatePearson(bool shift, int d0d1, Spline& splinesExp, Spline& splinesMod) {

    // If the norm of the experimental data is 0, the Pearson index can't be
    // calculated
    if (normD0Exp == 0 && d0d1 == 0) return -1;
    if (normD1Exp == 0 && d0d1 == 1) return -1;

    // If the norm of the model is 0, the Pearson index can't be calculated
    if (normD0Mod == 0 && d0d1 == 0) return -1;
    if (normD1Mod == 0 && d0d1 == 1) return -1;

    //int k = indexOfModelBeingCalculated;

    //std::vector<double> y;

    //if (shift == false) {

    //    if (d0d1 == 0)
    //        for (int a=0; a<numberOfTrapezoids+1; ++a)
    //            y.push_back(yD0[a]*splinesMod[k].D0(xPowers[a]));
    //    if (d0d1 == 1)
    //        for (int a=0; a<numberOfTrapezoids+1; ++a)
    //            y.push_back(yD1[a]*splinesMod[k].D1(xPowers[a]));

    //}

    //if (shift == true) {

    //    if (d0d1 == 0)
    //        for (int a=0; a<numberOfTrapezoids+1; ++a)
    //            y.push_back(yD0[a]*splinesMod[k].D0Shift(xPowers[a]));
    //    if (d0d1 == 1)
    //        for (int a=0; a<numberOfTrapezoids+1; ++a)
    //            y.push_back(yD1[a]*splinesMod[k].D1Shift(xPowers[a]));

    //}

    //double baseAlpha;
    //double baseBeta;
    //double result = 0;

    //baseAlpha = y[0];
    //for (int c=1; c<numberOfTrapezoids+1; ++c) {
    //    baseBeta = y[c];
    //    result += (baseAlpha+baseBeta) * halfHeight;
    //    baseAlpha = baseBeta;
    //}

    //double normExp = d0d1 == 0 ? normD0Exp : normD1Exp;
    //double normMod = d0d1 == 0 ? normD0Mod : normD1Mod;

    //return  result/normExp/normMod;

    std::vector<double> y2;

    if (shift == false) {

        if (d0d1 == 0)
            for (int a=0; a<numberOfTrapezoids+1; ++a) {
                double g = splinesMod.D0(xPowers[a]);
                double y = yD0[a]/normD0Exp - g/normD0Mod;
                y2.push_back(y*y);
            }
        if (d0d1 == 1)
            for (int a=0; a<numberOfTrapezoids+1; ++a) {
                double g = splinesMod.D1(xPowers[a]);
                double y = yD1[a]/normD1Exp - g/normD1Mod;
                y2.push_back(y*y);
            }

    }

    if (shift == true) {

        if (d0d1 == 0)
            for (int a=0; a<numberOfTrapezoids+1; ++a) {
                double g = splinesMod.D0Shift(xPowers[a]);
                double y = yD0[a]/normD0Exp - g/normD0Mod;
                y2.push_back(y*y);
            }
        if (d0d1 == 1)
            for (int a=0; a<numberOfTrapezoids+1; ++a) {
                double g = splinesMod.D1Shift(xPowers[a]);
                double y = yD1[a]/normD1Exp - g/normD1Mod;
                y2.push_back(y*y);
            }

    }

    double baseAlpha;
    double baseBeta;
    double result = 0;

    baseAlpha = y2[0];
    for (int c=1; c<numberOfTrapezoids+1; ++c) {
        baseBeta = y2[c];
        result += (baseAlpha+baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return  1.-sqrt(result)/2.;

}



double Indexes::calculateNorm(Spline& spline, bool shift, int d0d1) {

    std::vector<double> y2;
    // only if splineType is not model, so only for experiments!
    // since it did the yD0 shit already it does not repeat it for the experiments!

    if (spline.splineType != 1 /*Model*/) {

        if (d0d1 == 0)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y2.push_back(yD0[a]*yD0[a]);
        if (d0d1 == 1)
            for (int a=0; a<numberOfTrapezoids+1; ++a)
                y2.push_back(yD1[a]*yD1[a]);

    }

    else {

        if (shift == false) {

            if (d0d1 == 0)
                for (int a=0; a<numberOfTrapezoids+1; ++a) {
                    double y = spline.D0(xPowers[a]);
                    y2.push_back(y*y);
                }
            if (d0d1 == 1)
                for (int a=0; a<numberOfTrapezoids+1; ++a) {
                    double y = spline.D1(xPowers[a]);
                    y2.push_back(y*y);
                }

        }

        if (shift == true) {

            if (d0d1 == 0)
                for (int a=0; a<numberOfTrapezoids+1; ++a) {
                    double y = spline.D0Shift(xPowers[a]);
                    y2.push_back(y*y);
                }
            if (d0d1 == 1)
                for (int a=0; a<numberOfTrapezoids+1; ++a) {
                    double y = spline.D1Shift(xPowers[a]);
                    y2.push_back(y*y);
                }

        }

    }

    double baseAlpha;
    double baseBeta;
    double result = 0;

    baseAlpha = y2[0];
    for (int c=1; c<numberOfTrapezoids+1; ++c) {
        baseBeta = y2[c];
        if (d0d1 == 0)
            result += (baseAlpha+baseBeta) * halfHeight;
        if (d0d1 == 1)
            result += (baseAlpha+baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return sqrt(result);

}



void Indexes::sortScores(std::vector<double> scoresIn,
                         std::vector<double>& scoresOut,
                         std::vector<double> stdDevIn,
                         std::vector<double>& stdDevOut,
                         std::vector<std::string>& scoresOutNames,
                         std::vector<int>& scoresOutLineTypes,
                         std::vector<int>& scoresOutColorIndexes) {

    std::vector<std::string> names = modelNames;

    std::vector<int> lineTypes;
    for (int a=0; a<numberOfModels; ++a)
        lineTypes.push_back(a%3+1);

    std::vector<int> colorIndexes;
    for (int a=1; a<numberOfModels+1; ++a)
        colorIndexes.push_back(a);

    while (scoresIn.size() > 0) {

        double max = -1;
        int indexMax;

        for (int a=0; a<scoresIn.size(); ++a)
            if (scoresIn[a] > max) {
                max = scoresIn[a];
                indexMax = a;
            }

        if (max != -1) { // If at least one score is not NA
            scoresOut.push_back(max);
            stdDevOut.push_back(stdDevIn[indexMax]);
            scoresOutNames.push_back(names[indexMax]);
            scoresOutLineTypes.push_back(lineTypes[indexMax]);
            scoresOutColorIndexes.push_back(colorIndexes[indexMax]);
            scoresIn.erase(scoresIn.begin()+indexMax);
            stdDevIn.erase(stdDevIn.begin()+indexMax);
            names.erase(names.begin()+indexMax);
            lineTypes.erase(lineTypes.begin()+indexMax);
            colorIndexes.erase(colorIndexes.begin()+indexMax);
        }
        else {
            scoresOut.push_back(scoresIn[0]);
            stdDevOut.push_back(stdDevIn[0]);
            scoresOutNames.push_back(names[0]);
            scoresOutLineTypes.push_back(lineTypes[0]);
            scoresOutColorIndexes.push_back(colorIndexes[0]);
            scoresIn.erase(scoresIn.begin());
            stdDevIn.erase(stdDevIn.begin());
            names.erase(names.begin());
            lineTypes.erase(lineTypes.begin());
            colorIndexes.erase(colorIndexes.begin());
        }

    }

}



void Indexes::calculateShiftAmounts(Spline& splinesExp, Spline& splinesMod) {

    //int i = splineExpIndex;
    //int a = modelNumber;

    double expLeft = splinesExp.xMinNoAsymptotes;
    double expRight = splinesExp.xMaxNoAsymptotes;

    std::vector<double> acceptedShifts;

    bool possibleToLineUpMaxima = false;

    if (lineUpMaxima == true) {

        double shift;

        if (splinesExp.hasOneMaximum == true)
            if (splinesMod.hasOneMaximum == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp.xOfMax - splinesMod.xOfMax;
            }

        if (splinesExp.hasOneMaximum == true)
            if (splinesMod.hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                double xOfClosestPoint = splinesMod.xOfFirstMax;
                if (fabs(splinesExp.xOfMax-splinesMod.xOfSecondMax) <
                    fabs(splinesExp.xOfMax-splinesMod.xOfFirstMax))
                    xOfClosestPoint = splinesMod.xOfSecondMax;
                shift = splinesExp.xOfMax - xOfClosestPoint;
            }

        if (splinesExp.hasTwoMaxima == true)
            if (splinesMod.hasOneMaximum == true) {
                possibleToLineUpMaxima = true;
                double xOfClosestPoint = splinesExp.xOfFirstMax;
                if (fabs(splinesMod.xOfMax-splinesExp.xOfSecondMax) <
                    fabs(splinesMod.xOfMax-splinesExp.xOfFirstMax))
                    xOfClosestPoint = splinesExp.xOfSecondMax;
                shift = xOfClosestPoint - splinesMod.xOfMax;
            }

        if (splinesExp.hasTwoMaxima == true)
            if (splinesMod.hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp.xOfFirstMax-splinesMod.xOfFirstMax;
            }

        if (splinesExp.hasTwoMaxima == true)
            if (splinesMod.hasFirstOfTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp.xOfFirstMax-splinesMod.xOfMax;
            }

        if (splinesExp.hasFirstOfTwoMaxima == true)
            if (splinesMod.hasFirstOfTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp.xOfMax-splinesMod.xOfMax;
            }

        if (splinesExp.hasFirstOfTwoMaxima == true)
            if (splinesMod.hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp.xOfMax-splinesMod.xOfFirstMax;
            }

        if (possibleToLineUpMaxima == true) {

            acceptedShifts.push_back(shift-shiftAroundMaximum);
            double limit = shift + shiftAroundMaximum;
            while (acceptedShifts.back() < limit)
                acceptedShifts.push_back(
                    acceptedShifts.back()+distanceBetweenShiftedPoints);

        }

        if (acceptedShifts.size() == 0)
            possibleToLineUpMaxima = false;

    }

    if (possibleToLineUpMaxima == false) {

        std::vector<double> shiftsvector;

        double modLeft = splinesMod.xMinNoAsymptotes;
        double modRight = splinesMod.xMaxNoAsymptotes;

        shiftsvector.push_back(expLeft-modRight);
        while (modLeft+shiftsvector.back() < expRight)
            shiftsvector.push_back(
                shiftsvector.back()+distanceBetweenShiftedPoints);

        // If the model has either one or two maxima, the shifts which would
        // make its leftmost original abscissa < 0 are accepted only if at least
        // one maximum remains in the experimental data range after the shift

        std::vector<double> relevantShifts;

        if (shiftsvector.size() > 0)
            for (int b=0; b<shiftsvector.size(); ++b) {

                if (splinesMod.originalAbscissae[0]+shiftsvector[b] >= 0)
                    relevantShifts.push_back(shiftsvector[b]);

                if (splinesMod.originalAbscissae[0]+shiftsvector[b] < 0) {

                    if (splinesMod.hasOneMaximum == true) {
                        double max = splinesMod.xOfMax + shiftsvector[b];
                        if (max >= expLeft && max <= expRight)
                            relevantShifts.push_back(shiftsvector[b]);
                    }

                    if (splinesMod.hasTwoMaxima == true) {
                        bool atLeastOneMaximumWithinExpRange = false;
                        double maxLeft =
                            splinesMod.xOfFirstMax + shiftsvector[b];
                        double maxRight =
                            splinesMod.xOfSecondMax + shiftsvector[b];
                        if (maxLeft >= expLeft && maxLeft <= expRight)
                            atLeastOneMaximumWithinExpRange = true;
                        if (maxRight >= expLeft && maxRight <= expRight)
                            atLeastOneMaximumWithinExpRange = true;
                        if (atLeastOneMaximumWithinExpRange == true)
                            relevantShifts.push_back(shiftsvector[b]);
                    }

                    if (splinesMod.hasOneMaximum == false &&
                        splinesMod.hasTwoMaxima == false)
                        relevantShifts.push_back(shiftsvector[b]);

                }

            }

        // Checks whether the calculated shift amounts are compatible with the
        // limit on model extrapolation
        if (relevantShifts.size() > 0)
            for (int b=0; b<relevantShifts.size(); ++b) {
                double extrapolatedLeft =
                    splinesMod.originalAbscissae[0] +
                    relevantShifts[b] - expLeft;
                if (extrapolatedLeft < 0)
                    extrapolatedLeft = 0;
                double extrapolatedRight =
                    expRight - relevantShifts[b] -
                    splinesMod.originalAbscissae.back();
                if (extrapolatedRight < 0)
                    extrapolatedRight = 0;
                if (extrapolatedLeft+extrapolatedRight <= maxExtrapolatedLength)
                    acceptedShifts.push_back(relevantShifts[b]);
            }

    } // End of if (possibleToLineUpMaxima == false)

    if (acceptedShifts.size() > 0) {
        possibleToCalculateShifts[0][0] = true;
        shiftAmounts[0][0] = acceptedShifts;
    }

}



void Indexes::statisticalAnalysis(const std::vector<std::vector<double>>& bootstrap,
                                  int i, /*index of the model in modelNames*/
                                  double& mean,
                                  double& stdDeviation) {

    double n = (double)numberOfBootstrapVariations;

    // Calculates the mean
    double sum = 0;
    for (int a=0; a<numberOfBootstrapVariations; ++a)
        sum += bootstrap[a][i];
    mean = sum/n;

    // Calculates the standard deviation of the mean
    double sumOfSquares = 0;
    for (int a=0; a<numberOfBootstrapVariations; ++a) {
        double difference = bootstrap[a][i]-mean;
        sumOfSquares += difference*difference;
    }
    stdDeviation = sqrt(sumOfSquares/(n-1.));

}



void Indexes::saveIndexesToCsvFile() {

    std::string fileNameWithPath = "./Indexes/"+folderName+"_"+fileName+".csv";

    std::ofstream file(fileNameWithPath.c_str());

    for (int a=0; a<numberOfBootstrapVariations; ++a) {

        file << ",d0L2,d1L2,d0Pe,d1Pe,";
        useSumOfIndexesForAlignment == true ?
            file << "shift\n" :
            file << "shift_d0L2,shift_d1L2,shift_d0Pe,shift_d1Pe\n";

        for (int b=0; b<numberOfModels; ++b) {
            file << modelNames[b] << ",";
            for (int c=0; c<allIndexes[a][b].size()-1; ++c)
                file << allIndexes[a][b][c] << ",";
            file << allIndexes[a][b][allIndexes[a][b].size()-1] << "\n";
        }

        if (a != numberOfBootstrapVariations-1)
            file << "\n";

    }

    file.close();

}



std::string Indexes::roundToString(double number, int decimalPlaces) {

    if (number == -1) return "NA      ";

    double tens = 1.;
    for (int a=0; a<decimalPlaces; ++a)
        tens *= 10.;
    double roundedNumber = round(tens*number)/tens;

    std::string numberStr = std::to_string(roundedNumber);
    numberStr.erase(numberStr.begin()+(2+decimalPlaces),numberStr.end());

    return numberStr;

}
