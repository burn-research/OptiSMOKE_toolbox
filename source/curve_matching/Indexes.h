using namespace std;
using namespace std::chrono;

class Indexes {

public:

    /* Name of the _exp.txt and _mod.txt files with the data for the splines */
    string fileName;

    /* Number of models to be compared with the experimental data */
    int numberOfModels;

    /* Names of the models */
    vector<string> modelNames;

    /* Scores, before the shift. scores[i] contains the data relative to the
    comparison of modelNames[i] and the experimental data. scores[i] is equal to
    -1 in case it is not possible to calculate score i, or in case the
    calculated score is nan */
    vector<double> scores;

    /* Scores, after the shift. scores_shift[i] contains the data relative to
    the comparison of modelNames[i] and the experimental data. scores_shift[i]
    is equal to -1 in case it is not possible to calculate score i, or in case
    the calculated score is nan */
    vector<double> scores_shift;

    /* Standard deviations of the scores, before the shift. stdDev[i] contains
    the data relative to the comparison of modelNames[i] and the experimental
    data. stdDev[i] is equal to -1 in case it is not possible to calculate the
    corresponding score, or it is nan */
    vector<double> stdDev;

    /* Standard deviations of the scores, after the shift. stdDev_shift[i]
    contains the data relative to the comparison of modelNames[i] and the
    experimental data. stdDev_shift[i] is equal to -1 in case it is not possible
    to calculate the corresponding score, or it is nan */
    vector<double> stdDev_shift;

    /* Individual indexes for each variation of the experimental data obtained
    with the bootstrapping procedure. allIndexes[i][j][k]: i: variations of the
    experimental data; j: models; k: d0L2, d1L2, d0Pe, d1Pe, shift (or shifts if
    all different) */
    vector<vector<vector<double>>> allIndexes;

    ////////////////////////////////////////////////////////////////////////////

    /* Executes all the operations for the comparison of the experimental data
    with the models */
    void solve(bool calculateShift);

    ////////////////////////////////////////////////////////////////////////////

    /* Executes all the operations for the comparison of the experimental data
    with the models */
    void solvePython(bool CalculateShift,
                     vector<double> exp_data_x,
                     vector<double> exp_data_y,
                     vector<double> model_data_x,
                     vector<double> model_data_y,
                     vector<double> error_data);

private:

    /* Equivalent of 'scores' for every bootstrap variation */
    vector<vector<double>> bootstrapScores;

    /* Equivalent of scores_shift for every bootstrap variation */
    vector<vector<double>> bootstrapScores_shift;

    /* x-coordinates and y-coordinates of the experimental data and the models.
    inputData[0] and inputData[1] contain, respectively, the x-coordinates and
    the y-coordinates of the experimental data, inputData[i] and inputData[i+1]
    contain, respectively, the x-coordinates and the y-coordinates of
    modelNames[i/2-1] */
    vector<vector<double>> inputData;

    /* Splines for the experimental data */
    vector<Spline> splinesExp;

    /* Splines for the models. splinesMod[i] contains the spline for
    modelNames[i] */
    vector<Spline> splinesMod;

    /* Splines for the models, before being extended */
    vector<Spline> splinesMod_original;

    /* Path to the folder containing the input data */
    string folderPath;

    /* Name of the folder containing the input data */
    string folderName;

    /* Relative errors for each experimental data point */
    vector<double> relativeErrors;

    /* Specifies whether to calculate the indexes for the shifted splines */
    bool calculateShift;

    /* Amounts to be added to the abscissae of the models to find the shifted
    abscissae. shiftAmounts[i][j] refers to splinesExp[i] and modelNames[j] */
    vector<vector<vector<double>>> shiftAmounts;

    /* Positions in shiftAmounts, indicating the shifts that maximise the mean
    of d0L2, d1L2, d0Pe and d1Pe for the models. maxIndexes[i][j] refers to
    shiftAmounts[i][j] */
    vector<vector<int>> maxIndexes;

    /* Index in the splinesMod vector of the model being calculated */
    int indexOfModelBeingCalculated;

    /* Height of the trapezoids used used during integration, divided by 2 */
    double halfHeight;

    /* Powers of the abscissae of the sides of the trapezoids used for the
    numerical calculation of the indexes, distributed evenly between the
    extremes of the experimental data */
    vector<vector<double>> xPowers;

    /* Ordinates of the experimental data spline, with yD0[i] corresponding to
    xPowers[i][1] */
    vector<double> yD0;

    /* Ordinates of the first derivative of the experimental data spline, with
    yD1[i] corresponding to xPowers[i][1] */
    vector<double> yD1;

    /* Scores without shifts. Sorted from lowest to highest. Is equal to -1 in
    case it is not possible to calculate a score */
    vector<double> scoresSorted;

    /* Standard deviations without shifts. Sorted from lowest to highest. Is
    equal to -1 in case it is not possible to calculate a score */
    vector<double> stdDevSorted;

    /* Model names, sorted as the elements of scoresSorted */
    vector<string> modelNamesSorted;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted */
    vector<int> scoresSortedLineTypes;

    /* Color indexes in the .R graphs, sorted as the elements of scoresSorted */
    vector<int> scoresSortedColorIndexes;

    /* Scores with shifts. Sorted from lowest to highest. Is equal to -1 in case
    it is not possible to calculate a score */
    vector<double> scoresSorted_shift;

    /* Standard deviations with shifts. Sorted from lowest to highest. Is equal
    to -1 in case it is not possible to calculate a score */
    vector<double> stdDevSorted_shift;

    /* Model names, sorted as the elements of scoresSorted_shift */
    vector<string> modelNamesSorted_shift;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted_shift
    */
    vector<int> scoresSortedLineType_shift;

    /* Color indexes in the .R graphs, sorted as the elements of
    scoresSorted_shift */
    vector<int> scoresSortedColorIndexes_shift;

    /* Index in the splinesExp vector of the spline being considered */
    int splineExpIndex;

    /* Index in the splinesExp vector of the best fit of the experimental data
    */
    int indexBestSplineExp;

    /* Values of the d0L2, d1L2, d0Pe and d1Pe indexes, in this order, obtained
    before any shift. indexes[i][j] refers to splinesExp[i] and contains the
    data relative to the comparison of modelNames[j] and the experimental data.
    Is equal to -1 in case it is not possible to calculate an index */
    vector<vector<vector<double>>> indexes;

    /* Values of the d0L2_shift, d1L2_shift, d0Pe_shift, d1Pe_shift indexes, in
    this order, corresponding to their optimal values found while shifting.
    indexes_shift[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate an index */
    vector<vector<vector<double>>> indexes_shift;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximize the means of d0L2, d1L2, d0Pe and d1Pe for each
    model. shifts[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate a shift */
    vector<vector<double>> shifts;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximise each of d0L2, d1L2, d0Pe and d1Pe for each model.
    shifts[i][j] refers to splinesExp[i] and contains the data relative to the
    comparison of modelNames[j] and the experimental data. Is equal to -1 in
    case it is not possible to calculate a shift */
    vector<vector<vector<double>>> shifts_allDifferent;

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
    vector<vector<bool>> possibleToCalculateShifts;

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
    ////////////////////////////////////////////////////////////////////////////



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
    void calculateIndexes();

    /* Calculates the scores from the indexes */
    void calculateScores();

    /* Saves the data for creating the graph of each spline and the
    corresponding data points and knots to a .R file */
    void graphsFittingR(vector<Spline> &splines);

    /* Saves the data for creating the graph comparing the experimental data
    with the models (either shifted or non-shifted) to a .R file */
    void graphR(int derivativeOrder, bool shift);

    /* Saves the data necessary to plot the splines to .txt files */
    void saveGraphData();

    /* Given a vector containing the sums of d0L2, d1L2, d0Pe and d1Pe at
    various shift locations, returns the position of the maximum sum,
    corresponding to the shift in the same position in shiftAmounts. If there is
    more than one absolute maximum, returns the position of the maximum
    referring to the smallest shift. If the input vector contains only -1 (each
    meaning the mean could not be calculated), returns -1 */
    int findIndexOfMax(const vector<double> &sums);

    /* Calculates the L2 index */
    double calculateL2(bool shift, int d0d1);

    /* Calculates the Pearson index */
    double calculatePearson(bool shift, int d0d1);

    /* Calculates either the norm of a spline or of its first derivative, with
    or without a shift along the x-axis */
    double calculateNorm(Spline &spline, bool shift, int d0d1);

    /* Sorts the scoresIn vector from highest to lowest, and sorts stdDevIn
    accordingly. Saves the results to scoresOut and stdDevOut. scoresOutNames,
    scoresOutLineTypes and scoresOutColorIndexes are also sorted as the elements
    of scoresOut */
    void sortScores(vector<double> scoresIn,
                    vector<double> &scoresOut,
                    vector<double> stdDevIn,
                    vector<double> &stdDevOut,
                    vector<string> &scoresOutNames,
                    vector<int> &scoresOutLineTypes,
                    vector<int> &scoresOutColorIndexes);

    /* Calculates the shift amounts for each model */
    void calculateShiftAmounts(int modelNumber);

    /* Calculates the mean and the standard deviation of the mean of the
    bootstrap scores of a model, from either bootstrapScores or
    bootstrapScores_shift */
    void statisticalAnalysis(const vector<vector<double>> &bootstrap,
                             int i, /*index of the model in modelNames*/
                             double &mean,
                             double &stdDeviation);

    /* Saves the individual indexes for all the bootstrap variations to a .csv
    file */
    void saveIndexesToCsvFile();

    void readDataPython(
            vector<double> exp_data_x,
            vector<double> exp_data_y,
            vector<double> model_data_x,
            vector<double> model_data_y,
            vector<double> error_data);

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Indexes::solvePython(bool CalculateShift,
                          vector<double> exp_data_x,
                          vector<double> exp_data_y,
                          vector<double> model_data_x,
                          vector<double> model_data_y,
                          vector<double> error_data) {

    // Generates the Pascal's triangle
    pascalsTriangle = vector<vector<double>>(m, vector<double>(1, 1));

    pascalsTriangle[1].push_back(1);
    if (g > 1) {
        for (int a = 2; a < m; ++a) {
            for (int b = 0; b < a - 1; ++b)
                pascalsTriangle[a].push_back(
                        pascalsTriangle[a - 1][b] + pascalsTriangle[a - 1][b + 1]);
            pascalsTriangle[a].push_back(1);
        }
    }

    this->readDataPython(exp_data_x,
                         exp_data_y,
                         model_data_x,
                         model_data_y,
                         error_data);

    this->solve(CalculateShift);

//    for (int i = 0; i < numberOfBootstrapVariations; i++) {
//        for (int j = 0; j < numberOfModels; j++) {
//            for (int k = 0; k < 5; k++) {
//                int index = ((numberOfModels * i) + j) + (k * (numberOfBootstrapVariations * numberOfModels));
//                results_indexes[index] = allIndexes[i][j][k];
//            }
//        }
//    }

}


void Indexes::solve(bool CalculateShift) {


//    folderPath = FolderPath;
//    folderName = FolderName;
//    fileName = FileName;
    calculateShift = CalculateShift;

    // Processes the data in the input files
    // this->readData(); ADD to bring back


    // Creates a matrix with each line equal to a plausible version of the
    // ordinates of the experimental data, generated randomly by taking into
    // account the measurement error. The first line contains the ordinates
    // provided in the input files, and not randomly generated ordinates

    auto bootstrapExp = vector<vector<double>>(numberOfBootstrapVariations);
    for (int a = 0; a < numberOfBootstrapVariations; ++a)
        bootstrapExp[a] = inputData[1];

    for (int a = 0; a < inputData[1].size(); ++a) {

        // Calculates the standard deviation
        double stdDeviation = inputData[1][a] * relativeErrors[a];

        default_random_engine generator;
//        generator.seed(system_clock::now().time_since_epoch().count());
        generator.seed(0);

        normal_distribution<double> distribution(inputData[1][a], stdDeviation);

        for (int b = 1; b < numberOfBootstrapVariations; ++b) {
            double number = distribution(generator);
            if (number < 0 && possibleNegativeOrdinates == false) number = 0;
            bootstrapExp[b][a] = number;
        }

    }

    bootstrapScores = vector<vector<double>>(numberOfBootstrapVariations);
    bootstrapScores_shift = vector<vector<double>>(numberOfBootstrapVariations);

    // Calculates the splines for the models
    numberOfModels = inputData.size() / 2 - 1;
    splinesMod = vector<Spline>(numberOfModels);
    for (int a = 0; a < numberOfModels; ++a)
        splinesMod[a].solve(inputData[2 * (a + 1)], inputData[2 * (a + 1) + 1], 1, 0);

    splinesMod_original = splinesMod;


    for (int a = numberOfBootstrapVariations - 1; a > -1; --a) {

        bootstrapIndex = a;
        inputData[1] = bootstrapExp[bootstrapIndex];

        // Calculates the splines for the experimental data
        this->calculateSplines();


        // Initializes several elements used in the subsequent calculations

        possibleToCalculateShifts = vector<vector<bool>>(splinesExp.size(),
                                                         vector<bool>(numberOfModels, false));

        indexes = vector<vector<vector<double>>>(splinesExp.size(),
                                                 vector<vector<double>>(numberOfModels,
                                                                        vector<double>(4, -1)));

        indexes_shift = vector<vector<vector<double>>>(splinesExp.size(),
                                                       vector<vector<double>>(numberOfModels,
                                                                              vector<double>(4, -1)));

        shifts = vector<vector<double>>(splinesExp.size(),
                                        vector<double>(numberOfModels, -1));

        shifts_allDifferent = vector<vector<vector<double>>>(splinesExp.size(),
                                                             vector<vector<double>>(numberOfModels,
                                                                                    vector<double>(4, -1)));

        shiftAmounts = vector<vector<vector<double>>>(splinesExp.size(),
                                                      vector<vector<double>>(numberOfModels));

        maxIndexes = vector<vector<int>>(splinesExp.size(),
                                         vector<int>(numberOfModels, -1));


        // For each knot arrangement in the experimental data, calculates the
        // indexes for each model
        for (int b = 0; b < splinesExp.size(); ++b) {
            splineExpIndex = b;
            this->calculateIndexes();
        }

        // Calculates the scores from the indexes
        this->calculateScores();

    }


    scores = vector<double>(numberOfModels, -1);
    scores_shift = vector<double>(numberOfModels, -1);
    stdDev = vector<double>(numberOfModels, -1);
    stdDev_shift = vector<double>(numberOfModels, -1);

    for (int a = 0; a < numberOfModels; ++a) {

        // No shift:
        bool isNA = false;
        for (int b = 0; b < numberOfBootstrapVariations; ++b)
            if (bootstrapScores[b][a] == -1) {
                isNA = true;
                break;
            }
        if (isNA == false)
            statisticalAnalysis(bootstrapScores, a,
                                scores[a], stdDev[a]);

        // Shift:
        isNA = false;
        for (int b = 0; b < numberOfBootstrapVariations; ++b)
            if (bootstrapScores_shift[b][a] == -1) {
                isNA = true;
                break;
            }
        if (isNA == false)
            statisticalAnalysis(bootstrapScores_shift, a,
                                scores_shift[a], stdDev_shift[a]);

    }

    sortScores(scores,
               scoresSorted,
               stdDev,
               stdDevSorted,
               modelNamesSorted,
               scoresSortedLineTypes,
               scoresSortedColorIndexes);

    sortScores(scores_shift,
               scoresSorted_shift,
               stdDev_shift,
               stdDevSorted_shift,
               modelNamesSorted_shift,
               scoresSortedLineType_shift,
               scoresSortedColorIndexes_shift);

}


void Indexes::readDataPython(vector<double> exp_data_x,
                             vector<double> exp_data_y,
                             vector<double> model_data_x,
                             vector<double> model_data_y,
                             vector<double> error_data) {

    // The relative errors and the experimental and model data are given
    // Sorted and with no duplicates

    modelNames.push_back("tmp");

    // 4 Vectors: 2 from experiment, 2 from model
    vector<vector<double>> allInputData;
    vector<double> allRelativeErrors;

    vector<double> exp_data_x_vector = exp_data_x;
    vector<double> exp_data_y_vector = exp_data_y;

    vector<double> model_data_x_vector = model_data_x;
    vector<double> model_data_y_vector = model_data_y;

    allInputData.push_back(exp_data_x_vector);
    allInputData.push_back(exp_data_y_vector);
    allInputData.push_back(model_data_x_vector);
    allInputData.push_back(model_data_y_vector);

    for (int i = 0; i < exp_data_x.size(); i++) {
        allRelativeErrors.push_back(error_data[i]);
    }


    // Obtains inputData and relativeErrors by removing excess points with
    // ordinate 0 on the left and on the right of experimental data and models

    inputData = vector<vector<double>>(4);


    for (int a = 0; a < allInputData.size(); ++a)
        if (allInputData[a].size() > 1) {

            if (a == 0) {

                double maxOrdinate = allInputData[a + 1][0];
                double minOrdinate = allInputData[a + 1][0];
                for (int b = 1; b < allInputData[a].size(); ++b) {
                    if (allInputData[a + 1][b] > maxOrdinate)
                        maxOrdinate = allInputData[a + 1][b];
                    if (allInputData[a + 1][b] < minOrdinate)
                        minOrdinate = allInputData[a + 1][b];
                }
                double height = maxOrdinate - minOrdinate;
                double minHeight =
                        fractionOfOrdinateRangeForAsymptoteIdentification * height;

                int indexOfFirstPoint = 0;
                for (int b = 1; b < allInputData[a].size(); ++b)
                    if (fabs(allInputData[a + 1][b] - allInputData[a + 1][b - 1]) >
                        minHeight) {
                        indexOfFirstPoint = b - 1;
                        break;
                    }
                int indexOfLastPoint = allInputData[a].size() - 1;
                for (int b = allInputData[a].size() - 2; b > -1; --b)
                    if (fabs(allInputData[a + 1][b] - allInputData[a + 1][b + 1]) >
                        minHeight) {
                        indexOfLastPoint = b + 1;
                        break;
                    }

                inputData[a] =
                        vector<double>(indexOfLastPoint - indexOfFirstPoint + 1);
                inputData[a + 1] =
                        vector<double>(indexOfLastPoint - indexOfFirstPoint + 1);

                for (int b = 0; b < inputData[a].size(); ++b)
                    inputData[a][b] = allInputData[a][indexOfFirstPoint + b];

                for (int b = 0; b < inputData[a].size(); ++b)
                    inputData[a + 1][b] = allInputData[a + 1][indexOfFirstPoint + b];

                for (int b = indexOfFirstPoint; b <= indexOfLastPoint; ++b)
                    relativeErrors.push_back(allRelativeErrors[b]);

            }

            if (a > 0 && a % 2 == 0) {

                inputData[a] = allInputData[a];
                inputData[a + 1] = allInputData[a + 1];

            }

        }


}


void Indexes::calculateSplines() {

    if (inputData[0].size() < 3)
        splinesExp = vector<Spline>(1);
    else if (inputData[0].size() < 5)
        splinesExp = vector<Spline>(2);
    else
        splinesExp = vector<Spline>(3);

    splinesExp[0].solve(inputData[0], inputData[1], 0, 0);
    splinesExp[0].removeAsymptotes();

    if (splinesExp.size() > 1) {
        splinesExp[1].solve(inputData[0], inputData[1], 0, 2);
        splinesExp[1].removeAsymptotes();
    }

    if (splinesExp.size() > 2) {
        splinesExp[2].solve(inputData[0], inputData[1], 0, 5);
        splinesExp[2].removeAsymptotes();
    }

}


void Indexes::calculateIndexes() {

    int i = splineExpIndex;

    splinesMod = splinesMod_original;

    // Checks whether any model can be considered a flat line with ordinate = 0
    // when compared to the experimental data
    if (possibleNegativeOrdinates == false) {
        double minHeight =
                splinesExp[i].yD0Max * fractionOfExpHeightForZeroLineIdentification;
        for (int a = 0; a < numberOfModels; ++a)
            if (splinesMod[a].possibleToCalculateSpline == true)
                if (splinesMod[a].yD0Max < minHeight)
                    splinesMod[a].possibleToCalculateSpline = false;
    }

    double minLengthInCommon =
            splinesExp[i].xRange * fractionOfExpRangeForModelExtrapolation;

    maxExtrapolatedLength = splinesExp[i].xRange - minLengthInCommon;

    // Adds segments to the left and to the right of the models so that their
    // ordinates span the entirety of the experimental data range, regardless of
    // the shift value
    for (int a = 0; a < numberOfModels; ++a)
        if (splinesMod[a].possibleToCalculateSpline == true)
            splinesMod[a].extendSpline(maxExtrapolatedLength,
                                       maxExtrapolatedLength);

    // Normalizes the spline coefficients of the experimental data and of the
    // first derivative of the experimental data
    splinesExp[i].normalizeCoefficients(
            -splinesExp[i].yD0Min,
            splinesExp[i].yD0Max - splinesExp[i].yD0Min,
            splinesExp[i].yD1MaxAbs);

    // Calculates the values of shiftAmounts
    if (calculateShift == true) {

        shiftAroundMaximum =
                splinesExp[i].xRange * fractionOfExpRangeForShiftAroundMaximum;

        distanceBetweenShiftedPoints =
                splinesExp[i].xRange * fractionOfExpRangeForMinShift;

        for (int a = 0; a < numberOfModels; ++a)
            if (splinesMod[a].possibleToCalculateSpline == true)
                calculateShiftAmounts(a);

    }

    // Initializes several elements necessary for the calculation of the indexes

    double height = (splinesExp[i].knots.back() - splinesExp[i].knots[0]) /
                    (double) numberOfTrapezoids;
    halfHeight = height / 2.;

    vector<double> x;
    for (int a = 0; a < numberOfTrapezoids; ++a)
        x.push_back(splinesExp[i].knots[0] + (double) a * height);
    x.push_back(splinesExp[i].knots.back());

    xPowers = vector<vector<double>>(numberOfTrapezoids + 1, vector<double>(m, 1));
    for (int a = 0; a < numberOfTrapezoids + 1; ++a)
        for (int b = 1; b < m; ++b)
            xPowers[a][b] = xPowers[a][b - 1] * x[a];

    yD0 = vector<double>(numberOfTrapezoids + 1);
    for (int a = 0; a < numberOfTrapezoids + 1; ++a)
        yD0[a] = splinesExp[i].D0(xPowers[a]);

    yD1 = vector<double>(numberOfTrapezoids + 1);
    for (int a = 0; a < numberOfTrapezoids + 1; ++a)
        yD1[a] = splinesExp[i].D1(xPowers[a]);

    vector<double> xL2Means;
    for (int a = 0; a < numberOfTrapezoids; ++a)
        xL2Means.push_back((x[a] + x[a + 1]) / 2.);

    auto xL2MeansPowers =
            vector<vector<double>>(numberOfTrapezoids, vector<double>(m, 1));
    for (int a = 0; a < numberOfTrapezoids; ++a)
        for (int b = 1; b < m; ++b)
            xL2MeansPowers[a][b] = xL2MeansPowers[a][b - 1] * xL2Means[a];

    // Calculates the norms of the experimental data
    normD0Exp = calculateNorm(splinesExp[i], false, 0);
    normD1Exp = calculateNorm(splinesExp[i], false, 1);

    // Calculates the indexes for each model, if possible
    for (int a = 0; a < numberOfModels; ++a) {
        if (splinesMod[a].possibleToCalculateSpline == true) {

            indexOfModelBeingCalculated = a;

            // Normalizes the coefficients of the spline and of the first
            // derivative of the spline
            splinesMod[a].normalizeCoefficients(
                    -splinesExp[i].yD0Min,
                    splinesExp[i].yD0Max - splinesExp[i].yD0Min,
                    splinesExp[i].yD1MaxAbs);

            // Checks whether the experimental data and the model have enough
            // non-extrapolated lenght in common
            bool possibleToCalculateIndexes = true;
            double limitLeft = splinesMod[a].originalAbscissae[0];
            if (limitLeft < splinesExp[i].xMinNoAsymptotes)
                limitLeft = splinesExp[i].xMinNoAsymptotes;
            double limitRight = splinesMod[a].originalAbscissae.back();
            if (limitRight > splinesExp[i].xMaxNoAsymptotes)
                limitRight = splinesExp[i].xMaxNoAsymptotes;
            if (limitRight - limitLeft < minLengthInCommon)
                possibleToCalculateIndexes = false;

            if (possibleToCalculateIndexes == true) {

                // Calculates the norms for the model
                normD0Mod = calculateNorm(splinesMod[a], false, 0);
                normD1Mod = calculateNorm(splinesMod[a], false, 1);

                // Calculates d0L2, d1L2, d0Pe, d1Pe (no shift)
                indexes[i][a][0] = calculateL2(false, 0);
                indexes[i][a][1] = calculateL2(false, 1);
                indexes[i][a][2] = calculatePearson(false, 0);
                indexes[i][a][3] = calculatePearson(false, 1);

            }

            if (calculateShift == true &&
                possibleToCalculateShifts[i][a] == true) {

                // indexValues will contain the values for an index calculated
                // for each value of shiftAmounts[i][a].
                // indexValues[i][0] = d0L2, indexValues[i][1] = d1L2,
                // indexValues[i][2] = d0Pe, indexValues[i][3] = d1Pe
                auto indexValues =
                        vector<vector<double>>(4,
                                               vector<double>(shiftAmounts[i][a].size(), -1));

                auto sumsOfIndexValues =
                        vector<double>(shiftAmounts[i][a].size(), -1);

                // Calculates d0L2, d1L2, d0Pe and d1Pe for the various shift
                // amounts
                for (int b = 0; b < shiftAmounts[i][a].size(); ++b) {

                    splinesMod[a].calculateShift(shiftAmounts[i][a][b]);

                    // Calculates the norms for the shifted model
                    normD0Mod = calculateNorm(splinesMod[a], true, 0);
                    normD1Mod = calculateNorm(splinesMod[a], true, 1);

                    // Calculates d0L2, d1L2, d0Pe, d1Pe (with shift)
                    indexValues[0][b] = calculateL2(true, 0);
                    indexValues[1][b] = calculateL2(true, 1);
                    indexValues[2][b] = calculatePearson(true, 0);
                    indexValues[3][b] = calculatePearson(true, 1);

                } // End of the for cycle for the shifts

                for (int b = 0; b < shiftAmounts[i][a].size(); ++b) {

                    bool possibleToCalculateSum = true;

                    for (int c = 0; c < 4; ++c)
                        if (indexValues[c][b] == -1)
                            possibleToCalculateSum = false;

                    if (possibleToCalculateSum == true)
                        sumsOfIndexValues[b] =
                                indexValues[0][b] + indexValues[1][b] +
                                indexValues[2][b] + indexValues[3][b];

                }

                int indexOfMax = findIndexOfMax(sumsOfIndexValues);

                maxIndexes[i][a] = indexOfMax;

                if (useSumOfIndexesForAlignment == true)
                    if (indexOfMax != -1) {

                        indexes_shift[i][a][0] = indexValues[0][indexOfMax];
                        indexes_shift[i][a][1] = indexValues[1][indexOfMax];
                        indexes_shift[i][a][2] = indexValues[2][indexOfMax];
                        indexes_shift[i][a][3] = indexValues[3][indexOfMax];

                        double ratio = abs(shiftAmounts[i][a][indexOfMax]) /
                                       splinesExp[i].xRange;

                        if (ratio > 1.) ratio = 1.;

                        shifts[i][a] = 1. - ratio;

                    }

                if (useSumOfIndexesForAlignment == false)
                    for (int b = 0; b < 4; ++b) {
                        int c = findIndexOfMax(indexValues[b]);
                        if (c != -1) {
                            indexes_shift[i][a][b] = indexValues[b][c];

                            double ratio = abs(shiftAmounts[i][a][c]) /
                                           splinesExp[i].xRange;

                            if (ratio > 1.) ratio = 1.;

                            shifts_allDifferent[i][a][b] = 1. - ratio;
                        }
                    }

            } // End of if(calculateShift == true &&
            //           possibleToCalculateShifts[i][a] == true)

        } // End of if (splinesMod[a].possibleToCalculateSpline == true)
    } // End of the for cycle for each model

}


void Indexes::calculateScores() {

    auto matrix_scores = vector<vector<double>>(splinesExp.size(),
                                                vector<double>(numberOfModels, -1));
    auto matrix_scores_shift = vector<vector<double>>(splinesExp.size(),
                                                      vector<double>(numberOfModels, -1));

    // Calculates matrix_scores
    for (int i = 0; i < splinesExp.size(); ++i)
        for (int a = 0; a < numberOfModels; ++a) {

            bool scoreIsNA = false;

            for (int b = 0; b < 4; ++b)
                if (indexes[i][a][b] == -1)
                    scoreIsNA = true;

            for (int b = 0; b < 4; ++b)
                if (isnan(indexes[i][a][b]) == true)
                    scoreIsNA = true;

            if (scoreIsNA == false)
                matrix_scores[i][a] = (indexes[i][a][0] + indexes[i][a][1] +
                                       indexes[i][a][2] + indexes[i][a][3]) / 4.;

        }

    // Calculates matrix_scores_shift, without considering the shift, in order
    // to select the best spline from splinesExp
    if (calculateShift == true)
        for (int i = 0; i < splinesExp.size(); ++i)
            for (int a = 0; a < numberOfModels; ++a) {

                bool scoreIsNA = false;

                for (int b = 0; b < 4; ++b)
                    if (indexes_shift[i][a][b] == -1)
                        scoreIsNA = true;

                for (int b = 0; b < 4; ++b)
                    if (isnan(indexes_shift[i][a][b]) == true)
                        scoreIsNA = true;

                if (scoreIsNA == false)
                    matrix_scores_shift[i][a] =
                            (indexes_shift[i][a][0] + indexes_shift[i][a][1] +
                             indexes_shift[i][a][2] + indexes_shift[i][a][3]) / 4.;

            }

    // Selects the best spline of splinesExp

    auto meansOfScores = vector<double>(splinesExp.size(), -1);

    if (calculateShift == false)
        for (int i = 0; i < splinesExp.size(); ++i) {
            double sumOfScores = 0;
            int counter = 0;
            for (int a = 0; a < numberOfModels; ++a)
                if (matrix_scores[i][a] != -1) {
                    sumOfScores += matrix_scores[i][a];
                    ++counter;
                }
            if (counter != 0)
                meansOfScores[i] = sumOfScores / (double) counter;
        }

    if (calculateShift == true)
        for (int i = 0; i < splinesExp.size(); ++i) {
            double sumOfScores = 0;
            int counter = 0;
            for (int a = 0; a < numberOfModels; ++a)
                if (matrix_scores_shift[i][a] != -1) {
                    sumOfScores += matrix_scores_shift[i][a];
                    ++counter;
                }
            if (counter != 0)
                meansOfScores[i] = sumOfScores / (double) counter;
        }

    indexBestSplineExp = 0;
    double maxMeanOfScores = -1.;
    for (int i = 0; i < splinesExp.size(); ++i)
        if (meansOfScores[i] != -1)
            if (meansOfScores[i] > maxMeanOfScores) {
                maxMeanOfScores = meansOfScores[i];
                indexBestSplineExp = i;
            }

    bootstrapScores[bootstrapIndex] = matrix_scores[indexBestSplineExp];

    bootstrapScores_shift[bootstrapIndex] =
            matrix_scores_shift[indexBestSplineExp];
    if (calculateShift == true)
        for (int a = 0; a < numberOfModels; ++a)
            if (matrix_scores_shift[indexBestSplineExp][a] != -1)
                useSumOfIndexesForAlignment == true ?

                        bootstrapScores_shift[bootstrapIndex][a] =
                                (matrix_scores_shift[indexBestSplineExp][a] * 4. +
                                 (shifts[indexBestSplineExp][a] * 2.)) / 6. :

                        bootstrapScores_shift[bootstrapIndex][a] =
                                (matrix_scores_shift[indexBestSplineExp][a] * 4. +
                                 (shifts_allDifferent[indexBestSplineExp][a][0] +
                                  shifts_allDifferent[indexBestSplineExp][a][1] +
                                  shifts_allDifferent[indexBestSplineExp][a][2] +
                                  shifts_allDifferent[indexBestSplineExp][a][3]) / 2.) / 6.;

    // Rows: models, Columns: d0L2, d1L2, d0Pe, d1Pe, shift (or shifts if all
    // different)
    auto idx = vector<vector<double>>(numberOfModels);

    for (int a = 0; a < numberOfModels; ++a) {

        if (calculateShift == true)
            idx[a] = indexes_shift[indexBestSplineExp][a];
        else
            idx[a] = indexes[indexBestSplineExp][a];

        if (useSumOfIndexesForAlignment == true)
            idx[a].push_back(shifts[indexBestSplineExp][a]);
        else
            for (int b = 0; b < 4; ++b)
                idx[a].push_back(shifts_allDifferent[indexBestSplineExp][a][b]);
    }

    allIndexes.push_back(idx);

}


int Indexes::findIndexOfMax(const vector<double> &sums) {

    int positionMax = -1;
    double max = sums[0];
    int k = 0;
    if (sums.size() % 2 != 0) ++k;

    for (int a = 0; a < (sums.size() + k) / 2; ++a)
        if (sums[a] != -1)
            if (sums[a] >= max) {
                positionMax = a;
                max = sums[a];
            }
    for (int a = (sums.size() + k) / 2; a < sums.size(); ++a)
        if (sums[a] != -1)
            if (sums[a] > max) {
                positionMax = a;
                max = sums[a];
            }

    // Returns -1 if 'sums' contains only -1
    if (positionMax == -1) return -1;

    return positionMax;

}


double Indexes::calculateL2(bool shift, int d0d1) {

    int k = indexOfModelBeingCalculated;

    vector<double> y;

    if (shift == false) {

        if (d0d1 == 0)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y.push_back(yD0[a] - splinesMod[k].D0(xPowers[a]));
        if (d0d1 == 1)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y.push_back(yD1[a] - splinesMod[k].D1(xPowers[a]));

    }

    if (shift == true) {

        if (d0d1 == 0)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y.push_back(yD0[a] - splinesMod[k].D0Shift(xPowers[a]));
        if (d0d1 == 1)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y.push_back(yD1[a] - splinesMod[k].D1Shift(xPowers[a]));

    }

    double baseAlpha;
    double baseBeta;
    double integral = 0;

    baseAlpha = y[0] * y[0];
    for (int c = 1; c < numberOfTrapezoids + 1; ++c) {
        baseBeta = y[c] * y[c];
        integral += (baseAlpha + baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return 1. / (1. + sqrt(integral) / splinesExp[splineExpIndex].xRange);

}


double Indexes::calculatePearson(bool shift, int d0d1) {

    // If the norm of the experimental data is 0, the Pearson index can't be
    // calculated
    if (normD0Exp == 0 && d0d1 == 0) return -1;
    if (normD1Exp == 0 && d0d1 == 1) return -1;

    // If the norm of the model is 0, the Pearson index can't be calculated
    if (normD0Mod == 0 && d0d1 == 0) return -1;
    if (normD1Mod == 0 && d0d1 == 1) return -1;

    int k = indexOfModelBeingCalculated;

    //vector<double> y;

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

    vector<double> y2;

    if (shift == false) {

        if (d0d1 == 0)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                double g = splinesMod[k].D0(xPowers[a]);
                double y = yD0[a] / normD0Exp - g / normD0Mod;
                y2.push_back(y * y);
            }
        if (d0d1 == 1)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                double g = splinesMod[k].D1(xPowers[a]);
                double y = yD1[a] / normD1Exp - g / normD1Mod;
                y2.push_back(y * y);
            }

    }

    if (shift == true) {

        if (d0d1 == 0)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                double g = splinesMod[k].D0Shift(xPowers[a]);
                double y = yD0[a] / normD0Exp - g / normD0Mod;
                y2.push_back(y * y);
            }
        if (d0d1 == 1)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                double g = splinesMod[k].D1Shift(xPowers[a]);
                double y = yD1[a] / normD1Exp - g / normD1Mod;
                y2.push_back(y * y);
            }

    }

    double baseAlpha;
    double baseBeta;
    double result = 0;

    baseAlpha = y2[0];
    for (int c = 1; c < numberOfTrapezoids + 1; ++c) {
        baseBeta = y2[c];
        result += (baseAlpha + baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return 1. - sqrt(result) / 2.;

}


double Indexes::calculateNorm(Spline &spline, bool shift, int d0d1) {

    vector<double> y2;

    if (spline.splineType != 1 /*Model*/) {

        if (d0d1 == 0)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y2.push_back(yD0[a] * yD0[a]);
        if (d0d1 == 1)
            for (int a = 0; a < numberOfTrapezoids + 1; ++a)
                y2.push_back(yD1[a] * yD1[a]);

    } else {

        if (shift == false) {

            if (d0d1 == 0)
                for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                    double y = spline.D0(xPowers[a]);
                    y2.push_back(y * y);
                }
            if (d0d1 == 1)
                for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                    double y = spline.D1(xPowers[a]);
                    y2.push_back(y * y);
                }

        }

        if (shift == true) {

            if (d0d1 == 0)
                for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                    double y = spline.D0Shift(xPowers[a]);
                    y2.push_back(y * y);
                }
            if (d0d1 == 1)
                for (int a = 0; a < numberOfTrapezoids + 1; ++a) {
                    double y = spline.D1Shift(xPowers[a]);
                    y2.push_back(y * y);
                }

        }

    }

    double baseAlpha;
    double baseBeta;
    double result = 0;

    baseAlpha = y2[0];
    for (int c = 1; c < numberOfTrapezoids + 1; ++c) {
        baseBeta = y2[c];
        if (d0d1 == 0)
            result += (baseAlpha + baseBeta) * halfHeight;
        if (d0d1 == 1)
            result += (baseAlpha + baseBeta) * halfHeight;
        baseAlpha = baseBeta;
    }

    return sqrt(result);

}


void Indexes::sortScores(vector<double> scoresIn,
                         vector<double> &scoresOut,
                         vector<double> stdDevIn,
                         vector<double> &stdDevOut,
                         vector<string> &scoresOutNames,
                         vector<int> &scoresOutLineTypes,
                         vector<int> &scoresOutColorIndexes) {

    vector<string> names = modelNames;

    vector<int> lineTypes;
    for (int a = 0; a < numberOfModels; ++a)
        lineTypes.push_back(a % 3 + 1);

    vector<int> colorIndexes;
    for (int a = 1; a < numberOfModels + 1; ++a)
        colorIndexes.push_back(a);

    while (scoresIn.size() > 0) {

        double max = -1;
        int indexMax;

        for (int a = 0; a < scoresIn.size(); ++a)
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
            scoresIn.erase(scoresIn.begin() + indexMax);
            stdDevIn.erase(stdDevIn.begin() + indexMax);
            names.erase(names.begin() + indexMax);
            lineTypes.erase(lineTypes.begin() + indexMax);
            colorIndexes.erase(colorIndexes.begin() + indexMax);
        } else {
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


void Indexes::calculateShiftAmounts(int modelNumber) {

    int i = splineExpIndex;
    int a = modelNumber;

    double expLeft = splinesExp[i].xMinNoAsymptotes;
    double expRight = splinesExp[i].xMaxNoAsymptotes;

    vector<double> acceptedShifts;

    bool possibleToLineUpMaxima = false;

    if (lineUpMaxima == true) {

        double shift;

        if (splinesExp[i].hasOneMaximum == true)
            if (splinesMod[a].hasOneMaximum == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp[i].xOfMax - splinesMod[a].xOfMax;
            }

        if (splinesExp[i].hasOneMaximum == true)
            if (splinesMod[a].hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                double xOfClosestPoint = splinesMod[a].xOfFirstMax;
                if (fabs(splinesExp[i].xOfMax - splinesMod[a].xOfSecondMax) <
                    fabs(splinesExp[i].xOfMax - splinesMod[a].xOfFirstMax))
                    xOfClosestPoint = splinesMod[a].xOfSecondMax;
                shift = splinesExp[i].xOfMax - xOfClosestPoint;
            }

        if (splinesExp[i].hasTwoMaxima == true)
            if (splinesMod[a].hasOneMaximum == true) {
                possibleToLineUpMaxima = true;
                double xOfClosestPoint = splinesExp[i].xOfFirstMax;
                if (fabs(splinesMod[a].xOfMax - splinesExp[i].xOfSecondMax) <
                    fabs(splinesMod[a].xOfMax - splinesExp[i].xOfFirstMax))
                    xOfClosestPoint = splinesExp[i].xOfSecondMax;
                shift = xOfClosestPoint - splinesMod[a].xOfMax;
            }

        if (splinesExp[i].hasTwoMaxima == true)
            if (splinesMod[a].hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp[i].xOfFirstMax - splinesMod[a].xOfFirstMax;
            }

        if (splinesExp[i].hasTwoMaxima == true)
            if (splinesMod[a].hasFirstOfTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp[i].xOfFirstMax - splinesMod[a].xOfMax;
            }

        if (splinesExp[i].hasFirstOfTwoMaxima == true)
            if (splinesMod[a].hasFirstOfTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp[i].xOfMax - splinesMod[a].xOfMax;
            }

        if (splinesExp[i].hasFirstOfTwoMaxima == true)
            if (splinesMod[a].hasTwoMaxima == true) {
                possibleToLineUpMaxima = true;
                shift = splinesExp[i].xOfMax - splinesMod[a].xOfFirstMax;
            }

        if (possibleToLineUpMaxima == true) {

            acceptedShifts.push_back(shift - shiftAroundMaximum);
            double limit = shift + shiftAroundMaximum;
            while (acceptedShifts.back() < limit)
                acceptedShifts.push_back(
                        acceptedShifts.back() + distanceBetweenShiftedPoints);

        }

        if (acceptedShifts.size() == 0)
            possibleToLineUpMaxima = false;

    }

    if (possibleToLineUpMaxima == false) {

        vector<double> shiftsVector;

        double modLeft = splinesMod[a].xMinNoAsymptotes;
        double modRight = splinesMod[a].xMaxNoAsymptotes;

        shiftsVector.push_back(expLeft - modRight);
        while (modLeft + shiftsVector.back() < expRight)
            shiftsVector.push_back(
                    shiftsVector.back() + distanceBetweenShiftedPoints);

        // If the model has either one or two maxima, the shifts which would
        // make its leftmost original abscissa < 0 are accepted only if at least
        // one maximum remains in the experimental data range after the shift

        vector<double> relevantShifts;

        if (shiftsVector.size() > 0)
            for (int b = 0; b < shiftsVector.size(); ++b) {

                if (splinesMod[a].originalAbscissae[0] + shiftsVector[b] >= 0)
                    relevantShifts.push_back(shiftsVector[b]);

                if (splinesMod[a].originalAbscissae[0] + shiftsVector[b] < 0) {

                    if (splinesMod[a].hasOneMaximum == true) {
                        double max = splinesMod[a].xOfMax + shiftsVector[b];
                        if (max >= expLeft && max <= expRight)
                            relevantShifts.push_back(shiftsVector[b]);
                    }

                    if (splinesMod[a].hasTwoMaxima == true) {
                        bool atLeastOneMaximumWithinExpRange = false;
                        double maxLeft =
                                splinesMod[a].xOfFirstMax + shiftsVector[b];
                        double maxRight =
                                splinesMod[a].xOfSecondMax + shiftsVector[b];
                        if (maxLeft >= expLeft && maxLeft <= expRight)
                            atLeastOneMaximumWithinExpRange = true;
                        if (maxRight >= expLeft && maxRight <= expRight)
                            atLeastOneMaximumWithinExpRange = true;
                        if (atLeastOneMaximumWithinExpRange == true)
                            relevantShifts.push_back(shiftsVector[b]);
                    }

                    if (splinesMod[a].hasOneMaximum == false &&
                        splinesMod[a].hasTwoMaxima == false)
                        relevantShifts.push_back(shiftsVector[b]);

                }

            }

        // Checks whether the calculated shift amounts are compatible with the
        // limit on model extrapolation
        if (relevantShifts.size() > 0)
            for (int b = 0; b < relevantShifts.size(); ++b) {
                double extrapolatedLeft =
                        splinesMod[a].originalAbscissae[0] +
                        relevantShifts[b] - expLeft;
                if (extrapolatedLeft < 0)
                    extrapolatedLeft = 0;
                double extrapolatedRight =
                        expRight - relevantShifts[b] -
                        splinesMod[a].originalAbscissae.back();
                if (extrapolatedRight < 0)
                    extrapolatedRight = 0;
                if (extrapolatedLeft + extrapolatedRight <= maxExtrapolatedLength)
                    acceptedShifts.push_back(relevantShifts[b]);
            }

    } // End of if (possibleToLineUpMaxima == false)

    if (acceptedShifts.size() > 0) {
        possibleToCalculateShifts[i][a] = true;
        shiftAmounts[i][a] = acceptedShifts;
    }

}


void Indexes::statisticalAnalysis(const vector<vector<double>> &bootstrap,
                                  int i, /*index of the model in modelNames*/
                                  double &mean,
                                  double &stdDeviation) {

    double n = (double) numberOfBootstrapVariations;

    // Calculates the mean
    double sum = 0;
    for (int a = 0; a < numberOfBootstrapVariations; ++a)
        sum += bootstrap[a][i];
    mean = sum / n;

    // Calculates the standard deviation of the mean
    double sumOfSquares = 0;
    for (int a = 0; a < numberOfBootstrapVariations; ++a) {
        double difference = bootstrap[a][i] - mean;
        sumOfSquares += difference * difference;
    }
    stdDeviation = sqrt(sumOfSquares / (n - 1.)) / sqrt(n);

}


void Indexes::saveIndexesToCsvFile() {

    string fileNameWithPath = "./Indexes/" + folderName + "_" + fileName + ".csv";

    ofstream file(fileNameWithPath.c_str());

    for (int a = 0; a < numberOfBootstrapVariations; ++a) {

        file << ",d0L2,d1L2,d0Pe,d1Pe,";
        useSumOfIndexesForAlignment == true ?
        file << "shift\n" :
        file << "shift_d0L2,shift_d1L2,shift_d0Pe,shift_d1Pe\n";

        for (int b = 0; b < numberOfModels; ++b) {
            file << modelNames[b] << ",";
            for (int c = 0; c < allIndexes[a][b].size() - 1; ++c)
                file << allIndexes[a][b][c] << ",";
            file << allIndexes[a][b][allIndexes[a][b].size() - 1] << "\n";
        }

        if (a != numberOfBootstrapVariations - 1)
            file << "\n";

    }

    file.close();

}
