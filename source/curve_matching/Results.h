class Results {

public:

    /* List of the names of the files contained in the input folder, minus the
    _exp.txt and _mod.txt appendixes. For files 'Example_exp.txt' and
    'Example_mod.txt', nameList would contain the string 'Example' */
    vector<string> nameList;

    /* Number of files in the input folder */
    int numberOfFiles;

    /* Names of the models in the input files */
    vector<string> modelNames;

    /* Number of models */
    int numberOfModels;

    /* Score matrix, without shifts */
    vector<vector<double>> Kmatrix;

    /* Score matrix, with shifts */
    vector<vector<double>> Kmatrix_shift;

    /* Matrix with the values of the standard deviation, without shifts */
    vector<vector<double>> stdDevMatrix;

    /* Matrix with the values of the standard deviation, with shifts */
    vector<vector<double>> stdDevMatrix_shift;

    /* Score vector, without shifts */
    vector<double> Kint;

    /* Score vector, with shifts */
    vector<double> Kint_shift;

    /* Standard deviations of Kint */
    vector<double> KintStdDev;

    /* Standard deviations of Kint_shift */
    vector<double> KintStdDev_shift;

    /* Score vector, without shifts, sorted from largest to smallest */
    vector<double> KintSorted;

    /* Standard deviations of Kint, sorted as the elements of KintSorted */
    vector<double> KintStdDevSorted;

    /* Model names, sorted as the elements of KintSorted */
    vector<string> modelNamesSorted;

    /* Score vector, with shifts, sorted from smallest to largest */
    vector<double> KintSorted_shift;

    /* Standard deviations of Kint_shift, sorted as the elements of
    KintSorted_shift */
    vector<double> KintStdDevSorted_shift;

    /* Model names, sorted as the elements of KintSorted_shift */
    vector<string> modelNamesSorted_shift;

    ////////////////////////////////////////////////////////////////////////////

    /* Calculates Kmatrix, Kmatrix_shift, Kint and Kint_shift */
    void calculate(const vector<Indexes>& indexes);

////////////////////////////////////////////////////////////////////////////////

private:

    /* Name of the folder containing the input data */
    string folderName;

    ////////////////////////////////////////////////////////////////////////////

    /* Sorts KintIn from the highest to the lowest value, and sorts stdDevIn
    accordingly. Saves the results to KintOut and stdDevOut. namesOut is sorted
    in the same way */
    void sortVectors(vector<double> KintIn,
                     vector<double>& KintOut,
                     vector<double> stdDevIn,
                     vector<double>& stdDevOut,
                     vector<string>& namesOut);

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void Results::calculate(const vector<Indexes>& indexes) {

    numberOfFiles = indexes.size();

    for (int a=0; a<numberOfFiles; ++a)
        nameList.push_back(indexes[a].fileName);

    // Finds the names of the models in the input files
    modelNames = indexes[0].modelNames;
    if (indexes.size() > 1)
        for (int a=0; a<numberOfFiles; ++a)
            if (indexes[a].modelNames != indexes[0].modelNames)
                for (int b=0; b<indexes[a].numberOfModels; ++b) {
                    bool newName = true;
                    for (int c=0; c<modelNames.size(); ++c)
                        if (indexes[a].modelNames[b] == modelNames[c]) {
                            newName = false;
                            break;
                        }
                    if (newName == true)
                        modelNames.push_back(indexes[a].modelNames[b]);
                }

    numberOfModels = modelNames.size();

    // Initializes Kmatrix, Kmatrix_shift, stdDevMatrix and stdDevMatrix_shift
    Kmatrix = Kmatrix_shift = stdDevMatrix = stdDevMatrix_shift =
            vector<vector<double>>(numberOfModels,vector<double>(numberOfFiles,-1));

    // Fills Kmatrix with the scores and stdDevMatrix with the corresponding
    // standard deviations
    for (int a=0; a<numberOfFiles; ++a)
        for (int b=0; b<numberOfModels; ++b)
            for (int c=0; c<indexes[a].numberOfModels; ++c)
                if (modelNames[b] == indexes[a].modelNames[c]) {
                    Kmatrix[b][a] = indexes[a].scores[c];
                    stdDevMatrix[b][a] = indexes[a].stdDev[c];
                }

    // Fills Kmatrix_shift with the scores and stdDevMatrix_shift with the
    // corresponding standard deviations
    for (int a=0; a<numberOfFiles; ++a)
        for (int b=0; b<numberOfModels; ++b)
            for (int c=0; c<indexes[a].numberOfModels; ++c)
                if (modelNames[b] == indexes[a].modelNames[c]) {
                    Kmatrix_shift[b][a] = indexes[a].scores_shift[c];
                    stdDevMatrix_shift[b][a] = indexes[a].stdDev_shift[c];
                }

    // Initializes Kint, Kint_shift, KintStdDev and KintStdDev_shift
    Kint = Kint_shift = KintStdDev = KintStdDev_shift =
            vector<double>(numberOfModels,-1);

    // Calculates Kint and KintStdDev
    for (int a=0; a<numberOfModels; ++a) {
        double numerator = 0;
        double numeratorStdDev = 0;
        double denominator = 0;
        for (int b=0; b<numberOfFiles; ++b)
            if (Kmatrix[a][b] != -1) {
                numerator += Kmatrix[a][b];
                numeratorStdDev += stdDevMatrix[a][b];
                ++denominator;
            }
        if (denominator > 0) {
            Kint[a] = numerator / denominator;
            KintStdDev[a] = numeratorStdDev / denominator;
        }
    }

    // Calculates Kint_shift and KintStdDev_shift
    for (int a=0; a<numberOfModels; ++a) {
        double numerator = 0;
        double numeratorStdDev = 0;
        double denominator = 0;
        for (int b=0; b<numberOfFiles; ++b)
            if (Kmatrix_shift[a][b] != -1) {
                numerator += Kmatrix_shift[a][b];
                numeratorStdDev += stdDevMatrix_shift[a][b];
                ++denominator;
            }
        if (denominator > 0) {
            Kint_shift[a] = numerator / denominator;
            KintStdDev_shift[a] = numeratorStdDev / denominator;
        }
    }


    sortVectors(Kint,
                KintSorted,
                KintStdDev,
                KintStdDevSorted,
                modelNamesSorted);

    sortVectors(Kint_shift,
                KintSorted_shift,
                KintStdDev_shift,
                KintStdDevSorted_shift,
                modelNamesSorted_shift);

}



void Results::sortVectors(vector<double> KintIn,
                          vector<double>& KintOut,
                          vector<double> stdDevIn,
                          vector<double>& stdDevOut,
                          vector<string>& namesOut) {

    vector<string> names = modelNames;

    while (KintIn.size() > 0) {

        double max = -1;
        int indexMax;

        for (int a=0; a<KintIn.size(); ++a)
            if (KintIn[a] > max) {
                max = KintIn[a];
                indexMax = a;
            }

        if (max != -1) { // If at least one score is not NA
            KintOut.push_back(max);
            stdDevOut.push_back(stdDevIn[indexMax]);
            namesOut.push_back(names[indexMax]);
            KintIn.erase(KintIn.begin()+indexMax);
            stdDevIn.erase(stdDevIn.begin()+indexMax);
            names.erase(names.begin()+indexMax);
        }
        else {
            KintOut.push_back(KintIn[0]);
            stdDevOut.push_back(stdDevIn[0]);
            namesOut.push_back(names[0]);
            KintIn.erase(KintIn.begin());
            stdDevIn.erase(stdDevIn.begin());
            names.erase(names.begin());
        }

    }

}