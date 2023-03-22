//
// Created by Edoardo Ramalli on 24/09/2020.
//

#ifndef CURVEMATCHINGPYTHON_CURVEMATCHINGPYTHON_H
#define CURVEMATCHINGPYTHON_CURVEMATCHINGPYTHON_H
extern "C"
void curveMatching(int m_par,
                   int g_par,
                   double lambdaSearchInterval_par,
                   int numberOfStepsLambda_par,
                   double fractionOfOrdinateRangeForAsymptoteIdentification_par,
                   double fractionOfOrdinateRangeForMaximumIdentification_par,
                   double fractionOfExpHeightForZeroLineIdentification_par,
                   double fractionOfExpRangeForModelExtrapolation_par,
                   double fractionOfExpRangeForMinShift_par,
                   double fractionOfExpRangeForShiftAroundMaximum_par,
                   bool possibleNegativeOrdinates_par,
                   int numberOfBootstrapVariations_par,
                   bool lineUpMaxima_par,
                   bool useSumOfIndexesForAlignment_par,
                   int numberOfTrapezoids_par,
                   bool CalculateShift,
//                   double *results_indexes,
                   double *results_scores,
                   double *results_errors,
                   double *exp_data_x,
                   double *exp_data_y,
                   int exp_data_size,
                   double *model_data_x,
                   double *model_data_y,
                   int model_data_size,
                   double *error_data) {

    /* Settings */

    m = m_par;
    g = g_par;
    lambdaSearchInterval = lambdaSearchInterval_par;
    numberOfStepsLambda = numberOfStepsLambda_par;
    fractionOfOrdinateRangeForAsymptoteIdentification = fractionOfOrdinateRangeForAsymptoteIdentification_par;
    fractionOfOrdinateRangeForMaximumIdentification = fractionOfOrdinateRangeForMaximumIdentification_par;
    fractionOfExpHeightForZeroLineIdentification = fractionOfExpHeightForZeroLineIdentification_par;
    fractionOfExpRangeForModelExtrapolation = fractionOfExpRangeForModelExtrapolation_par;
    fractionOfExpRangeForMinShift = fractionOfExpRangeForMinShift_par;
    fractionOfExpRangeForShiftAroundMaximum = fractionOfExpRangeForShiftAroundMaximum_par;
    possibleNegativeOrdinates = possibleNegativeOrdinates_par;
    numberOfBootstrapVariations = numberOfBootstrapVariations_par;
    lineUpMaxima = lineUpMaxima_par;
    useSumOfIndexesForAlignment = useSumOfIndexesForAlignment_par;
    numberOfTrapezoids = numberOfTrapezoids_par;



    Indexes solver;
    solver.solvePython(CalculateShift,
//                       results_indexes,
                       exp_data_x, exp_data_y, exp_data_size,
                       model_data_x, model_data_y, model_data_size,
                       error_data);

    vector<Indexes> ind_vect;
    ind_vect.push_back(solver);

    Results r;
    r.calculate(ind_vect);

    vector<double> vector_scores;
    vector<double> vector_errors;
    vector<string> vector_names;


    if (CalculateShift) {
        vector_scores = r.KintSorted_shift;
        vector_errors = r.KintStdDevSorted_shift;
        vector_names = r.modelNamesSorted_shift;
    } else {
        vector_scores = r.KintSorted;
        vector_errors = r.KintStdDevSorted;
        vector_names = r.modelNamesSorted;
    }

    for (int i = 0; i < vector_scores.size(); i++) {
        results_scores[i] = vector_scores[i];
        results_errors[i] = vector_errors[i];
    }

}
#endif //CURVEMATCHINGPYTHON_CURVEMATCHINGPYTHON_H
