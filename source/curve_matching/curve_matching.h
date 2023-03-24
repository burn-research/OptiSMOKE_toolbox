#ifndef CURVE_MATCHING_H
#define CURVE_MATCHING_H

#include "GlobalVariables.h"
#include "BasisFunction.h"
#include "Spline.h"
#include "Indexes.h"
#include "Results.h"

std::vector<double> curveMatching(int numberOfBootstrapVariations_par,
								vector<double> exp_data_x,
								vector<double> exp_data_y,
								vector<double> model_data_x,
								vector<double> model_data_y,
								vector<double> error_data) 
{
    numberOfBootstrapVariations = numberOfBootstrapVariations_par;

    Indexes solver;
    solver.solvePython(false,
                       exp_data_x, 
					   exp_data_y,
                       model_data_x, 
					   model_data_y,
                       error_data);

    std::vector<Indexes> ind_vect;
    ind_vect.push_back(solver);

    Results r;
    r.calculate(ind_vect);

    std::vector<double> vector_scores;
    std::vector<double> vector_errors;
    std::vector<std::string> vector_names;

    vector_scores = r.KintSorted;
    vector_errors = r.KintStdDevSorted;
	vector_names = r.modelNamesSorted;
	// vector<double> results_scores;
	// vector<double> results_errors;

    // for (int i = 0; i < vector_scores.size(); i++) {
    //     results_scores.push_back(vector_scores[i]);
    //     results_errors.push_back(vector_errors[i]);
    // }
	// return results_scores;
    return ind_vect[0].allIndexes[0][0];
}
#endif // CURVE_MATCHING_H
