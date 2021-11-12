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

#include "GlobalVariables.h"
#include "BasisFunction.h"
class Spline {

public:

    /* Type of spline. 0: Experimental data;  1: Model;  2: Error spline */
    int splineType;

    /* Abscissae of the spline */
    std::vector<double> abscissae;

    /* Ordinates of the spline */
    std::vector<double> ordinates;

    /* Number of data points */
    int n;

    /* Specifies whether there are enough data points to calculate the spline,
    or whether the spline can be considered a flat line with ordinate = 0 when
    compared to the experimental data */
    bool possibleToCalculateSpline;

    /* Real knots of the spline */
    std::vector<double> knots;

    /* Number of real knots of the spline */
    int numberOfKnots;

    /* Number of polynomials of the spline */
    int numberOfPolynomials;

    /* Coefficients of the polynomials of the spline, excluding those
    corresponding to coincident knots at the end points. coeffD0[i][j] refers to
    polynomial i and the coefficient of x^j */
    std::vector<std::vector<double>> coeffD0;

    /* Coefficients of the first derivatives of the polynomials of the spline,
    excluding those corresponding to coincident knots at the end points.
    coeffD1[i][j] refers to polynomial i and the coefficient of x^j */
    std::vector<std::vector<double>> coeffD1;

    /* Coefficients of the second derivatives of the polynomials of the spline,
    excluding those corresponding to coincident knots at the end points.
    coeffD2[i][j] refers to polynomial i and the coefficient of x^j */
    std::vector<std::vector<double>> coeffD2;

    /* Abscissae of the spline, as initially obtained from the input file */
    std::vector<double> originalAbscissae;

    /* Ordinates of the spline, as initially obtained from the input file */
    std::vector<double> originalOrdinates;

    /* Distance between the biggest and the smallest abscissae of the spline */
    double xRange;

    /* Maximum abscissa of the spline if the aymptotes were removed */
    double xMaxNoAsymptotes;

    /* Minimum abscissa of the spline if the aymptotes were removed */
    double xMinNoAsymptotes;

    /* Maximum ordinate of the spline */
    double yD0Max;

    /* Minimum ordinate of the spline */
    double yD0Min;

    /* Distance between the biggest and the smallest ordinates of the spline */
    double yD0Range;

    /* Maximum ordinate of the first derivative of the spline */
    double yD1Max;

    /* Minimum ordinate of the first derivative of the spline */
    double yD1Min;

    /* Maximum ordinate of the spline before the addition of any segments */
    double yD0MaxOriginal;

    /* Minimum ordinate of the spline before the addition of any segments */
    double yD0MinOriginal;

    /* Maximum ordinate of the first derivative of the spline before the
    addition of any segments */
    double yD1MaxOriginal;

    /* Minimum ordinate of the first derivative of the spline before the
    addition of any segments */
    double yD1MinOriginal;

    /* Biggest between the absolute value of the maximum and the absolute value
    of the minimun of the first derivative of the spline */
    double yD1MaxAbs;

    /* Specifies whether the spline has a single well-defined maximum between
    the extremes */
    bool hasOneMaximum;

    /* Abscissa of the maximum of the spline, if a single well-defined maximum
    is present between the extremes */
    double xOfMax;

    /* Specifies whether the spline has two well-defined maxima between the
    extremes */
    bool hasTwoMaxima;

    /* Abscissa of the first maximum of the spline, if two well-defined maxima
    are present between the extremes */
    double xOfFirstMax;

    /* Abscissa of the second maximum of the spline, if two well-defined maxima
    are present between the extremes */
    double xOfSecondMax;

    /* Specifies whether the spline has a single well-defined maximum between
    the extremes, but would have had two had additional data been provided */
    bool hasFirstOfTwoMaxima;

    ////////////////////////////////////////////////////////////////////////////

    /* Calculates the spline */
    void solve(const std::vector<double>& abscissae,
               const std::vector<double>& ordinates,
               int splineType,
               int numberOfAbscissaeSeparatingConsecutiveKnots,
	       bool print_splines);
    
    bool print_splines;
    /* Removes any horizontal asymptotes on the left and on the right of the
    spline */
    void removeAsymptotes();

    /* Calculates the ordinate of the spline at position x on the x-axis */
    double D0(double x);

    /* Calculates the ordinate of the first derivative of the spline at position
    x on the x-axis */
    double D1(double x);

    /* Calculates the ordinate of the second derivative of the spline at
    position x on the x-axis */
    double D2(double x);

    /* Calculates the ordinate of the spline at position powersOfX[1] on the
    x-axis, using the normalized spline coefficients */
    double D0(const std::vector<double>& powersOfX);

    /* Calculates the ordinate of the first derivative of the spline at position
    powersOfX[1] on the x-axis, using the normalized spline coefficients */
    double D1(const std::vector<double>& powersOfX);

    /* Calculates the ordinate of the spline at position powersOfX[1] on the
    x-axis, using the normalized and shifted spline coefficients */
    double D0Shift(const std::vector<double>& powersOfX);

    /* Calculates the ordinate of the first derivative of the spline at position
    powersOfX[1] on the x-axis, using the normalized and shifted spline
    coefficients */
    double D1Shift(const std::vector<double>& powersOfX);

    /* Adds at the left of the spline a piece of length lengthLeftSide and at
    the right of the spline a piece of length lengthRightSide */
    void extendSpline(double lengthLeftSide, double lengthRightSide);

    /* Calculates coeffD0_shift_normalized, coeffD1_shift_normalized and
    knots_shift */
    void calculateShift(double shift);

    /* Normalizes the coefficients of the spline and of its first derivative, by
    translating them by verticalShift and dividing them by normalizationFactorD0
    and normalizationFactorD1 */
    void normalizeCoefficients(double verticalShift,
                               double normalizationFactorD0,
                               double normalizationFactorD1);

////////////////////////////////////////////////////////////////////////////////

private:

    /* Degrees of freedom of the spline */
    int K;

    /* Degrees of freedom of the spline minus 1 */
    int G;

    /* Smoothing parameter */
    double lambda;

    /* Base 10 logarithm of the smoothing parameter lambda */
    double log10lambda;

    /* Value of log10lambda for which the elements of the matrices FiTFi and R
    have the same order of magnitude, rounded to the nearest 0.5 */
    double log10lambdaForSameOrderOfMagnitude;

    /* Lowest value of the interval for the search of the minimum for
    log10lambda */
    double log10lambdaMin;

    /* Highest value of the interval for the search of the minimum for
    log10lambda */
    double log10lambdaMax;

    /* Spline coefficients for obtaining coeffD0, coeffD1 and coeffD2 */
    std::vector<double> splineCoefficients;

    /* Knots of the spline, including non-real ones at the end points */
    std::vector<double> knotsForCalculations;

    /* Number of polynomials which are asymptotes on the left of the spline */
    double numberOfAsymptotePolynomialsLeft;

    /* Number of polynomials which are asymptotes on the right of the spline */
    double numberOfAsymptotePolynomialsRight;

    /* Powers of the abscissa being considered */
    std::vector<double> powers;

    /* coeffD0, normalized */
    std::vector<std::vector<double>> coeffD0_normalized;

    /* coeffD1, normalized */
    std::vector<std::vector<double>> coeffD1_normalized;

    /* Shift value applied to coeffD0_normalized and coeffD1_normalized to
    obtain coeffD0_shift_normalized and coeffD1_shift_normalized, and to 'knots'
    to obtain knots_shift */
    double shift;

    /* coeffD0, normalized and shifted */
    std::vector<std::vector<double>> coeffD0_shift_normalized;

    /* coeffD1, normalized and shifted */
    std::vector<std::vector<double>> coeffD1_shift_normalized;

    /* Real knots of the spline, shifted */
    std::vector<double> knots_shift;

    ////////////////////////////////////////////////////////////////////////////

    /* Chooses the knots for the spline */
    void chooseKnots(int numberOfAbscissaeSeparatingConsecutiveKnots);

    /* Calculates the coefficients of the polynomials of the spline and the
    coefficients of the first derivative of the polynomials of the spline, for
    the real knots of the spline */
    void calculateCoefficients();

    /* Calculates the slopes at the end points of the spline, as the derivatives
    at the end points of the splines obtained by only fitting either a segment
    on the left side or a segment on the right side of the original spline */
    void calculateSlopeAtEndPoints(double& slopeLeftSide,
                                   double& slopeRightSide);

    /* Looks for horizontal asymptotes */
    void yAndAsymptoteAnalysis();
    /* Calculates the real different roots of the spline or of the first or of
    the second derivative of the spline. Returns a vector with the roots sorted
    from smallest to largest. Returns a size() = 0 vector if there are no real
    roots */
    std::vector<double> calculateRoots(double derivativeOrder);

    /* Replaces any negative segments of the spline with straight lines with
    ordinate 0 */
    void removeNegativeSegments();

    /* After changing the polynomials and/or the knots of the spline, updates
    the affected variables */
    void updateVariables();

    /* Finds whether the spline has either one or two well-defined maxima
    between the extremes */
    void findMaximaBetweenExtremes();

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void Spline::solve(const std::vector<double>& Abscissae,
                   const std::vector<double>& Ordinates,
                   int SplineType,
                   int numberOfAbscissaeSeparatingConsecutiveKnots,
		   bool print) {
    print_splines = print;
    
    abscissae = Abscissae;
    ordinates = Ordinates;
    splineType = SplineType;

    n = abscissae.size();

    if (originalAbscissae.size() == 0) {
        originalAbscissae = abscissae;
        originalOrdinates = ordinates;
    }
    //for (int i=0; i<originalAbscissae.size(); ++i)
    //{
    //	std::cout<<"originalAbscissae["<<i<<"] is "<<originalAbscissae[i]<<std::endl;
    //}
    //std::cout<<"originalAbscissae.back() is "<<originalAbscissae.back()<<std::endl;	
    possibleToCalculateSpline = abscissae.size() > 1 ? true : false;
    //std::cout<<"possibleToCalculateSpline is  "<<possibleToCalculateSpline<<std::endl;
    //std::cout <<"It's possible to calculate this spline? "<<possibleToCalculateSpline << std::endl;
    if (!possibleToCalculateSpline)
        return;
    //std::cout<<"it is possible to calculate the spline"<<std::endl;
    this->chooseKnots(numberOfAbscissaeSeparatingConsecutiveKnots);
    //std::cout<<"Knots chosen"<<std::endl;
    this->calculateCoefficients();
    //std::cout<<"Coefficients: done"<<std::endl;
    this->yAndAsymptoteAnalysis();
    //std::cout<<"Asymtotes: done"<<std::endl;
    if (!possibleNegativeOrdinates)
        this->removeNegativeSegments();
    //std::cout<<"Negative Ordinates: done"<<std::endl;
    this->findMaximaBetweenExtremes();
    //std::cout<<"MaximaBetweenExtremes: done"<<std::endl;
}



void Spline::removeAsymptotes() {

    // If the spline is completely flat, doesn't remove any segment
    if (numberOfAsymptotePolynomialsLeft == numberOfPolynomials)
        return;

    // If there are no horizontal asymptotes, doesn't remove any segment
    if (numberOfAsymptotePolynomialsLeft+numberOfAsymptotePolynomialsRight == 0)
        return;

    // If there are horizontal asymptotes, selects the new knots and the new
    // coefficients of the spline

    std::vector<double> newKnots;
    std::vector<std::vector<double>> newCoeffD0;
    std::vector<std::vector<double>> newCoeffD1;
    std::vector<std::vector<double>> newCoeffD2;

    for (int a=numberOfAsymptotePolynomialsLeft;
         a<numberOfPolynomials-numberOfAsymptotePolynomialsRight; ++a) {
        newKnots.push_back(knots[a]);
        newCoeffD0.push_back(coeffD0[a]);
        newCoeffD1.push_back(coeffD1[a]);
        newCoeffD2.push_back(coeffD2[a]);
    }
    newKnots.push_back(
        knots[numberOfPolynomials-numberOfAsymptotePolynomialsRight]);

    knots = newKnots;
    coeffD0 = newCoeffD0;
    coeffD1 = newCoeffD1;
    coeffD2 = newCoeffD2;

    this->updateVariables();

}



double Spline::D0(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<m; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D0(x)
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D1(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<g; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D1(x)
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D2(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<g-1; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D2(x)
    double y = 0;
    for (int i=0; i<g-1; ++i)
        y += coeffD2[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D0(const std::vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D0(powersOfX[1])
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D1(const std::vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D1(powersOfX[1])
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D0Shift(const std::vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots_shift[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D0(powersOfX[1])
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0_shift_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D1Shift(const std::vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots_shift[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D1(powersOfX[1])
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1_shift_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



void Spline::extendSpline(double lengthLeftSide, double lengthRightSide) {

    double D0Left = this->D0(knots[0]);
    double D0Right = this->D0(knots.back());

    double D1Left;
    double D1Right;
    calculateSlopeAtEndPoints(D1Left,D1Right);

    double knotLeft = knots[0];
    double knotRight = knots.back();

    auto segment = std::vector<double>(m,0);

    if (possibleNegativeOrdinates) {

        // Adds a diagonal line at the left of the spline
        if (lengthLeftSide > 0) {
            knots.insert(knots.begin(),knotLeft-lengthLeftSide);
            coeffD0.insert(coeffD0.begin(),segment);
            coeffD1.insert(coeffD1.begin(),segment);
            coeffD2.insert(coeffD2.begin(),segment);
            coeffD0[0][0] = D0Left-D1Left*knots[1];
            coeffD0[0][1] = D1Left;
            coeffD1[0][0] = D1Left;
        }

        // Adds a diagonal line at the right of the spline
        if (lengthRightSide > 0) {
            knots.push_back(knotRight+lengthRightSide);
            coeffD0.push_back(segment);
            coeffD1.push_back(segment);
            coeffD2.push_back(segment);
            coeffD0.back()[0] = D0Right-D1Right*knots[knots.size()-2];
            coeffD0.back()[1] = D1Right;
            coeffD1.back()[0] = D1Right;
        }

    }

    if (!possibleNegativeOrdinates) {

        // Adds a diagonal line at the left of the spline
        if (lengthLeftSide > 0)
            if (D0Left > 0) {
                knots.insert(knots.begin(),knotLeft-lengthLeftSide);
                coeffD0.insert(coeffD0.begin(),segment);
                coeffD1.insert(coeffD1.begin(),segment);
                coeffD2.insert(coeffD2.begin(), segment);
                coeffD0[0][0] = D0Left-D1Left*knots[1];
                coeffD0[0][1] = D1Left;
                coeffD1[0][0] = D1Left;
            }

        // Adds a diagonal line at the right of the spline
        if (lengthRightSide > 0)
            if (D0Right > 0) {
                knots.push_back(knotRight+lengthRightSide);
                coeffD0.push_back(segment);
                coeffD1.push_back(segment);
                coeffD2.push_back(segment);
                coeffD0.back()[0] = D0Right-D1Right*knots[knots.size()-2];
                coeffD0.back()[1] = D1Right;
                coeffD1.back()[0] = D1Right;
            }

        // Adds a horizontal line with ordinate 0 at the left of the spline
        if (lengthLeftSide > 0) {
            if (D0Left > 0) {
                double yFirstPoint = coeffD0[0][0]+coeffD0[0][1]*knots[0];
                if (yFirstPoint < 0) {
                    knots[0] = -coeffD0[0][0]/coeffD0[0][1];
                    knots.insert(knots.begin(),knotLeft-lengthLeftSide);
                    coeffD0.insert(coeffD0.begin(),segment);
                    coeffD1.insert(coeffD1.begin(),segment);
                    coeffD2.insert(coeffD2.begin(), segment);
                }
            }
            else if (D0Left == 0) {
                knots.insert(knots.begin(),knotLeft-lengthLeftSide);
                coeffD0.insert(coeffD0.begin(),segment);
                coeffD1.insert(coeffD1.begin(),segment);
                coeffD2.insert(coeffD2.begin(), segment);
            }
        }

        // Adds a horizontal line with ordinate 0 at the right of the spline
        if (lengthRightSide > 0) {
            if (D0Right > 0) {
                int i = coeffD0.size()-1;
                double yLastPoint = coeffD0[i][0]+coeffD0[i][1]*knots[i+1];
                if (yLastPoint < 0) {
                    knots[i+1] = -coeffD0[i][0]/coeffD0[i][1];
                    knots.push_back(knotRight+lengthRightSide);
                    coeffD0.push_back(segment);
                    coeffD1.push_back(segment);
                    coeffD2.push_back(segment);
                }
            }
            else if (D0Right == 0) {
                knots.push_back(knotRight+lengthRightSide);
                coeffD0.push_back(segment);
                coeffD1.push_back(segment);
                coeffD2.push_back(segment);
            }
        }

    }

    this->updateVariables();

}



void Spline::calculateShift(double Shift) {

    shift = Shift;

    auto powersShifts = std::vector<double>(m,1);

    for (int a=1; a<m; ++a)
        powersShifts[a] = powersShifts[a-1]*shift;

    coeffD0_shift_normalized =
        std::vector<std::vector<double>>(numberOfPolynomials,std::vector<double>(m,0));

    coeffD1_shift_normalized =
        std::vector<std::vector<double>>(numberOfPolynomials,std::vector<double>(m,0));

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<m; ++b)
            coeffD0_shift_normalized[a][b] = coeffD0_normalized[a][b];

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            coeffD1_shift_normalized[a][b] = coeffD1_normalized[a][b];

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            for (int c=1; c<=g-b; ++c) {
                if (c%2 != 0)
                    coeffD0_shift_normalized[a][b] -=
                    pascalsTriangle[b+c][b]*coeffD0_shift_normalized[a][b+c]*
                    powersShifts[c];
                else
                    coeffD0_shift_normalized[a][b] +=
                    pascalsTriangle[b+c][b]*coeffD0_shift_normalized[a][b+c]*
                    powersShifts[c];
            }

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            for (int c=1; c<=g-b; ++c) {
                if (c%2 != 0)
                    coeffD1_shift_normalized[a][b] -=
                    pascalsTriangle[b+c][b]*coeffD1_shift_normalized[a][b+c]*
                    powersShifts[c];
                else
                    coeffD1_shift_normalized[a][b] +=
                    pascalsTriangle[b+c][b]*coeffD1_shift_normalized[a][b+c]*
                    powersShifts[c];
            }

    knots_shift = std::vector<double>(numberOfKnots);

    for (int a=0; a<numberOfKnots; ++a)
        knots_shift[a] = knots[a]+shift;

}



void Spline::normalizeCoefficients(double verticalShift,
                                   double normalizationFactorD0,
                                   double normalizationFactorD1) {

    coeffD0_normalized = coeffD0;
    coeffD1_normalized = coeffD1;

    if (verticalShift != 0)
        for (int a=0; a<numberOfPolynomials; ++a)
            coeffD0_normalized[a][0] += verticalShift;

    if (normalizationFactorD0 != 1)
        for (int a=0; a<numberOfPolynomials; ++a)
            for (int b=0; b<m; ++b)
                coeffD0_normalized[a][b] /= normalizationFactorD0;

    if (normalizationFactorD1 != 1)
        for (int a=0; a<numberOfPolynomials; ++a)
            for (int b=0; b<g; ++b)
                coeffD1_normalized[a][b] /= normalizationFactorD1;

}



void Spline::chooseKnots(int numberOfAbscissaeSeparatingConsecutiveKnots) {

    int number = numberOfAbscissaeSeparatingConsecutiveKnots;

    double meanKnotDistance =
        (abscissae.back()-abscissae[0]) / (double)(abscissae.size()-1);
    //std::cout<< "Mean Knot Distance: "<< meanKnotDistance << std::endl;
    //only for the model
    if (splineType == 1 /*Model*/) {

        std::vector<double> newX;
        std::vector<double> newY;

        newX.push_back(abscissae[0]);
        newY.push_back(ordinates[0]);

        // If there are less than 30 points, adds enough points to the spline to
        // reach at least 30 points
        if (abscissae.size() < 30) {
            
            // span of the the domaine
            double abscissaeLength = (abscissae.back()-abscissae[0]);
            // he wants to have 30 points for the models
            int minPointsToAdd = 30-abscissae.size();
            // start from the position indexed 1, which is the second actually
            for (int a=1; a<abscissae.size(); ++a) 
            {
                // take the distance between two points
                double segmentLength = (abscissae[a]-abscissae[a-1]);
                // number of points to add between a-th and a-1-th position is a weighted number over the missing ones to add
                int numberOfPointstoAdd =
                    segmentLength/abscissaeLength*(double)(minPointsToAdd+1);
                // the points will be equidistant from eachother between a-th and a-1-th position
                double distanceBetweenPoints =
                    segmentLength/(double)(numberOfPointstoAdd+1);
                // it computes the slope of the straight line between a-th and a-1-th position
                double slope = (ordinates[a]-ordinates[a-1])/segmentLength;
                // loops over until it adds as many points as needed according to previous calculations
                for (int b=0; b<numberOfPointstoAdd; ++b) {
                    //pushes back to newX the new abscissa as the sum the back value in newX and distanceBetweenPoints 
                    newX.push_back(newX.back()+distanceBetweenPoints);
                    //pushes back to newY the new ordinate as the sum of ordinates[a-1] and the slope, multiplied with the
                    //difference between brand new newX we just added and abscissae[a-1]
                    newY.push_back(
                        ordinates[a-1]+slope*(newX.back()-abscissae[a-1]));
                }
            // i have done all the filling between previously empy spaces, we also add the values corresponding to a
            newX.push_back(abscissae[a]);
            newY.push_back(ordinates[a]);
            }
        }

        // Adds extra points between consecutive data points with a distance on
        // the x-axis greater than 3.*meanKnotDistance

        // this is for the case where we have many points (more or equal than 30), but some of them are to far from each other.
        if (abscissae.size() >= 30)
            for (int a=1; a<abscissae.size(); ++a) 
            {
                double segmentLength = (abscissae[a]-abscissae[a-1]);
                // if this segmentLength is huge, than it will populate the empty range in between the two points (same stuff as before)
                if (segmentLength > 3.*meanKnotDistance) 
                {
                    int numberOfNewPoints =
                        (int)(segmentLength/meanKnotDistance);
                    double distanceBetweenPoints =
                        segmentLength / (double)(numberOfNewPoints+1);
                    double slope = (ordinates[a]-ordinates[a-1])/segmentLength;
                    for (int b=0; b<numberOfNewPoints; ++b) 
                    {
                        newX.push_back(newX.back()+distanceBetweenPoints);
                        newY.push_back(
                            ordinates[a-1]+slope*(newX.back()-abscissae[a-1]));
                    }
                }
                newX.push_back(abscissae[a]);
                newY.push_back(ordinates[a]);
            }
        // assigns new values
        abscissae = newX;
        ordinates = newY;
        // recompute n
        n = abscissae.size();
        // recompute mean knot distance
		meanKnotDistance =
			(abscissae.back()-abscissae[0]) / (double)(abscissae.size()-1);

    }

    // finds max and min ordinate
	double maxOrdinate = ordinates[0];
	double minOrdinate = ordinates[0];
	for (int a=1; a<n; ++a) {
		if (ordinates[a] > maxOrdinate)
			maxOrdinate = ordinates[a];
		if (ordinates[a] < minOrdinate)
			minOrdinate = ordinates[a];
	}
//	std::cout<< "Max Ordinate: " <<  maxOrdinate << std::endl;
//    	std::cout<< "Min Ordinate: " <<  minOrdinate << std::endl;
     // compute the height
	double height = maxOrdinate - minOrdinate;
//    	std::cout<< "Height: " <<  height << std::endl;
    // knots is a vector of doubles, containing the abscissas of the knots
    knots.push_back(abscissae[0]);

    if (abscissae.size() > 2) {
        
        double y = ordinates[0];
        int k = 0;
		int l = 0;

        // loop over from the second abscissa till the penultimate position
        for (int a=1; a<abscissae.size()-1; ++a) {
            ++k;
            // compute the difference
            double difference = abscissae[a] - abscissae[a-1];

            // The following routine clearly serves the purpouse of erasing the outliers
            
            // CONDITIONS
            // 1 - k will always be > number, because number is set to 0 by default when calling the splines for both models and experimental data;
            // 2 - If the values are far from each other in terms of abscissa, by more than 2*meanKnotDistance;
            // 3 - If the distance in height is significant;
            // 4 - If the previous added point was significantly higher than the earlier one.
            if (k > number ||
				difference > 2.*meanKnotDistance ||
				fabs(ordinates[a]-ordinates[a-1]) > 0.1*height ||
                l > 0)
                // this will be more or less always satisfied for models, and experiments
                // but this is cool if you have outliers, because we will arrive here thanks to condition 3
                // and then finally we'll exclude the point from Knots because it's too close to the previous
                if (difference > 0.2*meanKnotDistance)
                    // only if ordinates[a] is different from ordinates[0]
                    // this way we get rid of the points with null derivative
                    // in this case k will not be set to 0 and l will remain equal to previous value
                    if (ordinates[a] != y) 
                    {
                        // in this case we push back the ordinate values knots
                        knots.push_back(abscissae[a]);
                        // we update y with the new value
                        y = ordinates[a];
                        // k is again set to 0
                        k = 0;
                        // l becomes negative at first
                        --l;
                        // then only if the difference in ordinate is significant it will be greater than 0
						if (fabs(ordinates[a]-ordinates[a-1]) > 0.1*height)
							l = number+1;
                    }
        }

        // Improves the positioning of the knots if the data ends with a
        // horizontal asymptote. Some of the intervals chosen with this approach
        // might not contain any data points
        if (y == ordinates.back() && k > 10) { 
            double knot = knots.back();
            knots.push_back(knot+(abscissae.back()-knot)*4./12.);
            knots.push_back(knot+(abscissae.back()-knot)*8./12.);
            knots.push_back(knot+(abscissae.back()-knot)*10./12.);
            knots.push_back(knot+(abscissae.back()-knot)*11./12.);
        }

    }

    // of course the last abscissa will always be part of the knots vector
    knots.push_back(abscissae.back());
    //for (int i=0; i<knots.size(); ++i)
    //{
    //    std::cout<<"KNOT  "<< i <<" is equal to"<<knots[i]<<std::endl;
    //}
    // Sets the values of numberOfKnots, numberOfPolynomials, K and G
    numberOfKnots = knots.size();
    numberOfPolynomials = numberOfKnots - 1;
    // m is the order of the basis function
    // K is the number of degrees of freedom of the spline
    K = numberOfKnots - 2 + m;
    // 
    G = K-1;

    // xRange is the domain, which is the same of the original abscissa, because they are added in the beginning and in the end
    // of the knot selection process
    xRange = knots.back() - knots[0];

    // Fills knotsForCalculations with the current knots plus additional knots
    // on the left and on the right of the spline, each at a distance from the
    // nearest knot equal to the mean distance of the other knots

    // OK
    double meanDistance = xRange / (double)numberOfPolynomials;

    // the number of non used knots, which will be added on the two sides of the spline is equal to 
    // to the degrees of the polynomials
    knotsForCalculations = std::vector<double>(numberOfKnots+2*g,0);
    for (int i=0; i<g; ++i)
        knotsForCalculations[i] = knots[0] + (double)(i-g)*meanDistance;
    for (int i=0; i<numberOfKnots; ++i)
        knotsForCalculations[i+g] = knots[i];
    for (int i=1; i<m; ++i)
        knotsForCalculations[i+g+numberOfPolynomials] =
            knots.back()+(double)i*meanDistance;
    
    //for (int i=0; i<knotsForCalculations.size(); ++i)
    //{
    //	std::cout<<"KNOT FOR CALCULATION "<< i <<" is equal to"<<knotsForCalculations[i]<<std::endl;	
    //}
}



void Spline::calculateCoefficients() {

    // Calculates the basis functions
    //std::cout<<"enters the calculate coeffiecient function"<<std::endl;
    auto basisFunctions = std::vector<BasisFunction>(K);
    for (int j=0; j<K; ++j)
        basisFunctions[j].calculateCoefficients(j,knotsForCalculations);
    //std::cout<<"Basis function: Initialised"<<std::endl;
    // Calculates the Fi matrix
    auto Fi = std::vector<std::vector<double>>(n,std::vector<double>(K,0));
    for (int i=0; i<n; ++i)
        for (int j=0; j<K; ++j)
	{
            	Fi[i][j] = basisFunctions[j].D0(abscissae[i]);
		//std::cout<<"y values from the Fi"<<i<<" "<<j<<" matrix is: "<<Fi[i][j]<<std::endl;
    	}
	// Finds the limits for the non-zero elements in Fi
    auto firstInFi = std::vector<int>(n,0);
    auto lastInFi = std::vector<int>(n,0);
    for (int i=0; i<n; ++i)
        for (int j=0; j<K; ++j)
            if (Fi[i][j] != 0) {
                firstInFi[i] = j;
                break;
            }
    for (int i=0; i<n; ++i)
        for (int j=G; j>-1; --j)
            if (Fi[i][j] != 0) {
                lastInFi[i] = j;
                break;
            }
    //std::cout<<"Here: 2"<<std::endl;
    // Finds the limits for the non-zero elements in FiT
    auto firstInFiT = std::vector<int>(K,0);
    auto lastInFiT = std::vector<int>(K,0);
    for (int j=0; j<K; ++j)
        for (int i=0; i<n; ++i)
            if (Fi[i][j] != 0) {
                firstInFiT[j] = i;
                break;
            }
    for (int j=0; j<K; ++j)
        for (int i=n-1; i>-1; --i)
            if (Fi[i][j] != 0) {
                lastInFiT[j] = i;
                break;
            }
    //std::cout<<"Here: 3"<<std::endl;
    // Finds the limits for the non-zero elements in M, R and FiTFi
    auto firstInBandMatrices = std::vector<int>(K,0);
    auto lastInBandMatrices = std::vector<int>(K,G);
    if (K > m)
        for (int i=m; i<K; ++i)
            firstInBandMatrices[i] = i-g;
    if (K > m)
        for (int i=0; i<K-m; ++i)
            lastInBandMatrices[i] = i+g;
    //std::cout<<"Here: 4"<<std::endl;
    // Calculates the FiTFi matrix, equal to the product of FiT and Fi
    auto FiTFi = std::vector<std::vector<double>>(K,std::vector<double>(K,0));
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            for (int k=0; k<n; ++k)
		{
                FiTFi[i][j] += Fi[k][i] * Fi[k][j];
    		//std::cout<<FiTFi[i][j]<<std::endl;
		}
    //std::cout<<"Here: 5"<<std::endl;
    // Calculates the R matrix
    auto R = std::vector<std::vector<double>>(K,std::vector<double>(K,0));
    for (int i=0; i<K; ++i)
	{
        for (int j=i; j<=lastInBandMatrices[i]; ++j)
		{
            	R[i][j] = basisFunctions[i].integralOfProductD2(basisFunctions[j]);
		//std::cout<<"integral [i][j]"<< R[i][j] <<std::endl;
    		}
	}    	
    for (int i=1; i<K; ++i)
    {
	for (int j=firstInBandMatrices[i]; j<i; ++j)
	{
            R[i][j] = R[j][i];
	    //std::cout<<"integral ["<<i<<"]["<<j<<"]"<< R[i][j] <<std::endl;
	}
    }

    // Calculates the FiTy vector
    auto FiTy = std::vector<double>(K,0);
    for (int i=0; i<K; ++i)
        for (int k=firstInFiT[i]; k<=lastInFiT[i]; ++k)
            FiTy[i] += Fi[k][i] * ordinates[k];
    // Estimates the first derivatives of the experimental data. Contains an
    // additional 0 at position 0
    auto estimatedD1 = std::vector<double>(n-1,0);
    for (int i=1; i<n-1; ++i)
        estimatedD1[i] =
            (ordinates[i+1]-ordinates[i-1]) / (abscissae[i+1]-abscissae[i-1]);
    // Calculates the square root of the sum of squares of the elements of FiTFi
    double indexFiTFi = 0;
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            indexFiTFi += FiTFi[i][j] * FiTFi[i][j];
    indexFiTFi = sqrt(indexFiTFi);
    // Calculates the square root of the sum of squares of the elements of R
    double indexR = 0;
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            indexR += R[i][j] * R[i][j];
    indexR = sqrt(indexR);
    // Calculates the value of log10lambda for which the elements of FiTFi and R
    // have the same order of magnitude, rounded to the nearest 0.5
    log10lambdaForSameOrderOfMagnitude =
        round(2.*(log10(indexFiTFi)-log10(indexR)))/2.;
    // Calculates the end points of the log10lambda minimization interval
    log10lambdaMin = log10lambdaForSameOrderOfMagnitude-lambdaSearchInterval/2.;
    log10lambdaMax = log10lambdaForSameOrderOfMagnitude+lambdaSearchInterval/2.;
    // Calculates the log10 of the distance between two consecutive steps in the
    // for cycle for minimizing log10lambda
    double log10lambdaStep =
        lambdaSearchInterval/(double)(numberOfStepsLambda-1);

    // Initializes the elements necessary for the minimization
    auto M = std::vector<std::vector<double>>(K,std::vector<double>(K,0));
    auto zed = std::vector<double>(K,0);
    auto Z = std::vector<std::vector<double>>(K,std::vector<double>(K,0));
    auto Minv = std::vector<std::vector<double>>(K,std::vector<double>(K,0));
    auto MinvFiT = std::vector<std::vector<double>>(K,std::vector<double>(n,0));
    auto splineCoefficientsForVariousLambdas =
        std::vector<std::vector<double>>(numberOfStepsLambda,std::vector<double>(K,0));
    auto GCV1 = std::vector<double>(numberOfStepsLambda,0);

    // Calculates the spline coefficients and GCV1 for each lambda in the for
    // cycle
    for (int a=0; a<numberOfStepsLambda; ++a) {

        // Obtains the value of lambda from that of a
        lambda = pow(10., log10lambdaMin + (double)a * log10lambdaStep);

        // Calculates the M matrix, sum of FiTFi and the product of Lambda and R
        for (int i=0; i<K; ++i)
            for (int j=firstInBandMatrices[i];j<=lastInBandMatrices[i];++j)
                M[i][j] = FiTFi[i][j] + lambda * R[i][j];

        // Uses the Doolittle decomposition to decompose M, and obtains the
        // triangular matrices L and U (M = L*U). Saves the values of the L and
        // U matrices, except for the main diagonal of L (consisting entirely of
        // ones), in place of the respective values in matrix M
        for (int i=0; i<K; ++i) {
            for (int j=i; j<=lastInBandMatrices[i]; ++j)
                for (int k=0; k<i; ++k)
                    M[i][j] -= M[i][k] * M[k][j];
            for (int j=i+1; j<=lastInBandMatrices[i]; ++j) {
                for (int k=0; k<i; ++k)
                    M[j][i] -= M[j][k] * M[k][i];
                M[j][i] /= M[i][i];
            }
        }

        // Calculates the zed vector, from the expression L*zed = FiTy, using
        // the forward substitution technique
        zed[0] = FiTy[0];
        for (int i=1; i<K; ++i) {
            zed[i] = FiTy[i];
            for (int k=firstInBandMatrices[i]; k<i; ++k)
                zed[i] -= M[i][k] * zed[k];
        }

        // Calculates the spline coefficients from R*coefficients = zed
        splineCoefficientsForVariousLambdas[a][G] = zed[G] / M[G][G];
        for (int i=G-1; i>-1; --i) {
            splineCoefficientsForVariousLambdas[a][i] = zed[i];
            for (int k=lastInBandMatrices[i]; k>i; --k)
                splineCoefficientsForVariousLambdas[a][i] -=
                M[i][k] * splineCoefficientsForVariousLambdas[a][k];
            splineCoefficientsForVariousLambdas[a][i] /= M[i][i];
        }

        // Solves L*Z = I
        for (int j=0; j<K; ++j) {
            Z[j][j] = 1.;
            for (int i=j+1; i<K; ++i) {
                Z[i][j] = 0;
                for (int k=firstInBandMatrices[i]; k<i; ++k)
                    Z[i][j] -= M[i][k] * Z[k][j];
            }
        }

        // Solves R*Minv = Z
        for (int j=G; j>-1; --j) {
            Minv[G][j] = Z[G][j] / M[G][G];
            for (int i=G-1; i>-1; --i) {
                if (j < i)
                    Minv[i][j] = Z[i][j];
                else if (j > i)
                    Minv[i][j] = 0;
                else
                    Minv[i][j] = 1.;
                for (int k=lastInBandMatrices[i]; k>i; --k)
                    Minv[i][j] -= M[i][k] * Minv[k][j];
                Minv[i][j] /= M[i][i];
            }
        }

        // Calculates the MinvFiT matrix, equal to the product of Minv and FiT,
        // in the locations necessary for the calculation of the trace of S
        for (int j=0; j<n; ++j)
            for (int i=firstInFi[j]; i<=lastInFi[j]; ++i) {
                MinvFiT[i][j] = 0;
                for (int k=firstInFi[j]; k<=lastInFi[j]; ++k)
                    MinvFiT[i][j] += Minv[i][k] * Fi[j][k];
            }

        // Calculates the numerator of GCV1(Lambda)
        double SSE1 = 0; // Sum of squared errors between yi' and f'(xi)
        for (int i=1; i<n-1; ++i) {
            double difference = estimatedD1[i];
            for (int j=0; j<K; ++j)
                difference -=
                splineCoefficientsForVariousLambdas[a][j] *
                basisFunctions[j].D1(abscissae[i]);
            SSE1 += difference * difference;
        }
        GCV1[a] = (double)n * SSE1;

        // Calculates the trace of matrix S
        double traceS = 0;
        for (int i=0; i<n; ++i)
            for (int k=firstInFi[i]; k<=lastInFi[i]; ++k)
                traceS += Fi[i][k] * MinvFiT[k][i];

        // Calculates GCV1(Lambda)
        GCV1[a] /= ((double)n-traceS)*((double)n-traceS);

    } // End of the for cycle for each lambda
    // Finds the minimum value of GCV1(lambda) in the GCV1 vector, and saves the
    // corresponding spline coefficients to splineCoefficients and the
    // corresponding lambda and log10(lambda) to 'lambda' and 'log10lambda'
    double GCVOne = GCV1[0];
    int index = 0;
    log10lambda = log10lambdaMin;
    lambda = pow(10.,log10lambda);
    for (int i=1; i<numberOfStepsLambda; ++i)
        if (GCV1[i] < GCVOne) {
            GCVOne = GCV1[i];
            index = i;
            log10lambda = log10lambdaMin + (double)i * log10lambdaStep;
            lambda = pow(10.,log10lambda);
        }
    splineCoefficients = splineCoefficientsForVariousLambdas[index];

    // Calculates the coefficients of the polynomials of the spline
    coeffD0 = std::vector<std::vector<double>>(numberOfPolynomials,std::vector<double>(m,0));
    int firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<m; ++c) {
                coeffD0[a-g][c] += basisFunctions[b].coeffD0[g+firstBasis-b][c]*
                                   splineCoefficients[b];
		//std::cout<<"coeffD0["<<a-g<<"]["<<c<<"] is equal to "<<coeffD0[a-g][c]<<std::endl;
            }
        }
        ++firstBasis;
    }

    // Calculates the coefficients of the first derivative of the polynomials of
    // the spline
    coeffD1 = std::vector<std::vector<double>>(numberOfPolynomials,std::vector<double>(m,0));
    firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<g; ++c) {
                coeffD1[a-g][c] += basisFunctions[b].coeffD1[g+firstBasis-b][c]*
                                   splineCoefficients[b];
            	//std::cout<<"coeffD1["<<a-g<<"]["<<c<<"] is equal to "<<coeffD1[a-g][c]<<std::endl;
		}
        }
        ++firstBasis;
    }

    // Calculates the coefficients of the second derivative of the polynomials
    // of the spline
    coeffD2 = std::vector<std::vector<double>>(numberOfPolynomials, std::vector<double>(m,0));
    firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<g-1; ++c) {
                coeffD2[a-g][c] += basisFunctions[b].coeffD2[g+firstBasis-b][c]*
                                   splineCoefficients[b];
            	//std::cout<<"coeffD2["<<a-g<<"]["<<c<<"] is equal to "<<coeffD2[a-g][c]<<std::endl;
		}
        }
        ++firstBasis;
    }

    // Initializes the 'powers' vector
    powers = std::vector<double>(m,1);
}



void Spline::calculateSlopeAtEndPoints(double& slopeLeftSide,
    double& slopeRightSide) {

    slopeLeftSide =
        (this->D0(knots[1])-this->D0(knots[0])) / (knots[1]-knots[0]);

    slopeRightSide =
        (this->D0(knots.back())-this->D0(knots[numberOfKnots-2])) /
        (knots.back()-knots[numberOfKnots-2]);

}



void Spline::yAndAsymptoteAnalysis() {

    std::vector<double> rootsD0 = this->calculateRoots(1);

    std::vector<double> pointsD0;

    if (rootsD0.size() > 0)
        for (int a=0; a<rootsD0.size(); ++a)
            pointsD0.push_back(this->D0(rootsD0[a]));
    for (int a=0; a<numberOfKnots; ++a)
        pointsD0.push_back(this->D0(knots[a]));

    std::vector<double> pointsD0Original;

    if (rootsD0.size() > 0)
        for (int a=0; a<rootsD0.size(); ++a)
            if (rootsD0[a] >= originalAbscissae[0])
                if (rootsD0[a] <= originalAbscissae.back())
                    pointsD0Original.push_back(this->D0(rootsD0[a]));
    
    
    for (int a=0; a<numberOfKnots; ++a)
    {
        	if (knots[a] >= originalAbscissae[0])
            	{
			if (knots[a] <= originalAbscissae.back())
            		{
				pointsD0Original.push_back(this->D0(knots[a]));
				if(splineType == 0)
    				{
    					if(print_splines == true)
    					{
						std::cout<< "Abscissa	"<< knots[a]<< "   Ordinate	"<<pointsD0Original.back()<<std::endl;
	    				}
				}
			}
		}
     }


    // Calculates yD0Max
    yD0Max = pointsD0[0];
    for (int a=1; a<pointsD0.size(); ++a)
        if (pointsD0[a] > yD0Max)
            yD0Max = pointsD0[a];

    // Calculates yD0Min
    yD0Min = pointsD0[0];
    for (int a=1; a<pointsD0.size(); ++a)
        if (pointsD0[a] < yD0Min)
            yD0Min = pointsD0[a];

    // Calculates yD0MaxOriginal
    yD0MaxOriginal = pointsD0Original[0];
    for (int a=1; a<pointsD0Original.size(); ++a)
        if (pointsD0Original[a] > yD0MaxOriginal)
            yD0MaxOriginal = pointsD0Original[a];

    // Calculates yD0MinOriginal
    yD0MinOriginal = pointsD0Original[0];
    for (int a=1; a<pointsD0Original.size(); ++a)
        if (pointsD0Original[a] < yD0MinOriginal)
            yD0MinOriginal = pointsD0Original[a];

    // Calculates yD0Range
    yD0Range = yD0Max - yD0Min;

    std::vector<double> rootsD1 = this->calculateRoots(2);

    std::vector<double> pointsD1;

    if (rootsD1.size() > 0)
        for (int a=0; a<rootsD1.size(); ++a)
            pointsD1.push_back(this->D1(rootsD1[a]));
    for (int a=0; a<numberOfKnots; ++a)
	{
        	pointsD1.push_back(this->D1(knots[a]));
		//std::cout<< "Abscissa'   "<< knots[a]<< "   Ordinate'     "<<pointsD1.back()<<std::endl;
	}

    std::vector<double> pointsD1Original;

    if (rootsD1.size() > 0)
        for (int a=0; a<rootsD1.size(); ++a)
            if (rootsD1[a] >= originalAbscissae[0])
                if (rootsD1[a] <= originalAbscissae.back())
                    pointsD1Original.push_back(this->D1(rootsD1[a]));
    for (int a=0; a<numberOfKnots; ++a)
        if (knots[a] >= originalAbscissae[0])
            if (knots[a] <= originalAbscissae.back())
                pointsD1Original.push_back(this->D1(knots[a]));

    // Calculates yD1Max
    yD1Max = pointsD1[0];
    for (int a=1; a<pointsD1.size(); ++a)
        if (pointsD1[a] > yD1Max)
            yD1Max = pointsD1[a];

    // Calculates yD1Min
    yD1Min = pointsD1[0];
    for (int a=1; a<pointsD1.size(); ++a)
        if (pointsD1[a] < yD1Min)
            yD1Min = pointsD1[a];

    // Calculates yD1MaxOriginal
    yD1MaxOriginal = pointsD1Original[0];
    for (int a=1; a<pointsD1Original.size(); ++a)
        if (pointsD1Original[a] > yD1MaxOriginal)
            yD1MaxOriginal = pointsD1Original[a];

    // Calculates yD1MinOriginal
    yD1MinOriginal = pointsD1Original[0];
    for (int a=1; a<pointsD1Original.size(); ++a)
        if (pointsD1Original[a] < yD1MinOriginal)
            yD1MinOriginal = pointsD1Original[a];

    // Calculates yD1MaxAbs
    yD1MaxAbs = fabs(yD1Max);
    if (fabs(yD1Min) > fabs(yD1Max))
        yD1MaxAbs = fabs(yD1Min);

    double minVariation =
        yD0Range*fractionOfOrdinateRangeForAsymptoteIdentification;

    // Finds the number of polynomials which count as asymptotes on the left of
    // the spline

    numberOfAsymptotePolynomialsLeft = 0;
    double yFront = this->D0(knots[0]);

    for (int a=1; a<numberOfKnots; ++a) {
        double y = this->D0(knots[a]);
        if (fabs(yFront-y) < minVariation)
            ++numberOfAsymptotePolynomialsLeft;
        else
            break;
    }

    // Finds the number of polynomials which count as asymptotes on the right of
    // the spline

    numberOfAsymptotePolynomialsRight = 0;
    double yBack = this->D0(knots.back());

    for (int a=numberOfKnots-2; a>-1; --a) {
        double y = this->D0(knots[a]);
        if (fabs(yBack-y) < minVariation)
            ++numberOfAsymptotePolynomialsRight;
        else
            break;
    }

    // Calculates xMaxNoAsymptotes and xMinNoAsymptotes

    xMaxNoAsymptotes =
        knots[numberOfPolynomials-numberOfAsymptotePolynomialsRight];

    xMinNoAsymptotes = knots[numberOfAsymptotePolynomialsLeft];
    if (xMinNoAsymptotes < 0)
        xMinNoAsymptotes = 0;

    // Calculates xRange
    xRange = knots.back() - knots[0];

}



std::vector<double> Spline::calculateRoots(double derivativeOrder) {

    std::vector<double> roots;

    for (int i=0; i<numberOfPolynomials; ++i) {

        // y = a0*x^3 + b0*x^2 + c0*x + d0
        double a0, b0, c0, d0;

        if (derivativeOrder == 0) {
            a0 = coeffD0[i][3];
            b0 = coeffD0[i][2];
            c0 = coeffD0[i][1];
            d0 = coeffD0[i][0];
        }

        if (derivativeOrder == 1) {
            a0 = coeffD1[i][3];
            b0 = coeffD1[i][2];
            c0 = coeffD1[i][1];
            d0 = coeffD1[i][0];
        }

        // Case: polynomial degree = 0. If the ordinate is 0, returns the knots
        // at the end points
        if (a0 == 0 && b0 == 0 && c0 == 0 && d0 == 0) {

            if (roots.size() == 0) {
                roots.push_back(knots[i]);
                if (knots[i+1] != knots[i])
                    roots.push_back(knots[i+1]);
            }
            else
                if (knots[i] != roots.back()) {
                    roots.push_back(knots[i]);
                    if (knots[i+1] != knots[i])
                        roots.push_back(knots[i+1]);
                }

        }

        // Case: polynomial degree = 1
        if (a0 == 0 && b0 == 0 && c0 != 0) {

            double zero = -d0/c0;

            if (zero >= knots[i] && zero <= knots[i+1])
                if (roots.size() == 0)
                    roots.push_back(zero);
                else
                    if (zero != roots.back())
                        roots.push_back(zero);

        }

        // Case: polynomial degree = 2
        if (a0 == 0 && b0 != 0) {

            // y = a*x^2 + b*x + c
            double a = b0;
            double b = c0;
            double c = d0;

            double delta = b*b-4*a*c;

            if (delta == 0) {
                double zero = -b/(2.*a);
                if (zero >= knots[i] && zero <= knots[i+1])
                    if (roots.size() == 0)
                        roots.push_back(zero);
                    else
                        if (zero != roots.back())
                            roots.push_back(zero);
            }

            if (delta > 0) {
                double x1 = (-b-sqrt(delta))/(2.*a);
                double x2 = (-b+sqrt(delta))/(2.*a);
                if (x1 >= knots[i] && x1 <= knots[i+1])
                    if (roots.size() == 0)
                        roots.push_back(x1);
                    else
                        if (x1 != roots.back())
                            roots.push_back(x1);
                if (x2 >= knots[i] && x2 <= knots[i+1])
                    if (roots.size() == 0)
                        roots.push_back(x2);
                    else
                        if (x2 != roots.back())
                            roots.push_back(x2);
            }

        }

        // Case: polynomial degree = 3
        if (a0 != 0) {

            // y = x^3 + a*x^2 + b*x + c
            double a = b0/a0;
            double b = c0/a0;
            double c = d0/a0;

            double Q = (a*a-3.*b)/9.;
            double R = (2.*a*a*a-9.*a*b+27.*c)/54.;

            double S = Q*Q*Q - R*R;

            double fi;
            double pi = atan(1.)*4.;

            if (S > 0) {

                fi = 1./3.*acos(R/sqrt(Q*Q*Q));

                double x1 = -2.*sqrt(Q)*cos(fi)-a/3.;
                double x2 = -2.*sqrt(Q)*cos(fi+2./3.*pi)-a/3.;
                double x3 = -2.*sqrt(Q)*cos(fi-2./3.*pi)-a/3.;

                std::vector<double> solutions;
                if (x1 >= knots[i] && x1 <= knots[i+1]) solutions.push_back(x1);
                if (x2 >= knots[i] && x2 <= knots[i+1]) solutions.push_back(x2);
                if (x3 >= knots[i] && x3 <= knots[i+1]) solutions.push_back(x3);
                if (solutions.size() > 0) {
                    sort(solutions.begin(),solutions.end());
                    for (int j=0; j<solutions.size(); ++j)
                        if (roots.size() == 0)
                            roots.push_back(solutions[j]);
                        else
                            if (solutions[j] != roots.back())
                                roots.push_back(solutions[j]);
                }

            }

            if (S < 0) {

                double signOfR;
                if (R == fabs(R)) signOfR = 1.;
                else signOfR = -1.;

                if (Q > 0) {
                    fi = 1./3.*acosh(fabs(R)/sqrt(Q*Q*Q));
                    double zero = -2.*signOfR*sqrt(Q)*cosh(fi)-a/3.;
                    if (zero >= knots[i] && zero <= knots[i+1])
                        if (roots.size() == 0)
                            roots.push_back(zero);
                        else
                            if (zero != roots.back())
                                roots.push_back(zero);
                }

                if (Q < 0) {
                    fi = 1./3.*asinh(fabs(R)/sqrt(fabs(Q*Q*Q)));
                    double zero = -2.*signOfR*sqrt(fabs(Q))*sinh(fi)-a/3.;
                    if (zero >= knots[i] && zero <= knots[i+1])
                        if (roots.size() == 0)
                            roots.push_back(zero);
                        else
                            if (zero != roots.back())
                                roots.push_back(zero);
                }

                if (Q == 0) {
                    double zero = -pow(c-a*a*a/27.,1./3.)-a/3.;
                    if (zero >= knots[i] && zero <= knots[i+1])
                        if (roots.size() == 0)
                            roots.push_back(zero);
                        else
                            if (zero != roots.back())
                                roots.push_back(zero);
                }

            }

            if (S == 0) {

                double signOfR;
                if (R == fabs(R)) signOfR = 1.;
                else signOfR = -1.;

                if (Q >= 0) {
                    double x1 = -2.*signOfR*sqrt(Q)-a/3.;
                    double x2 = signOfR*sqrt(Q)-a/3.;
                    std::vector<double> solutions;
                    if (x1 >= knots[i] && x1 <= knots[i+1])
                        solutions.push_back(x1);
                    if (x2 >= knots[i] && x2 <= knots[i+1])
                        solutions.push_back(x2);
                    if (solutions.size() > 0) {
                        sort(solutions.begin(),solutions.end());
                        for (int j=0; j<solutions.size(); ++j)
                            if (roots.size() == 0)
                                roots.push_back(solutions[j]);
                            else
                                if (solutions[j] != roots.back())
                                    roots.push_back(solutions[j]);
                    }
                }

            }

        } // End of Case: polynomial degree = 3

    } // End of the for cycle for each polynomial

    return roots;

}



void Spline::removeNegativeSegments() {

    // Calculates the real different roots of the spline
    std::vector<double> roots = calculateRoots(0);

    if (roots.size() == 0) return;

    // If required, adds new knots, corresponding to the roots, to the spline,
    // and adds the corresponding polynomials to the polynomials of the spline
    int i=knots.size()-1;
    int j=0;
    for (int a=0; a<i; ++a)
        for (int b=j; b<roots.size(); ++b)
            if (roots[b] > knots[a] && roots[b] < knots[a+1]) {
                ++a;
                ++i;
                ++j;
                knots.insert(knots.begin()+a,roots[b]);
                coeffD0.insert(coeffD0.begin()+a,coeffD0[a-1]);
                coeffD1.insert(coeffD1.begin()+a,coeffD1[a-1]);
                coeffD2.insert(coeffD2.begin()+a,coeffD2[a-1]);
            }

    // Updates the values of numberOfKnots, numberOfPolynomials, K and G
    numberOfKnots = knots.size();

    // Replaces any negative segments of the spline with straight lines with
    // ordinate 0
    auto segment = std::vector<double>(m,0);
    for (int a=0; a<knots.size()-1; ++a) {
        double midpoint = (knots[a]+knots[a+1])/2.;
        if (this->D0(midpoint) <= 0) {
            coeffD0[a] = segment;
            coeffD1[a] = segment;
            coeffD2[a] = segment;
        }
    }

    // Replaces any two consecutive segments with the same polynomial with a
    // single segment with that polynomial
    int k=knots.size();
    for (int a=1; a<k; ++a)
        if (coeffD0[a] == coeffD0[a-1]) {
            knots.erase(knots.begin()+a);
            coeffD0.erase(coeffD0.begin()+a);
            coeffD1.erase(coeffD1.begin()+a);
            coeffD2.erase(coeffD2.begin()+a);
            --a;
            --k;
        }

    this->updateVariables();

}



void Spline::updateVariables() {

    // Updates the values of numberOfKnots, numberOfPolynomials, K and G
    numberOfKnots = knots.size();
    numberOfPolynomials = numberOfKnots-1;
    K = numberOfKnots-2+m;
    G = K-1;

    this->yAndAsymptoteAnalysis();

    this->findMaximaBetweenExtremes();

}



void Spline::findMaximaBetweenExtremes() {

    std::vector<double> allRootsD1 = calculateRoots(1);

    if (allRootsD1.size() == 0) {
        hasOneMaximum = false;
        hasTwoMaxima = false;
        hasFirstOfTwoMaxima = false;
        return;
     }

    std::vector<double> rootsD1;

    for (int a=0; a<allRootsD1.size(); ++a)
        if (allRootsD1[a] != knots[0] && allRootsD1[a] != knots.back()) {
            if (rootsD1.size() == 0)
                rootsD1.push_back(allRootsD1[a]);
            if (rootsD1.size() > 0)
                if (allRootsD1[a] != rootsD1.back())
                    rootsD1.push_back(allRootsD1[a]);
        }

    std::vector<double> xMaxima, yMaxima, xMinima, yMinima;

    if (rootsD1.size() > 0)
        for (int a=0; a<rootsD1.size(); ++a) {
            if (this->D2(rootsD1[a]) < 0) {
                xMaxima.push_back(rootsD1[a]);
                yMaxima.push_back(this->D0(rootsD1[a]));
            }
            if (this->D2(rootsD1[a]) >= 0) {
                xMinima.push_back(rootsD1[a]);
                yMinima.push_back(this->D0(rootsD1[a]));
            }
        }

    if (xMaxima.size() == 0) {
        hasOneMaximum = false;
        hasTwoMaxima = false;
        hasFirstOfTwoMaxima = false;
        return;
    }

    double smallDistance =
        yD0Range*fractionOfOrdinateRangeForMaximumIdentification;

    double yLeft = this->D0(knots[0]);
    double yRight = this->D0(knots.back());

    // Checks whether any maximum found is not a well-defined maximum

    std::vector<double> xMaximaWellDefined;
    std::vector<double> yMaximaWellDefined;

    if (xMinima.size() == 0) {
        xMaximaWellDefined = xMaxima;
        yMaximaWellDefined = yMaxima;
    }

    if (xMinima.size() > 0)
        for (int a=0; a<xMaxima.size(); ++a) {

            double yMinimumLeft = yLeft;
            double yMinimumRight = yRight;

            for (int b=0; b<xMinima.size(); ++b)
                if (xMinima[b] < xMaxima[a])
                    yMinimumLeft = yMinima[b];

            for (int b=xMinima.size()-1; b>-1; --b) {
                if (xMinima[b] > xMaxima[a])
                    yMinimumRight = yMinima[b];
            }

            if (fabs(yMaxima[a]-yMinimumLeft) > smallDistance &&
                fabs(yMaxima[a]-yMinimumRight) > smallDistance) {
                xMaximaWellDefined.push_back(xMaxima[a]);
                yMaximaWellDefined.push_back(yMaxima[a]);
            }

        }

    if (xMaximaWellDefined.size() == 0) {
        hasOneMaximum = false;
        hasTwoMaxima = false;
        hasFirstOfTwoMaxima = false;
        return;
    }

    // Checks whether any maximum found is part of an asymptote

    double asymptoteRange =
        yD0Range*fractionOfOrdinateRangeForAsymptoteIdentification;

    bool asymptoteLeft = false;

    if (xMinima.size() == 0)
        if (fabs(yMaximaWellDefined[0]-yLeft) < asymptoteRange)
            asymptoteLeft = true;
    if (xMinima.size() > 0) {
        if (xMaximaWellDefined[0] < xMinima[0]) {
            if (fabs(yMaximaWellDefined[0]-yLeft) < asymptoteRange)
                asymptoteLeft = true;
        }
        else {
            int minimaBeforeFirstWellDefinedMaximum = 0;
            int numberOfMinimaWithinAsymptoteRange = 0;
            for (int a=0; a<xMinima.size(); ++a)
                if (xMinima[a] < xMaximaWellDefined[0]) {
                    ++minimaBeforeFirstWellDefinedMaximum;
                    if (fabs(yMinima[a]-yLeft) < asymptoteRange)
                        ++numberOfMinimaWithinAsymptoteRange;
                }
            if (numberOfMinimaWithinAsymptoteRange ==
                minimaBeforeFirstWellDefinedMaximum)
                if (fabs(yMaximaWellDefined[0]-yLeft) < asymptoteRange)
                    asymptoteLeft = true;
        }
    }

    bool asymptoteRight = false;

    if (xMinima.size() == 0)
        if (fabs(yMaximaWellDefined.back()-yRight) < asymptoteRange)
            asymptoteRight = true;
    if (xMinima.size() > 0) {
        if (xMaximaWellDefined.back() > xMinima.back()) {
            if (fabs(yMaximaWellDefined.back()-yRight) < asymptoteRange)
                asymptoteRight = true;
        }
        else {
            int minimaAfterLastWellDefinedMaximum = 0;
            int numberOfMinimaWithinAsymptoteRange = 0;
            for (int a=xMinima.size()-1; a>-1; --a)
                if (xMinima[a] > xMaximaWellDefined.back()) {
                    ++minimaAfterLastWellDefinedMaximum;
                    if (fabs(yMinima[a]-yRight) < asymptoteRange)
                        ++numberOfMinimaWithinAsymptoteRange;
                }
            if (numberOfMinimaWithinAsymptoteRange ==
                minimaAfterLastWellDefinedMaximum)
                if (fabs(yMaximaWellDefined.back()-yRight) < asymptoteRange)
                    asymptoteLeft = true;
        }
    }

    std::vector<double> xMaximaActual;
    std::vector<double> yMaximaActual;

    for (int a=0; a<xMaximaWellDefined.size(); ++a) {

        bool actualMaximum = true;

        if (a == 0 && asymptoteLeft == true)
            actualMaximum = false;

        if (a == xMaximaWellDefined.size()-1 && asymptoteRight == true)
            actualMaximum = false;

        if (actualMaximum == true) {
            xMaximaActual.push_back(xMaximaWellDefined[a]);
            yMaximaActual.push_back(yMaximaWellDefined[a]);
        }

    }
    //for (int a=0; a<xMaximaWellDefined.size(); ++a) {
    //	std::cout<< " xMaximaActual " << a << " is " << xMaximaActual[a] << std::endl;
    //	std::cout<< " yMaximaActual " << a << " is " << yMaximaActual[a] << std::endl;
    //}
    if (xMaximaActual.size() == 0 || xMaximaActual.size() > 2) {
        hasOneMaximum = false;
        hasTwoMaxima = false;
        hasFirstOfTwoMaxima = false;
    }

    // Checks whether the maxima found are higher than the corresponding
    // extreme(s) of the spline, and whether, if there is only one maximum, it
    // is because there would have been two maxima had further experimental data
    // been collected

    std::vector<double> xMaximaAccepted;

    if (xMaximaActual.size() == 1)
        if (yMaximaActual[0] > yLeft)
            xMaximaAccepted.push_back(xMaximaActual[0]);

    if (xMaximaActual.size() == 2)
        if (yMaximaActual[0] > yLeft && yMaximaActual[1] > yRight) {
            xMaximaAccepted.push_back(xMaximaActual[0]);
            xMaximaAccepted.push_back(xMaximaActual[1]);
        }

    if (xMaximaAccepted.size() == 0) {
        hasOneMaximum = false;
        hasTwoMaxima = false;
        hasFirstOfTwoMaxima = false;
    }

    // Saves the abscissae of the maximum or of the maxima

    if (xMaximaAccepted.size() == 1) {

        if (yMaximaActual[0] > yLeft && yMaximaActual[0] > yRight) {
            hasOneMaximum = true;
            hasTwoMaxima = false;
            hasFirstOfTwoMaxima = false;
            xOfMax = xMaximaAccepted[0];
        }

        if (yMaximaActual[0] > yLeft && yMaximaActual[0] <= yRight) {
            hasOneMaximum = false;
            hasTwoMaxima = false;
            hasFirstOfTwoMaxima = true;
            xOfMax = xMaximaAccepted[0];
        }

    }

    if (xMaximaAccepted.size() == 2) {
        hasOneMaximum = false;
        hasTwoMaxima = true;
        hasFirstOfTwoMaxima = false;
        xOfFirstMax = xMaximaAccepted[0];
        xOfSecondMax = xMaximaAccepted[1];
    }

}
