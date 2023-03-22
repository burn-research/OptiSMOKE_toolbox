class BasisFunction {

public:

    /* Coefficients of the basis function */
    vector<vector<double>> coeffD0;

    /* Coefficients of the first derivative of the basis function */
    vector<vector<double>> coeffD1;

    /* Coefficients of the second derivative of the basis function */
    vector<vector<double>> coeffD2;

    ////////////////////////////////////////////////////////////////////////////

    /* Calculates the coefficients of the basis function, and calculates the
    coefficients of its first and second derivative */
    void calculateCoefficients(int indexFirstKnot, const vector<double>& knots);

    /* Calculates the value of the basis function at position x on the x-axis */
    double D0(double x);

    /* Calculates the value of the first derivative of the basis function at
    position x on the x-axis */
    double D1(double x);

    /* Calculates the integral of the product of the second derivatives of the
    current BasisFunction and basisFunction. The basis functions must belong to
    the same spline */
    double integralOfProductD2(const BasisFunction& basisFunction);

////////////////////////////////////////////////////////////////////////////////

private:

    /* Knots of the basis function, including non-real ones */
    vector<double> knotsAll;

    /* Knots of the basis function, excluding non-real ones */
    vector<double> knotsReal;

    /* Index of the leftmost knot of the basis function in the knotsAll vector
    */
    int j;

    /* Index of the first real knot in the knotsAll vector */
    int indexOfFirstRealKnot;

    /* Maximum height the basis function might reach, if it were the leftmost or
    rightmost basis function of the spline */
    double maxHeight;

    /* Powers of the abscissa where the value of the basis function is to be
    calculated */
    vector<double> powers;

    ////////////////////////////////////////////////////////////////////////////

    /* Finds the knots in common, given two sets of knots */
    void findKnotsInCommon(const vector<double>& knotsAlpha,
                           const vector<double>& knotsBeta,
                           vector<double>& knotsInCommon);

    /* Sums the basisAlpha and basisBeta basis functions. The basis functions
    must be separated by a single knot */
    void sumBasisFunctions(vector<vector<double>> basisAlpha,
                           vector<vector<double>> basisBeta,
                           int degree /* of the basis functions */,
                           int index /* of the first knot */,
                           vector<vector<double>>& basisSum);

    // FOR TESTING PURPOSES ONLY:
    /* Saves the data for creating the graph of the basis function to a .txt
    file, and creates the .R script to create the graph */
    void printBasisR(int derivativeOrder);

    // FOR TESTING PURPOSES ONLY:
    /* Calculates the value of the second derivative of the basis function at
    position x on the x-axis */
    double D2(double x);

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void BasisFunction::calculateCoefficients(int indexFirstKnot,
                                          const vector<double>& knots) {

    j = indexFirstKnot;

    for (int i=j; i<=j+m; ++i)
        knotsAll.push_back(knots[i]);

    int numberOfPolynomials = knots.size() - 1;
    double firstRealKnot = knots[g];
    double lastRealKnot = knots[numberOfPolynomials-g];
    for (int i=0; i<=m; ++i)
        if (knotsAll[i] >= firstRealKnot)
            if (knotsAll[i] <= lastRealKnot)
                knotsReal.push_back(knotsAll[i]);

    indexOfFirstRealKnot = 0;
    for (int i=0; i<m; ++i)
        if (knotsAll[i] >= firstRealKnot) {
            indexOfFirstRealKnot = i;
            break;
        }

    maxHeight = (knots.back()-knots[0])/numberOfPolynomials*(double)m;

    // The 'pyramid' vector of vectors of vectors will contain the coefficients
    // of the basis function

    vector<vector<vector<vector<double>>>> pyramid;
    for (int a=0; a<g; ++a) {
        vector<vector<vector<double>>> pyramidLevel;
        auto basis = vector<vector<double>>(a+1,vector<double>(m,0));
        for (int b=0; b<m-a; ++b)
            pyramidLevel.push_back(basis);
        pyramid.push_back(pyramidLevel);
    }

    for (int b=0; b<m; ++b)
        pyramid[0][b][0][0] = maxHeight;

    for (int a=1; a<g; ++a)
        for (int b=0; b<m-a; ++b)
            sumBasisFunctions(pyramid[a-1][b],
                              pyramid[a-1][b+1],
                              a /*degree of the basis functions*/,
                              b /*index of the first knot*/,
                              pyramid[a][b]);

    // coefficientsD0 is equal to the top of 'pyramid'
    coeffD0 = vector<vector<double>>(m,vector<double>(m,0));
    sumBasisFunctions(pyramid[g-1][0],
                      pyramid[g-1][1],
                      g /*degree of the basis functions*/,
                      0 /*index of the first knot*/,
                      coeffD0);

    coeffD1 = vector<vector<double>>(m,vector<double>(m,0));
    for (int i=0; i<m; ++i)
        for (int a=1; a<m; ++a)
            coeffD1[i][a-1] = (double)a*coeffD0[i][a];

    coeffD2 = vector<vector<double>>(m,vector<double>(m,0));
    for (int i=0; i<m; ++i)
        for (int a=2; a<m; ++a)
            coeffD2[i][a-2] = (double)(a*(a-1))*coeffD0[i][a];

    powers = vector<double>(m,1);

}



double BasisFunction::D0(double x) {

    // If x is outside the knots or equal to the rightmost knot, the basis
    // function is equal to 0
    if (x < knotsAll[0] || x >= knotsAll[m]) return 0;

    int indexOfPolynomial = 0;
    for (int i=0; i<m; ++i)
        if (x < knotsAll[i+1]) {
            indexOfPolynomial = i;
            break;
        }

    // Calculates the powers of x
    for (int i=1; i<m; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D0(x)
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0[indexOfPolynomial][i]*powers[i];

    return y;

}



double BasisFunction::D1(double x) {

    // If x is outside the knots or equal to the rightmost knot, the basis
    // function is equal to 0
    if (x < knotsAll[0] || x >= knotsAll[m]) return 0;

    int indexOfPolynomial = 0;
    for (int i=0; i<m; ++i)
        if (x < knotsAll[i+1]) {
            indexOfPolynomial = i;
            break;
        }

    // Calculates the powers of x
    for (int i=1; i<g; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D1(x)
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1[indexOfPolynomial][i]*powers[i];

    return y;

}



double BasisFunction::integralOfProductD2(const BasisFunction& basisFunction) {

    // knotsInCommon will contain the real knots in common
    vector<double> knotsInCommon;

    findKnotsInCommon(knotsReal, basisFunction.knotsReal, knotsInCommon);

    int numberOfKnotsInCommon = knotsInCommon.size();

    // D2Basis1 and D2Basis2 will contain the coefficients of the second
    // derivatives of the basis functions for the segments between real knots in
    // common to both basis functions

    vector<vector<double>> D2Basis1;
    vector<vector<double>> D2Basis2;

    for (int a=indexOfFirstRealKnot; a<m; ++a)
        if (knotsAll[a] == knotsInCommon[0]) {
            for (int b=0; b<numberOfKnotsInCommon-1; ++b)
                D2Basis1.push_back(coeffD2[a+b]);
            break;
        }

    for (int a=basisFunction.indexOfFirstRealKnot; a<m; ++a)
        if (basisFunction.knotsAll[a] == knotsInCommon[0]) {
            for (int b=0; b<numberOfKnotsInCommon-1; ++b)
                D2Basis2.push_back(basisFunction.coeffD2[a+b]);
            break;
        }

    int numberOfPolynomialsInCommon = D2Basis1.size();
    int mTimesTwoMinusFive = m*2-5;

    // productD2 will contain the coefficients of the product of the second
    // derivatives of the basis functions
    auto productD2 =
        vector<vector<double>>(numberOfPolynomialsInCommon,
        vector<double>(mTimesTwoMinusFive,0));

    // For each segments calculates the coefficients of the product of the
    // second derivatives of the basis functions
    int mMinus2 = m-2;
    for (int i=0; i<numberOfPolynomialsInCommon; ++i)
        for (int a=0; a<mMinus2; ++a)
            for (int b=0; b<mMinus2; ++b)
                productD2[i][a+b] += D2Basis1[i][a] * D2Basis2[i][b];

    // Calculates the powers of each knot in common
    auto powersKnotsInCommon =
        vector<vector<double>>(numberOfKnotsInCommon,vector<double>(m*2-4,1));
    for (int i=0; i<numberOfKnotsInCommon; ++i)
        for (int a=1; a<=mTimesTwoMinusFive; ++a)
            powersKnotsInCommon[i][a] =
                powersKnotsInCommon[i][a-1]*knotsInCommon[i];

    // Calculates the integral, segment by segment
    double integral = 0;
    for (int i=0; i<numberOfPolynomialsInCommon; ++i)
        for (int a=mTimesTwoMinusFive; a>0; --a)
            integral +=
            (((productD2[i][a-1]/(double)a)*powersKnotsInCommon[i+1][a])-
            ((productD2[i][a-1]/(double)a)*powersKnotsInCommon[i][a]));

    return integral;

}



void BasisFunction::findKnotsInCommon(const vector<double>& knotsAlpha,
                                      const vector<double>& knotsBeta,
                                      vector<double>& knotsInCommon) {

    int knotsAlphaSize = knotsAlpha.size();
    int knotsBetaSize = knotsBeta.size();

    for (int a=0; a<knotsAlphaSize; ++a)
        for (int b=0; b<knotsBetaSize; ++b)
            if (knotsAlpha[a] == knotsBeta[b]) {
                if (knotsAlphaSize-a <= knotsBetaSize-b) {
                    for (int c=a; c<knotsAlphaSize; ++c)
                        knotsInCommon.push_back(knotsAlpha[c]);
                    return;
                }
                else {
                    for (int c=b; c<knotsBetaSize; ++c)
                        knotsInCommon.push_back(knotsBeta[c]);
                    return;
                }
            }

}



void BasisFunction::sumBasisFunctions(vector<vector<double>> basisAlpha,
                                      vector<vector<double>> basisBeta,
                                      int degree /* of the basis functions */,
                                      int index /* of the first knot */,
                                      vector<vector<double>>& basisSum) {

    int order = degree + 1;

    // basisAlpha * ( t - u_index ) / ( u_index+p - u_index )
    for (int a=0; a<degree; ++a) {
        // basisAlpha * t
        for (int b=degree; b>0; --b)
            basisAlpha[a][b] = basisAlpha[a][b-1];
        basisAlpha[a][0] = 0;
        // basisAlpha * -u_index
        for (int b=0; b<degree; ++b)
            basisAlpha[a][b] -= basisAlpha[a][b+1] * knotsAll[index];
    }
    // basisAlpha / ( u_index+p - u_index )
    for (int a=0; a<degree; ++a)
        for (int b=0; b<order; ++b)
            basisAlpha[a][b] /= (knotsAll[index+degree] - knotsAll[index]);

    // basisBeta * ( u_index+p+1 - t ) / ( u_index+p+1 - u_index+1 )
    for (int a=0; a<degree; ++a) {
        // basisBeta * -t
        for (int b=degree; b>0; --b)
            basisBeta[a][b] = -basisBeta[a][b-1];
        basisBeta[a][0] = 0;
        // basisBeta * u_index+p+1
        for (int b=0; b<degree; ++b)
            basisBeta[a][b] -= basisBeta[a][b+1] * knotsAll[index+order];
    }
    // basisBeta * ( u_index+p+1 - u_index+1 )
    for (int a=0; a<degree; ++a)
        for (int b=0; b<order; ++b)
            basisBeta[a][b] /= (knotsAll[index+order] - knotsAll[index+1]);

    // basisAlpha + basisBeta = basisSum
    for (int a=0; a<order; ++a) {
        if (a==0)
            basisSum[0] = basisAlpha[0];
        else if (a==degree)
            basisSum[degree] = basisBeta[degree-1];
        else
            for (int b=0; b<order; ++b)
                basisSum[a][b] = basisAlpha[a][b] + basisBeta[a-1][b];
    }

}



void BasisFunction::printBasisR(int derivativeOrder) {

    int numberPoints = 10000;

    double distance = (knotsAll.back()-knotsAll[0]) / (double)(numberPoints-1);

    // Saves the data for the plot to a .txt file

    string fileName = "DatiPlotR.txt";
    string fileNameWithPath = "./" + fileName;

    const char* fileNameWithPathChar = fileNameWithPath.c_str();
    ofstream file(fileNameWithPathChar);

    file << "x\ty" << endl;
    for (int i=0; i<numberPoints; ++i) {
        file << knotsAll[0]+(double)i*distance << "\t";
        if (derivativeOrder == 0)
            file << this->D0(knotsAll[0]+(double)i*distance) << endl;
        else if (derivativeOrder == 1)
            file << this->D1(knotsAll[0]+(double)i*distance) << endl;
        else if (derivativeOrder == 2)
            file << this->D2(knotsAll[0]+(double)i*distance) << endl;
    }
    file.close();

    // Writes the .R script to create the graph

    string scriptName = "Script.R";

    const char* scriptNameAndPathChar = scriptName.c_str();
    ofstream script(scriptNameAndPathChar);

    script << "Test = read.table(\"" << fileName
           << "\", sep=\"\\t\", header=TRUE)" << endl
           << "plot(Test[,1],Test[,2], cex=0.01, col=3)";

    script.close();

}



double BasisFunction::D2(double x) {

    // If x is outside the knots or equal to the rightmost knot, the basis
    // function is equal to 0
    if (x < knotsAll[0] || x >= knotsAll[m]) return 0;

    int indexOfPolynomial = 0;
    for (int i=0; i<m; ++i)
        if (x < knotsAll[i+1]) {
            indexOfPolynomial = i;
            break;
        }

    // Calculates the powers of x
    for (int i=1; i<m-2; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D2(x)
    double y = 0;
    for (int i=0; i<m-2; ++i)
        y += coeffD2[indexOfPolynomial][i]*powers[i];

    return y;

}