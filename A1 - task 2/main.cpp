#include "Polynomial.h"

int main()
{
    Polynomial p1 = Polynomial({ 1,-2.000000000001,1.000000000001 });
    Polynomial p2({ });     // Represents -1 + 4x

    //Polynomial sum = p1 + p2;
    //Polynomial difference = p1 - p2;
    //Polynomial product = p1 * p2;

    cout << "p1: " << p1 << endl;         // Expected: 3x^2 - 2x + 1
    //cout << "p2: " << p2 << endl;         // Expected: 4x - 1
    //cout << "p1 + p2: " << sum << endl;   // Expected: 3x^2 + 2x
    //cout << "p1 - p2: " << difference << endl; // Expected: 3x^2 - 6x + 2
    //cout << "p1 * p2: " << product << endl; // Expected: 12x^3 - 8x^2 + 10x - 1

    double x = 2;
    cout << "p1 evaluated at x = " << x << ": " << p1.evaluate(x) << endl;  // Expected: 13.25
    //cout << "p2 evaluated at x = " << x << ": " << p2.evaluate(x) << endl;  // Expected: 9

    Polynomial derivative = p1.derivative();
    cout << "Derivative of p1: " << derivative << endl; // Expected: 6x - 2

    Polynomial integral = p1.integral();
    cout << "Integral of p1: " << integral << endl; // Expected: 4x^2/2 - x + C -> 2x^2 - x

    // Corrected root finding line
    double root = p1.getRoot(1.1);  // Initial guess of 1.0 for finding root
    cout << "Root of p1 starting from 1.0: " << root << endl;  // Approximate root of p1

    // Composition: P(Q(x))
    Polynomial composed = p1.compose(p2);
    cout << "Composition of p1 and p2 (p1(p2(x))): " << composed << endl;

    // Definite integral of p1 from x = 0 to x = 1
    double integralValue = p1.definiteIntegral(0, 1);
    cout << "Definite integral of p1 from 0 to 1: " << integralValue << endl;

    return 0;
}
