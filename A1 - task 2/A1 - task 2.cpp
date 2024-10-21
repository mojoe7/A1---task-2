#include "Polynomial.h"

// Constructors
Polynomial::Polynomial() : coeffs(1, 0.0) {}

Polynomial::Polynomial(const vector<double>& coefficients) : coeffs(coefficients) {}

Polynomial::Polynomial(const Polynomial& other) : coeffs(other.coeffs) {}

// Destructor
Polynomial::~Polynomial() {}

// Assignment operator
Polynomial& Polynomial::operator=(const Polynomial& other) {
    if (this != &other) {
        coeffs = other.coeffs;
    }
    return *this;
}

// Arithmetic operators
Polynomial Polynomial::operator+(const Polynomial& other) const {
    vector<double> result(max(coeffs.size(), other.coeffs.size()), 0.0);
    for (size_t i = 0; i < coeffs.size(); ++i) result[i] += coeffs[i];
    for (size_t i = 0; i < other.coeffs.size(); ++i) result[i] += other.coeffs[i];
    return Polynomial(result);
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    vector<double> result(max(coeffs.size(), other.coeffs.size()), 0.0);
    for (size_t i = 0; i < coeffs.size(); ++i) result[i] += coeffs[i];
    for (size_t i = 0; i < other.coeffs.size(); ++i) result[i] -= other.coeffs[i];
    return Polynomial(result);
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    vector<double> result(coeffs.size() + other.coeffs.size() - 1, 0.0);
    for (size_t i = 0; i < coeffs.size(); ++i) {
        for (size_t j = 0; j < other.coeffs.size(); ++j) {
            result[i + j] += coeffs[i] * other.coeffs[j];
        }
    }
    return Polynomial(result);
}

// Equality operator
bool Polynomial::operator==(const Polynomial& other) const {
    return coeffs == other.coeffs;
}

// Output operator
ostream& operator<<(ostream& out, const Polynomial& poly) {
    for (int i = poly.coeffs.size() - 1; i >= 0; --i) {
        if (i < poly.coeffs.size() - 1 && poly.coeffs[i] >= 0) out << " + ";
        if (poly.coeffs[i] != 0) out << poly.coeffs[i] << "x^" << i;
    }
    return out;
}

// Utility functions
int Polynomial::degree() const {
    return coeffs.size() - 1;
}

double Polynomial::evaluate(double x) const {
    double result = 0;
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

Polynomial Polynomial::derivative() const {
    if (coeffs.size() <= 1) return Polynomial(vector<double>{0});
    vector<double> deriv(coeffs.size() - 1);
    for (size_t i = 1; i < coeffs.size(); ++i) {
        deriv[i - 1] = coeffs[i] * i;
    }
    return Polynomial(deriv);
}

Polynomial Polynomial::integral() const {
    vector<double> integral(coeffs.size() + 1);
    for (size_t i = 0; i < coeffs.size(); ++i) {
        integral[i + 1] = coeffs[i] / (i + 1);
    }
    return Polynomial(integral);
}

double Polynomial::integral(double x1, double x2) const {
    Polynomial integralPoly = integral();
    return integralPoly.evaluate(x2) - integralPoly.evaluate(x1);
}

double Polynomial::getRoot(double guess, double tolerance, int maxIter) {
    Polynomial deriv = derivative();
    double x = guess;
    for (int i = 0; i < maxIter; ++i) {
        double fx = evaluate(x);
        double fpx = deriv.evaluate(x);
        if (fabs(fx) < tolerance) return x;
        x = x - fx / fpx;
    }
    return x;
}

void Polynomial::setCoefficients(const vector<double>& coefficients) {
    coeffs = coefficients;
}

double Polynomial::getCoefficient(int degree) const {
    if (degree < 0 || degree >= coeffs.size()) return 0;
    return coeffs[degree];
}

// Method for composition
Polynomial Polynomial::compose(const Polynomial& other) const {
    vector<double> result;

    // Evaluate P(Q(x)) where P is the current polynomial and Q is the other polynomial
    for (int i = 0; i <= other.degree(); ++i) {
        double coeff = evaluate(other.getCoefficient(i)); // Evaluate the current polynomial at Q(x)
        // Add the evaluated coefficients to the resulting polynomial
        result.push_back(coeff);
    }

    return Polynomial(result);
}

// Method for definite integral
double Polynomial::definiteIntegral(double x1, double x2) const {
    return integral(x1, x2); // Reuse the existing integral method
}
