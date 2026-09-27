#include <iostream>
#include <cmath>
#include <memory>

struct Complex {
    virtual double getRe() const = 0;
    virtual double getIm() const = 0;
    virtual double getModule() const = 0;
    virtual double getAngle() const = 0;
    virtual ~Complex() = default;

    virtual std::unique_ptr<Complex> add(const Complex& z) const = 0;
    virtual std::unique_ptr<Complex> sub(const Complex& z) const = 0;
    virtual std::unique_ptr<Complex> mult(const Complex& z) const = 0;
    virtual std::unique_ptr<Complex> div(const Complex& z) const = 0;
    virtual std::unique_ptr<Complex> negative() const = 0;
};

using ComplexPtr = std::unique_ptr<Complex>;

class AlgComplex : public Complex {
    private:
        double Re, Im;
    public:
        AlgComplex(double Re, double Im) : Re(Re), Im(Im) {}
        double getRe() const override { return Re;}
        double getIm() const override { return Im;} 
        
        double getModule() const override {
            return std::sqrt(Re * Re + Im * Im);
        }

        double getAngle() const override {
            return std::atan2(Im, Re);
        }

        ComplexPtr add(const Complex& z) const override {
            double newRe = Re + z.getRe();
            double newIm = Im + z.getIm();
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr sub(const Complex& z) const override {
            double newRe = Re - z.getRe();
            double newIm = Im - z.getIm();
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr mult(const Complex& z) const override {
            double newModule = getModule() * z.getModule();
            double newAngle = getAngle() + z.getAngle();
            double newRe = newModule * std::cos(newAngle);
            double newIm = newModule * std::sin(newAngle);
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr div(const Complex& z) const override {
            if (z.getModule() < 1e-12) {
        throw std::invalid_argument("Division by zero complex number");
    }
            double newModule = getModule() / z.getModule();
            double newAngle = getAngle() - z.getAngle();
            double newRe = newModule * std::cos(newAngle);
            double newIm = newModule * std::sin(newAngle);
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr negative() const override {
           return std::make_unique<AlgComplex>(-Re, -Im);
        }
};

class TrigComplex : public Complex {
    private:
        double r, phi;

        static double valueR(double r) {
            if (r < 0) {
                throw std::invalid_argument("r не может быть отрицательным");
            }
            return r;
        }
    public:
        TrigComplex(double r, double phi) : r(valueR(r)), phi(phi) {}
        double getModule() const override { return r;}
        double getAngle() const override { return phi;}

        double getRe() const override { 
            return r * std::cos(phi); 
        }
        double getIm() const override { 
            return r * std::sin(phi); 
        }   

        ComplexPtr add(const Complex& z) const override {
            double newRe = getRe() + z.getRe();
            double newIm = getIm() + z.getIm();
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr sub(const Complex& z) const override {
            double newRe = getRe() - z.getRe();
            double newIm = getIm() - z.getIm();
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr mult(const Complex& z) const override {
            double newModule = getModule() * z.getModule();
            double newAngle = getAngle() + z.getAngle();
            double newRe = newModule * std::cos(newAngle);
            double newIm = newModule * std::sin(newAngle);
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr div(const Complex& z) const override {
            if (z.getModule() == 0) {
                throw std::invalid_argument("Нельзя делить на ноль");
            }
            double newModule = getModule() / z.getModule();
            double newAngle = getAngle() - z.getAngle();
            double newRe = newModule * std::cos(newAngle);
            double newIm = newModule * std::sin(newAngle);
            return std::make_unique<AlgComplex>(newRe, newIm);
        }

        ComplexPtr negative() const override {
           return std::make_unique<AlgComplex>(-getRe(), -getIm());
        }
};

ComplexPtr operator+(const Complex& z1, const Complex& z2) {
    return z1.add(z2);
}

ComplexPtr operator-(const Complex& z1, const Complex& z2) {
    return z1.sub(z2);
}

ComplexPtr operator*(const Complex& z1, const Complex& z2) {
    return z1.mult(z2);
}

ComplexPtr operator/(const Complex& z1, const Complex& z2) {
    return z1.div(z2);
}

ComplexPtr operator-(const Complex& z) {
    return z.negative();
}

std::ostream& operator<<(std::ostream& os, const Complex& z) {
        os << z.getRe() << " + " << z.getIm() << "i";
        return os;
    } 

std::pair<ComplexPtr, ComplexPtr> sqrtComplex(const Complex& z) {
    double newModule = std::sqrt(z.getModule());
    double angle1 = z.getAngle() / 2;
    double angle2 = angle1 + std::acos(-1);
    
    auto root1 = std::make_unique<AlgComplex>(newModule * std::cos(angle1), newModule * std::sin(angle1));
    auto root2 = std::make_unique<AlgComplex>(newModule * std::cos(angle2), newModule * std::sin(angle2));

    return std::make_pair(std::move(root1), std::move(root2));
}

std::pair<ComplexPtr, ComplexPtr> solveQuadratic (const Complex& a, const Complex& b, const Complex& c) {
    auto quadro_b = b.mult(b);
    auto ac = a.mult(c);
    auto four_ac = ac->mult(*std::make_unique<AlgComplex>(4, 0));
    auto D = quadro_b->sub(*four_ac);

    auto roots = sqrtComplex(*D);
    auto& root1 = roots.first;
    auto& root2 = roots.second;

    auto neg_b = b.negative();
    auto den = a.mult(*std::make_unique<AlgComplex>(2, 0));
    
    auto num1 = neg_b->add(*root1);
    auto num2 = neg_b->add(*root2);

    auto x1 = num1->div(*den);
    auto x2 = num2->div(*den);
    
    return std::make_pair(std::move(x1), std::move(x2));
}

int main() {
    AlgComplex a(1, 0);
    AlgComplex b(2, 0);
    AlgComplex c(3, 0);    

    auto roots = solveQuadratic(a, b, c);

    std::cout << "x1 = " << *roots.first << std::endl;
    std::cout << "x2 = " << *roots.second << std::endl;
} 
