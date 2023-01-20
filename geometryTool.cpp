/**
*
* Solution to course project # 3
* Introduction to programming course
* Faculty of Mathematics and Informatics of Sofia University
* Winter semester 2022/2023
*
* @author Anastasia Madzharova
* @idnumber 0MI0600156
* @compiler GCC
*
* main file
*
*/

#include <iostream>
#include <cmath>
#include <stdlib.h>
double const EPSILON = 0.0001;
unsigned int const NAMESIZE = 333;

double myFAbs(double num)
{
    if (num < 0.0000) {
        return -num;
    }
    return num;
}

void mySwap(int& num1, int& num2)
{
    double temp = num1;
    num1 = num2;
    num2 = temp;
}

double gcd(int& num1, int& num2)
{
    if (num1 < num2) {
        mySwap(num1, num2);
    }
    while (num2!=0) {
        int mod = num1%num2;
        num1 = num2;
        num2 = mod;
    }
    return num1;
}

void fixCoefficients(double& slCoeff, double& yIntercept, double& thirdCoeff)
{
    if (slCoeff < 0) {
        slCoeff *= -1;
        yIntercept *= -1;
        thirdCoeff *= -1;
    }
    int intSlCoeff = slCoeff;
    int intyIntercept = yIntercept;
    int intThirdCoeff = thirdCoeff;
    if (myFAbs(slCoeff - intSlCoeff) < EPSILON && myFAbs(yIntercept - intyIntercept) < EPSILON && myFAbs(thirdCoeff - intThirdCoeff) < EPSILON) {
        int gcdFirstSecond = gcd(intSlCoeff, intyIntercept);
        double gcdAll = gcd(gcdFirstSecond, intThirdCoeff);
        slCoeff /= gcdAll;
        yIntercept /= gcdAll;
        thirdCoeff /= gcdAll;
    }
}

void checkInput(double& num)
{
    while (true) {
        std::cin >> num;
        if (std::cin.fail()) {
            std::cout << "Invalid input. Please input a number. " << std::endl;

            std::cin.clear();
            std::cin.ignore(100000, '\n');
        }
        else {
            break;
        }
    }
}

void checkInput(unsigned int& num)
{
    while (true) {
        std::cin >> num;

        if (std::cin.fail()) {
            std::cout << "Invalid input. Please input a positive integer. " << std::endl;
            std::cin.clear();
            std::cin.ignore(100000, '\n');
        }
        else {
            break;
        }
    }
}

bool areStringsEqual(const char* str1, char* str2)
{
    while (*str1 && *str2) {
        if (*str1 != *str2) {
            return false;
        }
        str1++;
        str2++;
    }
    if (*str1=='\0' && *str2=='\0') {
        return true;
    }
    return false;
}

bool pointIsOnALine(double x, double y, double slopeCoeff, double yIntercept, double constant)
{
    return (myFAbs(slopeCoeff * x + yIntercept * y + constant) < EPSILON);
}

double thirdCoefficientOfALineGivenByASlopeCoefficientAndPoint(double x, double y, double lineSlopeCoefficient, double yIntercept)
{
    return -y*yIntercept - lineSlopeCoefficient*x;
}

void equationOfPerpendicularLinePassingThroughAPoint(double x, double y, double slopeCoefficient, double yIntercept, double lineThirdCoefficient, double& newLineSlopeCoefficient, double& newLineYIntercept, double& newLineThirdCoefficient)
{
    if (myFAbs(slopeCoefficient)>EPSILON) {
        if (myFAbs(yIntercept )< EPSILON) {
            newLineSlopeCoefficient = 0;
            newLineYIntercept = 1;
            newLineThirdCoefficient = -y;
        }
        else {
            newLineSlopeCoefficient = -yIntercept / slopeCoefficient;
            newLineYIntercept = 1;
            newLineThirdCoefficient = thirdCoefficientOfALineGivenByASlopeCoefficientAndPoint(x, y, newLineSlopeCoefficient, newLineYIntercept);
        }
    }
    else {
        newLineYIntercept = 0;
        newLineSlopeCoefficient = 1;
        newLineThirdCoefficient = -x;
    }
}

void intersectionOfTwoLines(double slCoeff1, double yIntercept1, double thirdCoeff1, double slCoeff2, double yIntercept2, double thirdCoeff2, bool& noIntersection, bool& linesMatch,double& pointX, double& pointY )
{
    if (myFAbs(yIntercept1)<EPSILON&& myFAbs(yIntercept2) < EPSILON) {
        if (myFAbs(thirdCoeff1 / slCoeff1 - thirdCoeff2 / slCoeff2) < EPSILON) {
            linesMatch = true;
        }
        else {
            noIntersection = true;
        }
    }
    else if (myFAbs(slCoeff1) < EPSILON&& myFAbs(slCoeff2)<EPSILON) {
        if (myFAbs(thirdCoeff1 / yIntercept1 - thirdCoeff2 / yIntercept2) < EPSILON) {
            linesMatch = true;
        }
        else {
            noIntersection = true;
        }
    }
    else if (myFAbs(slCoeff1) < EPSILON && myFAbs(yIntercept2) < EPSILON) {
        pointX = -thirdCoeff2 / slCoeff2;
        pointY = -thirdCoeff1 / yIntercept1;
    }
    else if (myFAbs(slCoeff2) <EPSILON&& myFAbs(yIntercept1) <EPSILON) {
        pointX = -thirdCoeff1 / slCoeff1;
        pointY = -thirdCoeff2 / yIntercept2;
    }
    else if (myFAbs(yIntercept1) <EPSILON) {
        pointX = -thirdCoeff1 / slCoeff1;
        pointY = -thirdCoeff2 / -slCoeff2 * pointX;
    }
    else if (myFAbs(yIntercept2) <EPSILON) {
        pointX = -thirdCoeff2 / slCoeff2;
        pointY = -thirdCoeff1 -slCoeff1 * pointX;
    }
    else if (myFAbs(slCoeff1) <EPSILON) {
        pointY = -thirdCoeff1 / yIntercept1;
        pointX = -thirdCoeff2 - yIntercept2 * pointY;
    }
    else {
        if (myFAbs(slCoeff1 / slCoeff2 - yIntercept1 / yIntercept2) < EPSILON) {
            if (myFAbs(slCoeff1 / slCoeff2 - thirdCoeff1 / thirdCoeff2) < EPSILON) {
                linesMatch = true;
            }
            else {
                noIntersection = true;
                //parallel lines
            }
        }
        else {
            slCoeff2 *= (-yIntercept1);
            yIntercept2 *= (-yIntercept1);
            thirdCoeff2 *= (-yIntercept1);
            slCoeff1 *= yIntercept2 / (-yIntercept1);
            thirdCoeff1*= yIntercept2 / (-yIntercept1);
            yIntercept1 *= yIntercept2 / (-yIntercept1);
            slCoeff1 += slCoeff2;
            thirdCoeff1 += thirdCoeff2;
            pointX = -thirdCoeff1 / slCoeff1;
            pointY = (-thirdCoeff2 - slCoeff2 * pointX) / yIntercept2;
        }
    }
}

void equationOfALineGivenByTwoPoints(double point1X, double point1Y, double point2X, double point2Y, double& lineSlopeCoefficient, double& lineYIntercept, double& lineThirdCoefficient)
{
        lineYIntercept = point2X - point1X;
        lineSlopeCoefficient = point1Y - point2Y;
        lineThirdCoefficient = (point1X - point2X) * point1Y + (point2Y - point1Y) * point1X;
}

void equationOfAltitudeOfATriangle(double AX,  double AY,double BX, double BY, double CX, double CY, double& altitudeSlopeCoefficient, double& altitudeYCoefficient, double& altitudeThirdCoefficient)
{
    //apex-C
    double lineABSlopeCoefficient = 0;
    double lineABYCoefficient=0;
    double lineABThirdCoefficient=0;
    equationOfALineGivenByTwoPoints(AX, AY, BX, BY, lineABSlopeCoefficient, lineABYCoefficient, lineABThirdCoefficient);
    equationOfPerpendicularLinePassingThroughAPoint(CX, CY, lineABSlopeCoefficient, lineABYCoefficient, lineABThirdCoefficient, altitudeSlopeCoefficient, altitudeYCoefficient, altitudeThirdCoefficient);
}

void coordinatesOfAMiddlePointOfASegment(double x1, double y1, double x2, double y2, double& XMiddle, double& YMiddle)
{
    XMiddle = (x1 + x2) / 2;
    YMiddle = (y1 + y2) / 2;
}

void equationOfAMedianOfATriangle(double AX, double AY, double BX, double BY, double CX, double CY, double& medianSlopeCoefficient, double& medianYCoefficient, double& medianThirdCoefficient)
{
    //apex-C
    double XMiddle;
    double YMiddle;
    coordinatesOfAMiddlePointOfASegment(AX, AY, BX, BY, XMiddle, YMiddle);
    equationOfALineGivenByTwoPoints(CX, CY, XMiddle, YMiddle, medianSlopeCoefficient, medianYCoefficient, medianThirdCoefficient);
}

void equationOfBisectorOfATriangle(double AX, double AY, double BX, double BY, double CX, double CY, double& bisectorSlopeCoefficient, double& bisectorYCoefficient, double& bisectorThirdCoefficient)
{
    //apex-C
    double XMiddle=0;
    double YMiddle=0;
    coordinatesOfAMiddlePointOfASegment(AX, AY, BX, BY, XMiddle, YMiddle);
    double ABSlopeCoefficient =0 ;
    double ABYCoefficient=0;
    double ABThirdCoefficient=0;
    equationOfALineGivenByTwoPoints(AX, AY, BX, BY, ABSlopeCoefficient, ABYCoefficient, ABThirdCoefficient);
    equationOfPerpendicularLinePassingThroughAPoint(XMiddle, YMiddle, ABSlopeCoefficient, ABYCoefficient, ABThirdCoefficient, bisectorSlopeCoefficient, bisectorYCoefficient, bisectorThirdCoefficient);
}

void derivativeOfQuadraticEquation(double x2Coeff, double xCoeff, double constant, double& derivativeXCoeff, double& derivativeConst)
{
    derivativeXCoeff = 2 * x2Coeff;
    derivativeConst = xCoeff;
}

void tangentToParabolaAtPoint(double x1, double y1, double x2Coeff, double xCoeff, double constant, double& tanSlopeCoeff, double& tanYCoeff, double& tanThirdCoeff)
{
    tanYCoeff = -1;
    double derivativeXCoeff;
    double derivativeConst;
    derivativeOfQuadraticEquation(x2Coeff, xCoeff, constant, derivativeXCoeff, derivativeConst);
    double derValueAtPoint = derivativeXCoeff * x1 + derivativeConst;
    double quadrEqValueAtPoint = x2Coeff * x1 * x1 + xCoeff * x1 + constant;
    tanSlopeCoeff = derValueAtPoint;
    tanThirdCoeff = quadrEqValueAtPoint - x1 * derValueAtPoint;
}

void tangentToParabolaFromPoint(double a, double b, double c, double x1, double y1, double& slope1, double& yIntercept1, double& thirdCoeff1, double& slope2, double& yIntercept2, double& thirdCoeff2, unsigned int& linesCount)
{
    double discriminant = a * a * x1 * x1 + a * b * x1 + c;
    double sqrtDiscriminant = sqrt(discriminant);
    double tanPoint1X = (a * x1 - sqrtDiscriminant) / a;
    double tanPoint1Y = -2.0 * x1 * ( sqrtDiscriminant -b) - b * sqrtDiscriminant / a + 2.0 * a * x1 * x1+2.0*c;
    double tanPoint2X = (sqrtDiscriminant + a * x1) / a;
    double tanPoint2Y = 2.0 * x1 * (sqrtDiscriminant + b) + b * sqrtDiscriminant / a + 2.0 * a * x1 * x1 + 2.0 * c;
    if (myFAbs(a * x1 * x1 + b * x1 + c) < EPSILON) {
        tangentToParabolaAtPoint(x1, y1, a, b, c, slope1, yIntercept1, thirdCoeff1);
        linesCount = 1;
    }
    else {
        slope1 = 2.0 * a * tanPoint1X + b;
        thirdCoeff1 = tanPoint1Y - tanPoint1X * slope1;
        yIntercept1 = -1;
        slope2 = 2.0 * a * tanPoint2X + b;
        thirdCoeff2 = tanPoint2Y - tanPoint2X * slope2;
        yIntercept2 = -1;
        linesCount = 2;
    }
}

void intersectionCheck(const double* linesCoefficients, double* linesIntersectionPoints, bool* intersections, unsigned int numLine1, unsigned int numLine2, unsigned int& line1IntersectionPoints, unsigned int& line2IntersectionPoints, bool& linesMatch, bool noIntersection, unsigned int boolIndex, unsigned int intersPointIndex)
{
    unsigned int line1StartCoeff = numLine1 * 3 - 3;
    unsigned int line2StartCoeff = numLine2 * 3 - 3;
    intersectionOfTwoLines(linesCoefficients[line1StartCoeff], linesCoefficients[line1StartCoeff + 1], linesCoefficients[line1StartCoeff + 2], linesCoefficients[line2StartCoeff], linesCoefficients[line2StartCoeff + 1], linesCoefficients[line2StartCoeff + 2], noIntersection, linesMatch, linesIntersectionPoints[intersPointIndex], linesIntersectionPoints[intersPointIndex + 1]);
    if (!noIntersection) {
        intersections[boolIndex] = true;
        line1IntersectionPoints++;
        line2IntersectionPoints++;
    }
}

bool isQuadrilateral(const double* linesCoefficients, double* linesIntersectionPoints, bool* intersections, unsigned int& line1IntersectionPoints, unsigned int& line2IntersectionPoints, unsigned int& line3IntersectionPoints, unsigned int& line4IntersectionPoints)
{
    bool noIntersection = false;
    bool linesMatch = false;
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 1, 2, line1IntersectionPoints, line2IntersectionPoints, linesMatch, noIntersection, 0, 0);
    if (linesMatch) {
        return false;
    }
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 1, 3, line1IntersectionPoints, line3IntersectionPoints, linesMatch, noIntersection, 1, 2);
    if (linesMatch) {
        return false;
    }
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 1, 4, line1IntersectionPoints, line4IntersectionPoints, linesMatch, noIntersection, 2, 4);
    if (linesMatch) {
        return false;
    }
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 2, 3, line2IntersectionPoints, line3IntersectionPoints, linesMatch, noIntersection, 3, 6);
    if (linesMatch) {
        return false;
    }
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 2, 4, line2IntersectionPoints, line4IntersectionPoints, linesMatch, noIntersection, 4, 8);
    if (linesMatch) {
        return false;
    }
    intersectionCheck(linesCoefficients, linesIntersectionPoints, intersections, 3, 4, line3IntersectionPoints, line4IntersectionPoints, linesMatch, noIntersection, 5, 10);
    if (linesMatch) {
        return false;
    }
    if (line1IntersectionPoints < 2 || line2IntersectionPoints < 2 || line3IntersectionPoints < 2 || line4IntersectionPoints < 2) {
        return false;
    }
    if (line1IntersectionPoints == 3 && line2IntersectionPoints == 3 && line3IntersectionPoints == 3 && line4IntersectionPoints == 3) {
        unsigned int equal = 0;
        for (int i = 1; i < 7; i+=2) {
            for (int m = i+1; m < 11; m+=2) {
                if (myFAbs(linesIntersectionPoints[i] - linesIntersectionPoints[m + 1]) < EPSILON && myFAbs(linesIntersectionPoints[m] - linesIntersectionPoints[i - 1]) < EPSILON) {
                    equal += 1;
                    if (equal == 3)
                        return false;
                }
            }
            equal = 0;
        }
    }
    return true;
}

bool areLinesPerpendicular(double slCoeff1, double yIntr1, double slCoeff2, double yIntr2)
{
    if (myFAbs(slCoeff1) > EPSILON) {
        if (myFAbs(slCoeff2) > EPSILON) {
            return myFAbs(slCoeff1 + slCoeff2) < EPSILON;
        }
        else {
            return myFAbs(yIntr1) < EPSILON;
        }
    }
    else {
        return myFAbs(yIntr2) < EPSILON;
    }
}

std::string findTypeOfQuadrilateral(double* linesCoefficients, double* intrPoints, bool* intersections, unsigned int line1IntrPoints, unsigned int line2IntrPoints, unsigned int line3IntrPoints, unsigned int line4IntrPoints)
{
    if ((line1IntrPoints + line2IntrPoints + line3IntrPoints + line4IntrPoints) == 10) {
        return "trapezium";
    }
    else if ((line1IntrPoints + line2IntrPoints + line3IntrPoints + line4IntrPoints) == 8) {
        //check which line is parallel to line1
        if (intersections[0] && intersections[1]) {
            //1 and 4 are parallel
            //2 and 3 are parallel
            double line1Length = sqrt((intrPoints[0] - intrPoints[2]) * (intrPoints[0] - intrPoints[2]) + (intrPoints[1] - intrPoints[3]) * (intrPoints[1] - intrPoints[3]));
            double line2Length = sqrt((intrPoints[8] - intrPoints[0]) * (intrPoints[8] - intrPoints[0]) + (intrPoints[9] - intrPoints[1]) * (intrPoints[9] - intrPoints[1]));
            if (myFAbs(line1Length - line2Length) < EPSILON) {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[3], linesCoefficients[4])) {
                    return "square";
                }
                return "rhombus";
            }
            else {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[3], linesCoefficients[4])) {
                    return "rectangle";
                }
                return "parallelogram";
            }
        }
        else if (intersections[0] && intersections[2]) {
            //1 and 3 are parallel
            //2 and 4 are parallel
            double line1Length = sqrt((intrPoints[0] - intrPoints[6]) * (intrPoints[0] - intrPoints[6]) + (intrPoints[1] - intrPoints[7]) * (intrPoints[1] - intrPoints[7]));
            double line2Length = sqrt((intrPoints[0] - intrPoints[4]) * (intrPoints[0] - intrPoints[4]) + (intrPoints[1] - intrPoints[5]) * (intrPoints[1] - intrPoints[5]));
            if (myFAbs(line1Length - line2Length) < EPSILON) {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[3], linesCoefficients[4])) {
                    return "square";
                }
                return "rhombus";
            }
            else {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[3], linesCoefficients[4])) {
                    return "rectangle";
                }
                return "parallelogram";
            }
        }
        else {
            //1 and 2 are parallel
            //3 and 4 are parallel
            double line1Length = sqrt((intrPoints[4] - intrPoints[2]) * (intrPoints[4] - intrPoints[2]) + (intrPoints[5] - intrPoints[3]) * (intrPoints[5] - intrPoints[3]));
            double line3Length = sqrt((intrPoints[2] - intrPoints[6]) * (intrPoints[2] - intrPoints[6]) + (intrPoints[3] - intrPoints[7]) * (intrPoints[3] - intrPoints[7]));
            if (myFAbs(line1Length - line3Length) < EPSILON) {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[6], linesCoefficients[7])) {
                    return "square";
                }
                return "rhombus";
            }
            else {
                if (areLinesPerpendicular(linesCoefficients[0], linesCoefficients[1], linesCoefficients[6], linesCoefficients[7])) {
                    return "rectangle";
                }
                return "parallelogram";
            }

        }
    }
    else {
        return "irregular quadrilateral";
    }
}

void intersectionPointsLineAndParabola(double parA, double parB, double parC, double lineSlCoeff, double lineYInt, double lineThirdCoeff,double& intersectionPoint1X, double& intersectionPoint1Y, double& intersectionPoint2X, double& intersectionPoint2Y, unsigned int& pointsCount)
{
    if (myFAbs(lineYInt) > EPSILON) {
        lineSlCoeff /= -lineYInt;
        lineThirdCoeff /= -lineYInt;
        parB -= lineSlCoeff;
        parC -= lineThirdCoeff;
        double discriminant = sqrt(parB * parB) - 4 * parA * parC;
    
        if (myFAbs(discriminant) <EPSILON) {
            pointsCount = 1;
            intersectionPoint1X = -parB / (2 * parA);
            intersectionPoint1Y = intersectionPoint1X * lineSlCoeff + lineYInt;
        }
        if (discriminant<0) {
            pointsCount = 0;
        }
        else {
            pointsCount = 2;
            intersectionPoint1X = (-parB + sqrt(discriminant)) / 2 * parA;
            intersectionPoint1Y = intersectionPoint1X * lineSlCoeff + lineYInt;
            intersectionPoint2X = (-parB - sqrt(discriminant)) / 2 * parA;
            intersectionPoint2Y = intersectionPoint2X * lineSlCoeff + lineYInt;
        }
    }
    else {
        pointsCount = 1;
        intersectionPoint1X = -lineThirdCoeff/lineSlCoeff;
        intersectionPoint1Y = parA * intersectionPoint1X * intersectionPoint1X + parB * intersectionPoint1X + parC;
    }
}

bool isNameValid(const char* name)
{
    while (*name)
    {
        bool currentSymbolValid = false;
        if ((*name >= '0' && *name <= '9') || (*name >= 'a' && *name <= 'z') || (*name >= 'A' && *name <= 'Z') || *name == '_')
        {
            currentSymbolValid = true;
        }
        if (!currentSymbolValid)
        {
            return false;
        }
        name++;
    }
    return true;
}

void nameOfAlineOrPoint(char* name)
{
    std::cin >> name;
    if (!isNameValid(name))
    {
        std::cout << "Invalid name. Please enter new name";
        std::cin >> name;
    }
}

void defineALine(char* name,double& slopeCoeff, double& yIntercept, double& constant)
{
    std::cout << "Please enter the coefficients A,B,C in the equation A*x + B*y + C = 0 to define your line." << std::endl;
    std::cout << "A = ";
    checkInput(slopeCoeff);
    std::cout << std::endl;
    std::cout << "B = ";
    checkInput(yIntercept);
    std::cout << std::endl;
    std::cout << "C = ";
    checkInput(constant);
    std::cout << std::endl;
    std::cout << "Do you want to enter a name for your line? yes/no" << std::endl;
    char answer[12] = "an";
    char answerYes[4] = "yes";
    char answerNo[3] = "no";
    while (!areStringsEqual(answer, answerYes) && !areStringsEqual(answer, answerNo)) {
        std::cin >> answer;
        if (areStringsEqual(answer, answerYes)) {
            nameOfAlineOrPoint(name);
        }
        else if (areStringsEqual(answer, answerNo)) {
            break;
        }
        else {
            std::cout << "Please enter valid answer: yes/no";
        }
    }
}

void defineAPoint(char* name, double& x, double& y)
{
    std::cout << "Please enter x and y to define a point." << std::endl;
    std::cout << "x = ";
    checkInput(x);
    std::cout << std::endl;
    std::cout << "y = ";
    checkInput(y);
    std::cout << std::endl;
    std::cout << "Do you want to enter a name for your point? yes/no " << std::endl;
     char answer[12] = "an";
     char answerYes[4] = "yes";
     char answerNo[3] = "no";
     while (!areStringsEqual(answer, answerYes) && !areStringsEqual(answer, answerNo)) {
         std::cin >> answer;
         if (areStringsEqual(answer, answerYes)) {
             nameOfAlineOrPoint(name);
         }
         else if (areStringsEqual(answer, answerNo)) {
             break;
         }
         else {
             std::cout << "Please enter valid answer: yes/no";
         }
     }
}

void printALine(double slopeCoeff, double yIntercept, double constant)
{
    fixCoefficients(slopeCoeff, yIntercept, constant);
    if (myFAbs(slopeCoeff) < EPSILON) {
        if (constant > 0) {
            std::cout << yIntercept << "y + " << constant << " = 0." << std::endl;
        }
        else if (constant < EPSILON) {
            std::cout << yIntercept << "y " << constant << " = 0." << std::endl;
        }
        else {
            std::cout << yIntercept << "y = 0." << std::endl;
        }
    }
    else if (myFAbs(yIntercept)<EPSILON) {
        if (constant >0) {
            std::cout << slopeCoeff << "x + " << constant << " = 0." << std::endl;
        }
        else if (constant<0) {
            std::cout <<slopeCoeff << "x " << constant << " = 0." << std::endl;
        }
        else {
            std::cout  << slopeCoeff << "x = 0." << std::endl;
        }
    }
    else {
        if (constant > 0) {
            if (yIntercept > 0) {
                std::cout << slopeCoeff << "x + " << yIntercept << "y + " << constant << " = 0." << std::endl;
            }
            else {
                std::cout <<  slopeCoeff << "x " << yIntercept << "y + " << constant << " = 0." << std::endl;
            }
        }
        else if (constant < 0) {
            if (yIntercept > 0) {
                std::cout <<  slopeCoeff << "x + " << yIntercept << "y " << constant << " = 0." << std::endl;
            }
            else {
                std::cout << slopeCoeff << "x " << yIntercept << "y " << constant << " = 0." << std::endl;
            }
        }
        else {
            if (yIntercept > 0) {
                std::cout << slopeCoeff << "x + " << yIntercept << "y = 0." << std::endl;
            }
            else {
                std::cout << slopeCoeff << "x " << yIntercept << "y = 0." << std::endl;
            }
        }
    }
}

void printApoint(double x, double y)
{
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
}

void printMessageChoice1()
{
    double slopeCoeff=0, yIntercept = 0, constant = 0;
    char name1[333] = "g";
    defineALine(name1,slopeCoeff, yIntercept, constant);
    while (slopeCoeff == 0 && yIntercept == 0) {
        std::cout << "At least one of the coefficients A and B must be non-zero. Please try again.";
        defineALine(name1,slopeCoeff, yIntercept, constant);
    }
    double x = 0, y =0;
    char name2[333] = "p";
    defineAPoint(name2,x, y);
    if (pointIsOnALine(x, y, slopeCoeff, yIntercept, constant)) {
        std::cout << "The point "<<name2<<" is on the line "<<name1 << std::endl;
    }
    else {
        std::cout << "The point "<<name2<<" is not on the line "<<name1 << std::endl;
    }
}

void printMessageChoice2()
{
    double slopeCoeff = 0, yIntercept = 0, constant = 0;
    char name1[333] = "g";
    defineALine(name1, slopeCoeff, yIntercept, constant);
    while (slopeCoeff == 0 && yIntercept == 0) {
        std::cout << "At least one of the coefficients A and B must be non-zero. Please try again.";
        defineALine(name1,slopeCoeff, yIntercept, constant);
    }
    double x = 0, y=0;
    char name2[333] = "p";
    defineAPoint(name2, x, y);
    double thirdCoeff =thirdCoefficientOfALineGivenByASlopeCoefficientAndPoint(x, y, slopeCoeff, yIntercept);
    std::cout << "The equation of the parallel line is: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
}

void printMessageChoice3()
{
    double slopeCoeff = 0, yIntercept=0, constant = 0;
    char name1[NAMESIZE] = "g";
    defineALine(name1,slopeCoeff, yIntercept, constant);
    while (slopeCoeff == 0 && yIntercept == 0)
    {
        std::cout << "At least one of the coefficients A and B must be non-zero. Please try again."<<std::endl;
        defineALine(name1,slopeCoeff, yIntercept, constant);
    }
    double x = 0, y=0;
    char name2[NAMESIZE] = "g";
    defineAPoint(name2, x, y);
   while(!pointIsOnALine(x, y, slopeCoeff, yIntercept, constant))
    {
        std::cout << "The point is not on the line. Plese try again." << std::endl;
        defineAPoint(name2,x, y);
    }
    double perpLineSlopeCoeff = 0, perpLineYIntercept = 0, perpLineConst = 0;
    equationOfPerpendicularLinePassingThroughAPoint(x, y, slopeCoeff, yIntercept, constant, perpLineSlopeCoeff, perpLineYIntercept, perpLineConst);
    std::cout << "The equation of the line, which is perpendicular to g with a fifth at p is: ";
    printALine(perpLineSlopeCoeff, perpLineYIntercept, perpLineConst);
}

void printMessageChoice4()
{
    double slopeCoeff1 = 0,yIntercept1 = 0,constant1 = 0;
    char name1[NAMESIZE] = "g";
    defineALine(name1, slopeCoeff1, yIntercept1, constant1);
    while (slopeCoeff1 == 0 && yIntercept1 == 0)
    {
        std::cout << "At least one of the coefficients A and B must be non-zero. Please try again." << std::endl;
        defineALine(name1, slopeCoeff1, yIntercept1, constant1);
    }
    double slopeCoeff2 = 0, yIntercept2 = 0, constant2=0;
    char name2[NAMESIZE] = "p";
    defineALine(name2,slopeCoeff2, yIntercept2, constant2);
    while (slopeCoeff2 == 0 && yIntercept2 == 0)
    {
        std::cout << "At least one of the coefficients A and B must be non-zero. Please try again." << std::endl;
        defineALine(name2,slopeCoeff2, yIntercept2, constant2);
    }
    bool noIntersection = false;
    bool linesMatch = false;
    double interY = 0, interX=0;
    intersectionOfTwoLines(slopeCoeff1, yIntercept1, constant1, slopeCoeff2, yIntercept2, constant2, noIntersection, linesMatch, interX, interY);
    if (noIntersection)
    {
        std::cout << "No intersection point between lines "<< name1<<" and "<<name2<<" ." << std::endl;
    }
    else if (linesMatch)
    {
        std::cout << "Lines "<<name1<<" and "<<name2<<" match." << std::endl;
    }
    else
    {
        std::cout << "The intersection point of "<<name1<< " and "<<name2<< " is x= " << interX << " , y= " << interY << std::endl;
    }
}

bool isValidTriangle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    if ((myFAbs(x1 - x2) < EPSILON && myFAbs(y1 - y2) < EPSILON) || (myFAbs(x1 - x3) < EPSILON && myFAbs(y1 - y3) < EPSILON) || (myFAbs(x2 - x3) < EPSILON && myFAbs(y2 - y3) < EPSILON)) {
        return false;
    }
    double lineSlopeCoeff = 0;
    double lineYIntercept = 0;
    double lineThirdCoeff = 0;
    equationOfALineGivenByTwoPoints(x1, y1, x2, y2, lineSlopeCoeff, lineYIntercept, lineThirdCoeff);
    if (pointIsOnALine(x3, y3, lineSlopeCoeff, lineYIntercept, lineThirdCoeff)) {
        return false;
    }
    equationOfALineGivenByTwoPoints(x2, y2, x3, y3, lineSlopeCoeff, lineYIntercept, lineThirdCoeff);
    if (pointIsOnALine(x1, y1, lineSlopeCoeff, lineYIntercept, lineThirdCoeff)) {
        return false;
    }
    equationOfALineGivenByTwoPoints(x1, y1, x3, y3, lineSlopeCoeff, lineYIntercept, lineThirdCoeff);
    if (pointIsOnALine(x2, y2, lineSlopeCoeff, lineYIntercept, lineThirdCoeff)) {
        return false;
    }
    return true;
}

void printEquationsOfAltitudes(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double slopeCoeff=0, yIntercept=0, thirdCoeff = 0;
    equationOfAltitudeOfATriangle(x1, y1, x2, y2, x3, y3, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "hc: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfAltitudeOfATriangle(x3, y3, x2, y2, x1, y1, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "ha: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfAltitudeOfATriangle(x3, y3, x1, y1, x2, y2, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "hb: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
}

void printEquationsOfMedians(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double slopeCoeff = 0, yIntercept = 0, thirdCoeff = 0;
    equationOfAMedianOfATriangle(x1, y1, x2, y2, x3, y3, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "mc: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfAMedianOfATriangle(x3, y3, x2, y2, x1, y1, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "ma: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfAMedianOfATriangle(x1, y1, x3, y3, x2, y2, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "mb: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    
}

void printEquationsOfBisectors(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double slopeCoeff = 0, yIntercept = 0, thirdCoeff = 0;
    equationOfBisectorOfATriangle(x1, y1, x2, y2, x3, y3, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "sc: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfBisectorOfATriangle(x3, y3, x2, y2, x1, y1, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "sa: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
    equationOfBisectorOfATriangle(x3, y3, x1, y1, x2, y2, slopeCoeff, yIntercept, thirdCoeff);
    std::cout << "sa: ";
    printALine(slopeCoeff, yIntercept, thirdCoeff);
}

void printMessageChoice5()
{
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0, x3 = 0, y3= 0;
    char name1[NAMESIZE] = "g";
    defineAPoint(name1,x1, y1);
    char name2[NAMESIZE] = "p";
    defineAPoint(name2, x2, y2);
    char name3[NAMESIZE] = "e";
    defineAPoint(name3, x3, y3);
    while (!isValidTriangle)
    {
        std::cout << "This is not valid triangle. Please enter new points." << std::endl;
        defineAPoint(name1,x1, y1);
        defineAPoint(name2, x2, y2);
        defineAPoint(name3, x3, y3);
    }
    printEquationsOfAltitudes(x1, y1, x2, y2, x3, y3);
    printEquationsOfMedians(x1, y1, x2, y2, x3, y3);
    printEquationsOfBisectors(x1, y1, x2, y2, x3, y3);
}

void defineAParabola(double& a, double& b, double& c)
{
    std::cout << "To define a parabola, enter the coefficients a, b and c (ax2 + bx + c = 0)" << std::endl;
    std::cout << "a= ";
    checkInput(a);
    std::cout << std::endl;
    std::cout << "b= ";
    checkInput(b);
    std::cout << std::endl;
    std::cout << "c= ";
    checkInput(c);
}

void printMessageChoice6()
{
    double a = 0, b = 0, c = 0;
    defineAParabola(a,b,c);
    while (myFAbs(a) <EPSILON) {
        std::cout << "a must be non-zero. Please enter new coefficients";
        defineAParabola(a, b, c);
    }
    std::cout << std::endl;
    double x = 0, y = 0;
    char name[NAMESIZE]="g";
    defineAPoint(name,x, y);
    while (myFAbs(y) > EPSILON) {
        std::cout << "The point " << name << " is not on the real line. Please enter a new point.";
        defineAPoint(name, x, y);
    }
    if (a * (a * x * x + b * x + c) < 0) {
        std::cout << "The point " << name << " is between the solutions of the equation a*x*x + b*x + c = 0 (x1 and x2). There is no tangent.\n";
    }
    else {
        double slope1 = 0, yIntercept1 = 0, thirdCoeff1 = 0, slope2 = 0, yIntercept2 = 0, thirdCoeff2 = 0;
        unsigned int linesCounter = 0;
        tangentToParabolaFromPoint(a, b, c, x, y, slope1, yIntercept1, thirdCoeff1, slope2, yIntercept2, thirdCoeff2, linesCounter);
        std::cout << "Equation of the tangent to the parabola from point " << name << std::endl;
        printALine(slope1, yIntercept1, thirdCoeff1);
        if (linesCounter == 2) {
            printALine(slope2, yIntercept2, thirdCoeff2);
        }
    }
}

void printMessageChoice7()
{
    double a = 0, b = 0, c = 0;
    defineAParabola(a, b, c);
    double lineSlCoeff = 0, lineYInt = 0, lineThirdCoeff = 0;
    char name1[NAMESIZE] = "g";
    defineALine(name1,lineSlCoeff, lineYInt, lineThirdCoeff);
    double point1X = 0, point1Y = 0, point2X = 0, point2Y = 0;
    unsigned int intersectionPointsCounter = 0;
    intersectionPointsLineAndParabola(a, b, c, lineSlCoeff, lineYInt, lineThirdCoeff, point1X, point1Y, point2X, point2Y, intersectionPointsCounter);
    if (intersectionPointsCounter == 0) {
        std::cout << "No intersection points." << std::endl;
    }
    else if (intersectionPointsCounter == 1) {
        std::cout << "Intersection point :" << std::endl;
        printApoint(point1X, point1Y);
    }
    else {
        std::cout << "Intersection point 1:" << std::endl;
        printApoint(point1X, point1Y);
        std::cout << "Intersection point 2:" << std::endl;
        printApoint(point2X, point2Y);
    }
}

void inIt(bool* arr, size_t length, bool value)
{
    for (size_t i = 0; i < length; i++) {
        arr[i] = value;
    }
}

void printMessageChoice8()
{
    double linesCoefficients[12];
    //index 0-line1slCoeff, 1-line1yIntercept, 2-line3thirdCoeff, from 3 to 5- line2 Coefficients, from 6 to 8 line3 coefficients and from 9 to 11- line4 coeff
    for (int i = 0; i <=9; i+=3) {
        std::cout<<"Enter coefficients to define a line (Ax + By + C = 0)"<<std::endl;
        std::cin >> linesCoefficients[i] >> linesCoefficients[i + 1] >> linesCoefficients[i + 2];
    }
    double intersectionPoints[12];
    //index 0 and 1 inters. point coordinates(0-x, 1-y)lines 1 and 2(if it exist), index 2 and 3(inter. point lines 1 and 3), index 4 and 5(inter. point lines 1 and 4)
    //index 6 and 7(inters. point lines 2 and 3), index 8 and 9 (inters. point lines 2 and 4), index 10 and 11(inters. point lines 3 and 4)
    bool intersections[6];
    // index 0 -lines 1 and 2 have intr. point?, index 1 - lines 1 and 3 have intr. point?, index 2 - lines 1 and 4..., index 3 - lines 2 and 3...,  index 4 - lines 2 and 4..
    //index 5 - lines 3 and 4
    inIt(intersections, 6, false);
    char typeOfQuadrilateral[20];
    unsigned int line1IntersectionPoints = 0, line2IntersectionPoints = 0, line3IntersectionPoints = 0, line4IntersectionPoints = 0;
    if (isQuadrilateral(linesCoefficients, intersectionPoints, intersections, line1IntersectionPoints, line2IntersectionPoints, line3IntersectionPoints, line4IntersectionPoints)) {
        std::cout << findTypeOfQuadrilateral(linesCoefficients, intersectionPoints, intersections, line1IntersectionPoints, line2IntersectionPoints, line3IntersectionPoints, line4IntersectionPoints) << std::endl;
    }
    else {
        std::cout << "Not quadrilateral" << std::endl;
    }
}

void clearConsole() {
    std::cout << "\033[;H"; // Moves cursor to the top left
    std::cout << "\033[J"; // Clears the console
}
void printChoices( unsigned int choice)
{
    while (choice < 1 || choice>9) {
        checkInput(choice);
        switch (choice)
        {
        case 1: printMessageChoice1();break;
        case 2: printMessageChoice2();break;
        case 3: printMessageChoice3();break;
        case 4: printMessageChoice4();break;
        case 5: printMessageChoice5();break;
        case 6: printMessageChoice6();break;
        case 7: printMessageChoice7();break;
        case 8:printMessageChoice8();break;
        default:
            std::cout << "Incorrect input! Please try again and enter a number between 1 and 9.\n";
            break;
        }
    }
}
int main()
{
    unsigned int choice=11;
    char answer[12] = "yes";
    char yes[4] = "yes";
    char no[3] = "no";
    while (areStringsEqual(answer, yes))
    {
        std::cout << "Hello! What would be useful to you? \n"<<
            "1) Check whether a point is on a line;\n" <<
            "2) From a line g and point p, derive an equation of a line that is parallel to g and passes through p; \n" <<
            "3) Given a line g and a point p lying on it, derive an equation of a line perpendicular to g with a fifth at p \n" <<
            "4) For two lines, find their intersection if it exists; \n" <<
            "5) By triangle (set by three points) constructs equations of:\n    *)heights;\n    *)medians;\n    *)bisectors.\n" <<
            "6) From a given equation of a parabola(ax2 + bx + c = 0) and a point (on the real line) derive an equation of the tangent to the parabola from the given point;\n" <<
            "7) From a given equation of a parabola and a line, find their points of intersection;\n" <<
            "8) Give four equations of lines, determine the type of quadrilateral they form when they intersect.\n" <<
            "Please enter the number of your choice!" << std::endl;
        
        printChoices(choice);
        std::cout << "Do you want to run the program again? yes/no"<<std::endl;
        std::cin >> answer;
        choice = 10;
        while (!(areStringsEqual(yes, answer) || areStringsEqual(no, answer)) ){
            std::cout << "Invalid answer. Please enter your answer again. yes/no";
            std::cin >> answer;
        }
        clearConsole();
    }
}

