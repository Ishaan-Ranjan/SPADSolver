#pragma once

#include <iostream>
#include <cmath>
#include <vector>

double alpha(double AlFrac, double Field, double Temp);

double beta(double AlFrac, double Field, double Temp);

double gammaLowEnergy(double AlFrac, double Energy);

double gammaHighEnergy(double AlFrac, double Energy);

double gamma(double AlFrac, double Energy);