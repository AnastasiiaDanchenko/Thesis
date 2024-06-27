#include "..\headers\Kernel.h"

// Compute cubic spline kernel function
double CubicSplineKernel(const Eigen::Vector3d& r) {
    double alpha = 1 / (M_PI * pow(SPACING, 3));
    double q = r.norm() / SPACING;
    double result = 0.0f;

    if (q >= 0 && q < 1) {
		result = alpha * (1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3));
	} else if (q >= 1 && q < 2) {
		result = alpha * 0.25 * pow(2 - q, 3);
	}
    return result;
}

// Compute cubic spline kernel gradient
Eigen::Vector3d CubicSplineKernelGradient(const Eigen::Vector3d& r) {
    double alpha = 1 / (M_PI * pow(SPACING, 3));
    double q = r.norm() / SPACING;
    double derivative = 0.0f;
    
    if (q >= 0 && q < 1) {
        derivative = alpha / SPACING * (-3 * q + 2.25 * pow(q, 2));
    } else if (q >= 1 && q < 2) {
        derivative = alpha / SPACING * (-0.75 * pow(2 - q, 2));
	}
    return derivative * r.normalized();
}

double CubicSplineKernel2D(Eigen::Vector2d r) {
    double alpha = 5 / (14 * M_PI * pow(SPACING, 2));
    double q = r.norm() / SPACING;
    double result = 0.0;

    if (q >= 0 && q < 1) {
        result = alpha * (pow(2 - q, 3) - 4 * pow(1 - q, 3));
    }
    else if (q >= 1 && q < 2) {
        result = alpha * pow(2 - q, 3);
    }
    return result;
}

// Compute cubic spline kernel gradient
Eigen::Vector2d CubicSplineKernelGradient2D(Eigen::Vector2d r) {
    double alpha = 10 / (7 * M_PI * pow(SPACING, 2));
    double q = r.norm() / SPACING;
    double derivative = 0.0;

    if (q >= 0 && q < 1) {
        derivative = alpha / SPACING * (-3 * q + 2.25 * pow(q, 2));
    }
    else if (q >= 1 && q < 2) {
        derivative = alpha / SPACING * (-0.75 * pow(2 - q, 2));
    }
    return derivative * r.normalized();
}

double CohesionSpline2D(Eigen::Vector2d r) {
    double alpha = 32 / (M_PI * pow(SPACING, 9));
    double q = r.norm();
    double result = 0.0f;

    if (2 * q > SPACING && q <= SPACING) {
        result = alpha * pow(SPACING - q, 3) * pow(q, 3);
    }
    else if (q > 0 && 2 * q <= SPACING) {
        result = alpha * (2 * pow(SPACING - q, 3) * pow(q, 3) - pow(SPACING, 6) / 64);
	}

    if (FIRST_STEP) {
        result *= FIRST_STEP_CORRECTION;
    }

    return result;
}

