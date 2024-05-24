#include "..\headers\Kernel.h"

// Compute cubic spline kernel function
float CubicSplineKernel(const Eigen::Vector3f& r) {
    float alpha = 1 / (M_PI * pow(SPACING, 3));
    float q = r.norm() / SPACING;
    float result = 0.0f;

    if (q >= 0 && q < 1) {
		result = alpha * (1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3));
	} else if (q >= 1 && q < 2) {
		result = alpha * 0.25 * pow(2 - q, 3);
	}
    return result;
}

// Compute cubic spline kernel gradient
Eigen::Vector3f CubicSplineKernelGradient(const Eigen::Vector3f& r) {
    float alpha = 1 / (M_PI * pow(SPACING, 3));
    float q = r.norm() / SPACING;
    float derivative = 0.0f;
    
    if (q >= 0 && q < 1) {
        derivative = alpha / SPACING * (-3 * q + 2.25 * pow(q, 2));
    } else if (q >= 1 && q < 2) {
        derivative = alpha / SPACING * (-0.75 * pow(2 - q, 2));
	}
    return derivative * r.normalized();
}

float CubicSplineKernel2D(Eigen::Vector2f r) {
    float alpha = 10 / (7 * M_PI * pow(SPACING, 2));
    float q = r.norm() / SPACING;

    if (q >= 0 && q < 1) {
        return alpha * (1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3));
    }
    else if (q >= 1 && q < 2) {
        return alpha * 0.25 * pow(2 - q, 3);
    }
    return 0.0f;
}

// Compute cubic spline kernel gradient
Eigen::Vector2f CubicSplineKernelGradient2D(Eigen::Vector2f r) {
    float alpha = 10 / (7 * M_PI * pow(SPACING, 2));
    float q = r.norm() / SPACING;
    float derivative = 0.0f;

    if (q >= 0 && q < 1) {
        derivative = alpha / SPACING * (-3 * q + 2.25 * pow(q, 2));
    }
    else if (q >= 1 && q < 2) {
        derivative = alpha / SPACING * (-0.75 * pow(2 - q, 2));
    }
    return derivative * r.normalized();
}

float CohesionSpline2D(Eigen::Vector2f r) {
    float alpha = 32 / (M_PI * pow(SPACING, 9));
    float q = r.norm();
    float result = 0.0f;

    if (2 * q > SPACING && q <= SPACING) {
        result = alpha * pow(SPACING - q, 3) * pow(q, 3);
    }
    else if (q > 0 && 2 * q <= SPACING) {
        result = alpha * (2 * pow(SPACING - q, 3) * pow(q, 3) - pow(SPACING, 6) / 64);
	}

    return result;
}

