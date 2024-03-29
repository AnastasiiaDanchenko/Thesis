#include "..\headers\Kernel.h"

// Compute cubic spline kernel function
float CubicSplineKernel(const Eigen::Vector3f& r) {
    float alpha = 4 / (4 * M_PI * pow(SPACING, 3));
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
    float alpha = 4 / (4 * M_PI * pow(SPACING, 3));
    float q = r.norm() / SPACING;
    float derivative = 0.0f;
    
    if (q >= 0 && q < 1) {
        derivative = alpha / SPACING * (-3 * q + 2.25 * pow(q, 2));
    } else if (q >= 1 && q < 2) {
        derivative = alpha / SPACING * (-0.75 * pow(2 - q, 2));
	}
    return derivative * r.normalized();
}