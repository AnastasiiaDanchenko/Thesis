#include "..\headers\ComputationsSPH.h"

// Loop through all particles and compute density
void ComputeDensityPressure() {
    for (auto& p : particles) {

        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Compute density
        float density = 0.0f;
        for (auto neighbor : p.neighbors) {
            const Eigen::Vector3f r = p.position - neighbor->position;
            density += neighbor->mass * CubicSplineKernel(r);
        }
        p.density = density;

        // Compute pressure
        float pressure = STIFFNESS * (p.density / REST_DENSITY - 1);
        p.pressure = std::max(0.0f, pressure);
    }
}

void ComputeAcceleration() {
    for (auto& p : particles) {
        if (p.isFluid == false) { // Skip boundary particles
			continue;
		}

        // Gravity force
        Eigen::Vector3f acceleration = GRAVITY;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector3f r = p.position - neighbor->position;
            const float rSquaredNorm = r.squaredNorm();
            const Eigen::Vector3f kernel = CubicSplineKernelGradient(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2)) 
                * kernel;

            // Viscosity force
            const Eigen::Vector3f v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) / 
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
        }

        p.acceleration = acceleration;
    }
}

// Euler integration
void UpdateParticles() {
    for (auto& p : particles) {
        if (p.isFluid == false) { continue; } // Skip boundary particles
		// Update velocity
		p.velocity += TIME_STEP * p.acceleration;
		// Update position
		p.position += TIME_STEP * p.velocity;
    }
}
