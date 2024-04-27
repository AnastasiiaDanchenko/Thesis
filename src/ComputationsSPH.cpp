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
		p.velocity += TIME_STEP * (p.acceleration + p.pressureAcceleration);
		// Update position
		p.position += TIME_STEP * p.velocity;
    }
}

// IISPH density computation
void ComputeDensity() {
    for (auto& p : particles) {
		if (p.isFluid == false) { continue; }

		float density = 0.0f;
        for (auto neighbor : p.neighbors) {
			const Eigen::Vector3f r = p.position - neighbor->position;
			density += neighbor->mass * CubicSplineKernel(r);
		}
		p.density = density;
	}
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity() {
    for (auto& p : particles) {
        // Compute non-pressure acceleration
		if (p.isFluid == false) { continue; }

        // Acceleration due to gravity
		Eigen::Vector3f acceleration = GRAVITY;
        for (auto neighbor : p.neighbors) {
			if (neighbor == &p) { continue; }
			const Eigen::Vector3f r = p.position - neighbor->position;
			const float rSquaredNorm = r.squaredNorm();
			const Eigen::Vector3f kernel = CubicSplineKernelGradient(r);

            // Viscosity acceleration
            const Eigen::Vector3f v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
		}

		p.acceleration = acceleration; 
        p.predictedVelocity = p.velocity + TIME_STEP * p.acceleration;
	}
}

// IISPH density error computation - source term
void ComputeDensityError() {
    for (auto& p : particles) {
		if (p.isFluid == false) { continue; }

        float sourceTerm = REST_DENSITY - p.density;
        for (auto neighbor : p.neighbors) { 
            if (neighbor->isFluid) {
                sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->predictedVelocity).
                              dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
            else {
				sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->velocity).
                              dot(CubicSplineKernelGradient(p.position - neighbor->position));
			}
        }

        p.sourceTerm = sourceTerm;
	}
}

// IISPH diagonal value computation
void ComputeLaplacian() {
    for (auto& p : particles) {
		if (p.isFluid == false) { continue; }

        Eigen::Vector3f coefficient = Eigen::Vector3f::Zero();
        float diagonal = 0.0f;

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                coefficient -= neighbor->mass / pow(REST_DENSITY, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
			} 
            else {
                coefficient -= 2 * GAMMA * neighbor->mass / pow(REST_DENSITY, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
            }
		}

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
				diagonal += neighbor->mass * (coefficient + p.mass / pow(REST_DENSITY, 2) * 
                            CubicSplineKernelGradient(neighbor->position - p.position)).
                            dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
            else {
                diagonal += neighbor->mass * coefficient.dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
        }

        diagonal *= pow(TIME_STEP, 2);
        p.diagonal = diagonal;

        // initialize pressure value with 0
        p.pressure = 0.0f;
	}
}

// IISPH compression convergence
void CompressionConvergence() {
    float densityError = REST_DENSITY;

    for (int l = 0; l < 3; l++) {
        densityError = 0.0f;

        // first loop: compute pressure acceleration
        for (auto& p : particles) {
            if (p.isFluid == false) { continue; }

            Eigen::Vector3f acceleration = Eigen::Vector3f::Zero();
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2)) *
                        CubicSplineKernelGradient(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * GAMMA * neighbor->mass * p.pressure / pow(p.density, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
                }
            }

            p.pressureAcceleration = acceleration;
        }

        // second loop: update pressure
        for (auto& p : particles) {
            if (p.isFluid == false) { continue; }

            float divergenceVel = 0.0f;
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    divergenceVel += neighbor->mass * (p.pressureAcceleration - neighbor->pressureAcceleration).
                        dot(CubicSplineKernelGradient(p.position - neighbor->position));
                }
                else {
                    divergenceVel += neighbor->mass * p.pressureAcceleration.dot(CubicSplineKernelGradient(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(TIME_STEP, 2);

            if (p.neighbors.size() == 0) {
                p.diagonal = 0;
            }
            else if (p.diagonal != 0) {
                p.pressure = std::max(0.0f, p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal);
            }

            densityError += divergenceVel - p.sourceTerm;
        }

        // Compute average density error
        densityError /= particles.size();
    }
}
