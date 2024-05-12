#include "..\headers\ComputationsSPH.h"

// Loop through all particles and compute density
void ComputeDensityPressure() {
    float averageDensity = 0.0f;
    int numFluidParticles = 0;

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
        averageDensity += density;
        numFluidParticles++;

        // Compute pressure
        float pressure = STIFFNESS * (p.density / REST_DENSITY - 1);
        p.pressure = std::max(0.0f, pressure);
    }

    std::cout << "Average density: " << averageDensity / numFluidParticles << std::endl;
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
    float avgDensity = 0.0f;
    int numFluidParticles = 0;

    for (auto& p : particles) {
		if (p.isFluid == false) { continue; }

		float density = 0.0f;
        for (auto neighbor : p.neighbors) {
			density += neighbor->mass * CubicSplineKernel(p.position - neighbor->position);
		}
		p.density = density;
        avgDensity += density;
        numFluidParticles++;
	}

    avgDensity /= numFluidParticles;
    //std::cout << "Average density: " << avgDensity << std::endl;
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity() {
    for (auto& p : particles) {
        // Compute non-pressure acceleration
		if (p.isFluid == false) { continue; }

        // Acceleration due to gravity
		Eigen::Vector3f acceleration = GRAVITY;

        // Viscosity acceleration
        for (auto neighbor : p.neighbors) {
			if (neighbor == &p) { continue; }
			const Eigen::Vector3f r = p.position - neighbor->position;
			const float rSquaredNorm = r.squaredNorm();
			const Eigen::Vector3f kernel = CubicSplineKernelGradient(r);

            const Eigen::Vector3f v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
		}

		p.acceleration = acceleration; 
        p.predictedVelocity = p.velocity + TIME_STEP * p.acceleration;
	}
}

// IISPH current density error computation - source term
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
                // Boundary velocity is 0
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

        float diagonal = 0.0f;
        Eigen::Vector3f coefficient = Eigen::Vector3f::Zero();

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                coefficient -= neighbor->mass / pow(p.density, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
			} 
            else {
                coefficient -= 2 * GAMMA * neighbor->mass / pow(p.density, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
            }
		}

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                diagonal += neighbor->mass * (coefficient + p.mass / pow(p.density, 2) *
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
    int l = 0;

    //while (densityError > 0.01 * REST_DENSITY || l < 3) {
    while (l < 10) {
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
                    acceleration -= 2 * GAMMA * neighbor->mass * p.pressure / pow(p.density, 2) * 
                                    CubicSplineKernelGradient(p.position - neighbor->position);
                }
            }

            p.pressureAcceleration = acceleration;
        }

        // second loop: update pressure
        for (auto& p : particles) {
            if (p.isFluid == false) { continue; }

            float divergenceVel = 0.0f;
            float pressure = p.pressure;

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
                pressure = std::max(0.0f, p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal);
            }

            p.pressure = pressure;
            densityError += abs(divergenceVel - p.sourceTerm);
        }

        densityError /= 1300;
        std::cout << "Average density error: " << densityError << std::endl;
        l++;
    }
    std::cout << std::endl;
}

void ComputeDensityPressure2D() {
    float averageDensity = 0.0f;
    int numFluidParticles = 0;

    for (auto& p : particles2D) {

        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Compute density
        float density = 0.0f;
        for (auto neighbor : p.neighbors) {
            const Eigen::Vector2f r = p.position - neighbor->position;
            density += neighbor->mass * CubicSplineKernel2D(r);
        }
        p.density = density;
        averageDensity += density;
        numFluidParticles++;

        // Compute pressure
        float pressure = STIFFNESS * (p.density / REST_DENSITY - 1);
        p.pressure = std::max(0.0f, pressure);
    }

    //std::cout << "Average density: " << averageDensity / numFluidParticles << std::endl;
}

void ComputeAcceleration2D() {
    for (auto& p : particles2D) {
        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Gravity force
        Eigen::Vector2f acceleration = GRAVITY2D;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector2f r = p.position - neighbor->position;
            const float rSquaredNorm = r.squaredNorm();
            const Eigen::Vector2f kernel = CubicSplineKernelGradient2D(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2))
                * kernel;

            // Viscosity force
            const Eigen::Vector2f v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
        }

        p.acceleration = acceleration;
    }
}

// Euler integration
void UpdateParticles2D() {
    //#pragma omp parallel for
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; } // Skip boundary particles
        // Update velocity
        p.velocity = p.predictedVelocity + TIME_STEP * p.pressureAcceleration;
        // Update position
        p.position += TIME_STEP * p.velocity;
    }
}

// IISPH density computation 2D
void ComputeDensity2D() {
    float avgDensity = 0.0f;

    //#pragma omp parallel for
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }

        float density = 0.0f;
        for (auto neighbor : p.neighbors) {
            density += neighbor->mass * CubicSplineKernel2D(p.position - neighbor->position);
        }
        p.density = density;
        avgDensity += density;
    }

    AVG_DENSITY = avgDensity / (PARTICLES_X * PARTICLES_Y);
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity2D() {
    //#pragma omp parallel for
    for (auto& p : particles2D) {
        // Compute non-pressure acceleration
        if (p.isFluid == false) { continue; }

        // Acceleration due to gravity
        Eigen::Vector2f acceleration = GRAVITY2D;

        // Viscosity acceleration
        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; }
            const Eigen::Vector2f r = p.position - neighbor->position;
            const float rSquaredNorm = r.squaredNorm();
            const Eigen::Vector2f kernel = CubicSplineKernelGradient2D(r);

            const Eigen::Vector2f v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
        }

        p.acceleration = acceleration;
        p.predictedVelocity = p.velocity + TIME_STEP * p.acceleration;
    }
}

// IISPH current density error computation - source term
void ComputeDensityError2D() {
    //#pragma omp parallel for
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }

        float sourceTerm = REST_DENSITY - p.density;
        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->predictedVelocity).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                // Boundary velocity is 0
                sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->velocity).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
        }

        p.sourceTerm = sourceTerm;
    }
}

// IISPH diagonal value computation
void ComputeLaplacian2D() {
    //#pragma omp parallel for
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }

        float diagonal = 0.0f;
        Eigen::Vector2f coefficient = Eigen::Vector2f::Zero();

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                coefficient -= neighbor->mass / pow(p.density, 2) * CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
            else {
                coefficient -= 2 * GAMMA * neighbor->mass / pow(p.density, 2) * CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
        }

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                diagonal += neighbor->mass * (coefficient + p.mass / pow(p.density, 2) *
                    CubicSplineKernelGradient2D(neighbor->position - p.position)).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                diagonal += neighbor->mass * coefficient.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
        }

        diagonal *= pow(TIME_STEP, 2);
        p.diagonal = diagonal;

        // initialize pressure value with 0
        //p.pressure = 0.0f;
        if (p.diagonal != 0) {
            p.pressure = std::max(0.0f, OMEGA * p.sourceTerm / p.diagonal);
        }
        else {
			p.pressure = 0.0f;
        }
    }
}

// IISPH compression convergence
void CompressionConvergence2D() {
    float densityError = REST_DENSITY;
    int l = 0;

    while (densityError > ERR_THRESHOLD * REST_DENSITY || l < 2) {
    //while (l < 10) {
        densityError = 0.0f;

        // first loop: compute pressure acceleration
        //#pragma omp parallel for
        for (auto& p : particles2D) {
            if (p.isFluid == false) { continue; }

            Eigen::Vector2f acceleration = Eigen::Vector2f::Zero();
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2)) *
                        CubicSplineKernelGradient2D(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * GAMMA * neighbor->mass * p.pressure / pow(p.density, 2) *
                        CubicSplineKernelGradient2D(p.position - neighbor->position);
                }
            }

            p.pressureAcceleration = acceleration;
        }

        // second loop: update pressure        
        //#pragma omp parallel for
        for (auto& p : particles2D) {
            if (p.isFluid == false) { continue; }

            float divergenceVel = 0.0f;
            float pressure = p.pressure;

            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    divergenceVel += neighbor->mass * (p.pressureAcceleration - neighbor->pressureAcceleration).
                        dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
                }
                else {
                    divergenceVel += neighbor->mass * p.pressureAcceleration.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(TIME_STEP, 2);

            if (p.diagonal != 0) {
                pressure = std::max(0.0f, p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal);
                //std::cout << "Pressure: " << pressure << std::endl;
                //std::cout << "Second term: " << p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal << std::endl;
            }

            densityError += abs(divergenceVel - p.sourceTerm);

            p.pressure = pressure;
        }

        densityError /= PARTICLES_X * PARTICLES_Y;
        DENSITY_ERR = densityError / 1000 * 100;

        l++;
    }
}
