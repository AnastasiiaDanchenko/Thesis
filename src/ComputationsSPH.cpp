#include "..\headers\ComputationsSPH.h"

// Loop through all particles and compute density
void ComputeDensityPressure() {
    double averageDensity = 0.0f;
    int numFluidParticles = 0;

    for (auto& p : particles) {

        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Compute density
        double density = 0.0f;
        for (auto neighbor : p.neighbors) {
            const Eigen::Vector3d r = p.position - neighbor->position;
            density += neighbor->mass * CubicSplineKernel(r);
        }
        p.density = density;
        averageDensity += density;
        numFluidParticles++;

        // Compute pressure
        double pressure = parameters.stiffness * (p.density / parameters.restDensity - 1);
        p.pressure = std::max(0.0, pressure);
    }

    std::cout << "Average density: " << averageDensity / numFluidParticles << std::endl;
}

void ComputeAcceleration() {
    for (auto& p : particles) {
        if (p.isFluid == false) { // Skip boundary particles
			continue;
		}

        // Gravity force
        Eigen::Vector3d acceleration = parameters.gravity;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector3d r = p.position - neighbor->position;
            const double rSquaredNorm = r.squaredNorm();
            const Eigen::Vector3d kernel = CubicSplineKernelGradient(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / 
            pow(neighbor->density, 2)) * kernel;

            // Viscosity force
            const Eigen::Vector3d v = p.velocity - neighbor->velocity;
            acceleration += 2 * parameters.viscosity * neighbor->mass * v / neighbor->density * r.dot(kernel) / 
            (rSquaredNorm + 0.01f * pow(parameters.spacing, 2));
        }

        p.acceleration = acceleration;
    }
}

// Euler integration
void UpdateParticles() {
    double maxVelocity = 0.0;

#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];

        if (p.isFluid == true) {
            p.velocity = p.predictedVelocity + parameters.timeStep * p.pressureAcceleration;
        }
        p.position += parameters.timeStep * p.velocity;
        maxVelocity = std::max(maxVelocity, p.velocity.norm());

        // Delete particles that are out of bounds
        if (p.position.x() < - parameters.windowSize.width || p.position.x() > parameters.windowSize.width  * 2 || 
            p.position.y() < -parameters.windowSize.height || p.position.y() > parameters.windowSize.height * 2) {
			particles.erase(particles.begin() + i);
		}
    }

    parameters.timeStep = (maxVelocity < std::numeric_limits<double>::epsilon()) ? parameters.maxTimeStep : 
        std::min(parameters.maxTimeStep, 0.4 * parameters.spacing / maxVelocity);
}

// IISPH density computation
void ComputeDensity() {
    double avgDensity = 0.0;

#pragma omp parallel for reduction(+:avgDensity)
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];

        if (p.isFluid == false) {
            continue;
        }

        double density = 0.0;
        for (auto neighbor : p.neighbors) {
            density += neighbor->mass * CubicSplineKernel(p.position - neighbor->position);
        }
        p.density = density;
        avgDensity += density;
    }

    parameters.avgDensity = avgDensity / (parameters.particlesPerDimension.x * parameters.particlesPerDimension.y * 
    parameters.particlesPerDimension.z);
}

Eigen::Vector3d ViscosityAcceleration(Particle* p) {
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();

    for (auto neighbor : p->neighbors) {
        if (neighbor == p) { continue; }

        const Eigen::Vector3d r = p->position - neighbor->position;
        const Eigen::Vector3d kernel = CubicSplineKernelGradient(r);

        const Eigen::Vector3d v = p->velocity - neighbor->velocity;
        acceleration += 2 * parameters.viscosity * neighbor->mass * v / neighbor->density * r.dot(kernel) /
        (r.squaredNorm() + 0.01f * pow(parameters.spacing, 2));
    }

    return acceleration;
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        p.acceleration = parameters.gravity + ViscosityAcceleration(&p);
        p.predictedVelocity = p.velocity + parameters.timeStep * p.acceleration;
    }
}

// IISPH current density error computation - source term
void ComputeDensityError() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        double sourceTerm = parameters.restDensity - p.density;
        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                sourceTerm -= parameters.timeStep * neighbor->mass * 
                (p.predictedVelocity - neighbor->predictedVelocity).
                dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
            else {
                sourceTerm -= parameters.timeStep * neighbor->mass * 
                p.predictedVelocity.dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
        }

        p.sourceTerm = sourceTerm;
    }
}

// IISPH diagonal value computation
void ComputeLaplacian() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        double diagonal = 0.0;
        Eigen::Vector3d coefficient = Eigen::Vector3d::Zero();

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                coefficient -= neighbor->mass / pow(parameters.restDensity, 2) * 
                CubicSplineKernelGradient(p.position - neighbor->position);
            }
            else {
                coefficient -= 2 * parameters.gamma * neighbor->mass / pow(parameters.restDensity, 2) * 
                CubicSplineKernelGradient(p.position - neighbor->position);
            }
        }

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                diagonal += neighbor->mass * (coefficient + p.mass / pow(parameters.restDensity, 2) *
                CubicSplineKernelGradient(neighbor->position - p.position)).
                dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
            else {
                diagonal += neighbor->mass * 
                coefficient.dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
        }

        diagonal *= pow(parameters.timeStep, 2);
        p.diagonal = diagonal;

        if (p.diagonal != 0) {
            p.pressure = std::max(0.0, parameters.omega * p.sourceTerm / p.diagonal);
        }
        else {
            p.pressure = 0.0;
        }
    }
}

// IISPH compression convergence
void CompressionConvergence() {
    double densityError = parameters.restDensity;
    parameters.nbIterations = 0;
    parameters.firstErr = 0.0;

    while ((parameters.densityErr > parameters.errThreshold * parameters.restDensity || parameters.nbIterations < 2) 
        && parameters.nbIterations < 100) {
        densityError = 0.0;

        // first loop: compute pressure acceleration
#pragma omp parallel for
        for (int i = 0; i < particles.size(); i++) {
            Particle& p = particles[i];
            if (p.isFluid == false) { continue; }

            Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    acceleration -= neighbor->mass * (p.pressure / pow(parameters.restDensity, 2) + 
                    neighbor->pressure / pow(parameters.restDensity, 2)) *
                    CubicSplineKernelGradient(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * parameters.gamma * neighbor->mass * p.pressure / 
                    pow(parameters.restDensity, 2) * CubicSplineKernelGradient(p.position - neighbor->position);
                }
            }

            p.pressureAcceleration = acceleration;
        }

#pragma omp parallel for reduction(+:densityError)
        for (int i = 0; i < particles.size(); i++) {
            Particle& p = particles[i];
            if (p.isFluid == false) { continue; }

            double divergenceVel = 0.0;

            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    divergenceVel += neighbor->mass * (p.pressureAcceleration - neighbor->pressureAcceleration).
                    dot(CubicSplineKernelGradient(p.position - neighbor->position));
                }
                else {
                    divergenceVel += neighbor->mass * 
                    p.pressureAcceleration.dot(CubicSplineKernelGradient(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(parameters.timeStep, 2);

            if (abs(p.diagonal) >= 1e-6) {
                p.pressure = std::max(0.0, 
                p.pressure + parameters.omega * (p.sourceTerm - divergenceVel) / p.diagonal);
            }
            else {
                p.pressure = 0.0;
            }

            densityError += divergenceVel - p.sourceTerm;
        }

        densityError /= (parameters.particlesPerDimension.x * parameters.particlesPerDimension.y * 
            parameters.particlesPerDimension.z);
        parameters.densityErr = densityError;

        if (parameters.nbIterations == 0) {
            parameters.firstErr = parameters.densityErr;
        }

        parameters.nbIterations++;
    }
}

void ComputeDensityPressure2D() {
    double averageDensity = 0.0f;

    for (auto& p : particles2D) {

        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Compute density
        double density = 0.0f;
        for (auto neighbor : p.neighbors) {
            density += neighbor->mass * CubicSplineKernel2D(p.position - neighbor->position);
        }
        p.density = density;
        averageDensity += density;

        // Compute pressure
        double pressure = parameters.stiffness * (p.density / parameters.restDensity - 1);
        p.pressure = std::max(0.0, pressure);
    }

    parameters.avgDensity = averageDensity / (parameters.particlesPerDimension.x * parameters.particlesPerDimension.y);
}

void ComputeAcceleration2D() {
    for (auto& p : particles2D) {
        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Gravity force
        Eigen::Vector2d acceleration = parameters.gravity2D;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector2d r = p.position - neighbor->position;
            const double rSquaredNorm = r.squaredNorm();
            const Eigen::Vector2d kernel = CubicSplineKernelGradient2D(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure /
            pow(neighbor->density, 2)) * kernel;

            // Viscosity force
            const Eigen::Vector2d v = p.velocity - neighbor->velocity;
            acceleration += 2 * parameters.viscosity * neighbor->mass * v / neighbor->density * r.dot(kernel) /
            (rSquaredNorm + 0.01f * pow(parameters.spacing, 2));
        }

        p.acceleration = acceleration;
    }
}

void Update2D() {
    for (auto& p : particles2D) {
		if (p.isFluid == false) { continue; }
		p.velocity += parameters.timeStep * p.acceleration;
		p.position += parameters.timeStep * p.velocity;
	}
}

// Euler integration
void UpdateParticles2D() {
    double maxVelocity = 0.0;

#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];

        if (p.isFluid == true) {
            p.velocity = p.predictedVelocity + parameters.timeStep * p.pressureAcceleration;
        }
        p.position += parameters.timeStep * p.velocity;
        maxVelocity = std::max(maxVelocity, p.velocity.norm());
    }

    parameters.timeStep = (maxVelocity < std::numeric_limits<double>::epsilon()) ? parameters.maxTimeStep : 
        std::min(parameters.maxTimeStep, 0.4 * parameters.spacing / maxVelocity);
}

// IISPH density computation 2D
void ComputeDensity2D() {
    double avgDensity = 0.0;

#pragma omp parallel for reduction(+:avgDensity)
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];

        if (p.isFluid == false) { 
            continue; 
        }

        double density = 0.0;
        for (auto neighbor : p.neighbors) {
            density += neighbor->mass * CubicSplineKernel2D(p.position - neighbor->position);
        }
        p.density = density;
        avgDensity += density;
    }

    parameters.avgDensity = avgDensity / (parameters.particlesPerDimension.x * parameters.particlesPerDimension.y);

    /*if (0.4 * parameters.spacing / maxVelocity <= parameters.maxTimeStep) {
        parameters.timeStep = 0.4 * parameters.spacing / maxVelocity;
    }*/
}

void ComputeSurface2D() {
    for (auto& p : particles2D) {
		
        Eigen::Vector2d normal = Eigen::Vector2d::Zero();
        for (auto neighbor : p.neighbors) {
            normal += neighbor->mass / neighbor->density * CubicSplineKernelGradient2D(p.position - neighbor->position);
		}

		p.normal = parameters.spacing * normal;

        // Define surface particles
        if (p.normal.squaredNorm() > 0.05) {
			p.isSurface = true;
		}
        else {
			p.isSurface = false;
		}
	}
}

Eigen::Vector2d ViscosityAcceleration2D(Particle2D* p) {
    Eigen::Vector2d acceleration = Eigen::Vector2d::Zero();

    for (auto neighbor : p->neighbors) {
        if (neighbor == p) { continue; }

        const Eigen::Vector2d r = p->position - neighbor->position;
        const Eigen::Vector2d kernel = CubicSplineKernelGradient2D(r);

        const Eigen::Vector2d v = p->velocity - neighbor->velocity;
        acceleration += 2 * parameters.viscosity * neighbor->mass * v / neighbor->density * r.dot(kernel) /
            (r.squaredNorm() + 0.01f * pow(parameters.spacing, 2));
    }

    return acceleration;
}

Eigen::Vector2d SurfaceTensionAcceleration2D(Particle2D* p) {
	Eigen::Vector2d surfaceTension = Eigen::Vector2d::Zero();

    for (auto neighbor : p->neighbors) {
		if (neighbor == p || (p->isSurface == false && neighbor->isSurface == false)) { continue; }

        const Eigen::Vector2d r = p->position - neighbor->position;
        const Eigen::Vector2d n = p->normal - neighbor->normal;

        Eigen::Vector2d cohesion = Eigen::Vector2d::Zero();
        Eigen::Vector2d curvature = Eigen::Vector2d::Zero();

        cohesion = -parameters.cohesion * p->mass * neighbor->mass * CohesionSpline2D(r) / r.norm() * r;
        curvature = -parameters.cohesion * p->mass * n;

        surfaceTension += 2 * parameters.restDensity / (p->density + neighbor->density) * (cohesion + curvature);
	}

	return surfaceTension;
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) { continue; }

        if (parameters.surfaceTension) {
			p.acceleration = parameters.gravity2D + ViscosityAcceleration2D(&p) + SurfaceTensionAcceleration2D(&p);
		}
        else {
			p.acceleration = parameters.gravity2D + ViscosityAcceleration2D(&p);
		}

        p.predictedVelocity = p.velocity + parameters.timeStep * p.acceleration;
    }
}

// IISPH current density error computation - source term
void ComputeDensityError2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) { continue; }

        double sourceTerm = parameters.restDensity - p.density;
        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                sourceTerm -= parameters.timeStep * neighbor->mass * (p.predictedVelocity - 
                    neighbor->predictedVelocity).dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                sourceTerm -= parameters.timeStep * neighbor->mass * 
                    p.predictedVelocity.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
        }

        p.sourceTerm = sourceTerm;
    }
}

// IISPH diagonal value computation
void ComputeLaplacian2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) { continue; }

        double diagonal = 0.0;
        Eigen::Vector2d coefficient = Eigen::Vector2d::Zero();

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                coefficient -= neighbor->mass / pow(parameters.restDensity, 2) * 
                    CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
            else {
                coefficient -= 2 * parameters.gamma * neighbor->mass / pow(parameters.restDensity, 2) * 
                    CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
        }

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                diagonal += neighbor->mass * (coefficient + p.mass / pow(parameters.restDensity, 2) *
                    CubicSplineKernelGradient2D(neighbor->position - p.position)).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                diagonal += neighbor->mass * 
                    coefficient.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
        }

        diagonal *= pow(parameters.timeStep, 2);
        p.diagonal = diagonal;

        if (p.diagonal != 0) {
            p.pressure = std::max(0.0, parameters.omega * p.sourceTerm / p.diagonal);
        }
        else {
            p.pressure = 0.0;
        }
    }
}

// IISPH compression convergence
void CompressionConvergence2D() {
    double densityError = parameters.restDensity;
    parameters.nbIterations = 0;
    parameters.firstErr = 0.0;

    while ((parameters.densityErr > parameters.errThreshold * parameters.restDensity || parameters.nbIterations < 2) 
        && parameters.nbIterations < 1000) {
        double densityErrorPrev = densityError;
        densityError = 0.0;

        // first loop: compute pressure acceleration
#pragma omp parallel for
        for (int i = 0; i < particles2D.size(); i++) {
            Particle2D& p = particles2D[i];
            if (p.isFluid == false) { continue; }

            Eigen::Vector2d acceleration = Eigen::Vector2d::Zero();
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    acceleration -= neighbor->mass * (p.pressure / pow(parameters.restDensity, 2) + 
                    neighbor->pressure / pow(parameters.restDensity, 2)) *
                    CubicSplineKernelGradient2D(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * parameters.gamma * neighbor->mass * p.pressure / 
                    pow(parameters.restDensity, 2) * CubicSplineKernelGradient2D(p.position - neighbor->position);
                }
            }

            p.pressureAcceleration = acceleration;
        }

#pragma omp parallel for reduction(+:densityError)
        for (int i = 0; i < particles2D.size(); i++) {
            Particle2D& p = particles2D[i];
            if (p.isFluid == false) { continue; }

            double divergenceVel = 0.0;

            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    divergenceVel += neighbor->mass * (p.pressureAcceleration - neighbor->pressureAcceleration).
                        dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
                }
                else {
                    divergenceVel += neighbor->mass * 
                        p.pressureAcceleration.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(parameters.timeStep, 2);

            if (abs(p.diagonal) >= 1e-6) {
                p.pressure = std::max(0.0, p.pressure + parameters.omega * (p.sourceTerm - divergenceVel) / 
                    p.diagonal);
            }
            else {
				p.pressure = 0.0;
			}
           
            densityError += divergenceVel - p.sourceTerm;
        }

        densityError /= (parameters.particlesPerDimension.x * parameters.particlesPerDimension.y);
        parameters.densityErr = densityError;

        if (parameters.nbIterations == 0) {
            parameters.firstErr = parameters.densityErr;
        }

        parameters.nbIterations++;
    }
}

void BoundaryMassUpdate2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) {
			double kernelSum = 0.0;
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid == false) {
					kernelSum += CubicSplineKernel2D(p.position - neighbor->position);
				}
			}
			p.mass = parameters.restDensity * 0.75 / kernelSum;
		}
	}
}

void BoundaryMassUpdate() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) {
            double kernelSum = 0.0;
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid == false) {
                    kernelSum += CubicSplineKernel(p.position - neighbor->position);
                }
            }
            p.mass = parameters.restDensity * 0.8 / kernelSum;
        }
    }
}

void RotateBoundary2D() {
    double maxVelocity = 0.0;

    int width_max = parameters.windowSize.width / parameters.spacing / 5;
    int width_avg = parameters.windowSize.width / parameters.spacing / 10;
    int height_avg = parameters.windowSize.height / parameters.spacing / 7;

#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];

        if (p.isFluid == true) {
            p.velocity = p.predictedVelocity + parameters.timeStep * p.pressureAcceleration;
            p.position += parameters.timeStep * p.velocity;
        }
        else if (p.ID >= parameters.boundaryTestID) {
            const Eigen::Vector2d center = Eigen::Vector2d((width_max + width_avg) * parameters.spacing, height_avg * 
                parameters.spacing);

            double radius = (p.position - center).norm();
            double angularVelocity = 0.2;
            
            double angle = std::atan2(p.position.y() - center.y(), p.position.x() - center.x()) + angularVelocity * 
                parameters.timeStep;

            p.position = center + Eigen::Vector2d(radius * std::cos(angle), radius * std::sin(angle));
        }

        maxVelocity = std::max(maxVelocity, p.velocity.norm());
    }

    parameters.timeStep = (maxVelocity < std::numeric_limits<double>::epsilon()) ? parameters.maxTimeStep : 
        std::min(parameters.maxTimeStep, 0.4 * parameters.spacing / maxVelocity);
}

void UpdateGhosts() {
#pragma omp parallel for
    for (int i = 0; i < ghostParticles.size(); i++) {
		Particle& g = ghostParticles[i];
		
        double density = 0, pressure = 0, weight = 0;
        Eigen::Vector3d velocity = Eigen::Vector3d::Zero();

        for (auto p : g.neighbors) {
            // weighted sum of the density of the neighbors depending on the distance
            double distance = (g.position - p->position).norm();
            weight += 1 / distance;

            density += 1 / distance * p->density;
            pressure += 1 / distance * p->pressure;
            velocity += 1 / distance * p->velocity;
		}
        
		g.density = density / weight;
        g.pressure = pressure / weight;
		g.velocity.x() = velocity.x() / weight;
        g.velocity.y() = velocity.y() / weight;

        g.position += parameters.timeStep * g.velocity;
	}
}
