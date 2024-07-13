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
        double pressure = STIFFNESS * (p.density / REST_DENSITY - 1);
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
        Eigen::Vector3d acceleration = GRAVITY;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector3d r = p.position - neighbor->position;
            const double rSquaredNorm = r.squaredNorm();
            const Eigen::Vector3d kernel = CubicSplineKernelGradient(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2)) 
                * kernel;

            // Viscosity force
            const Eigen::Vector3d v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) / 
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
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
            p.velocity = p.predictedVelocity + TIME_STEP * p.pressureAcceleration;
        }
        p.position += TIME_STEP * p.velocity;
        maxVelocity = std::max(maxVelocity, p.velocity.norm());
    }

    TIME_STEP = (maxVelocity < std::numeric_limits<double>::epsilon()) ? MAX_TIME_STEP : std::min(MAX_TIME_STEP, 0.4 * SPACING / maxVelocity);
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

    AVG_DENSITY = avgDensity / (PARTICLES_X * PARTICLES_Y * PARTICLES_Z);
}

Eigen::Vector3d ViscosityAcceleration(Particle* p) {
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();

    for (auto neighbor : p->neighbors) {
        if (neighbor == p) { continue; }

        const Eigen::Vector3d r = p->position - neighbor->position;
        const Eigen::Vector3d kernel = CubicSplineKernelGradient(r);

        const Eigen::Vector3d v = p->velocity - neighbor->velocity;
        acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
            (r.squaredNorm() + 0.01f * pow(SPACING, 2));
    }

    return acceleration;
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        p.acceleration = GRAVITY + ViscosityAcceleration(&p);
        p.predictedVelocity = p.velocity + TIME_STEP * p.acceleration;
    }
}

// IISPH current density error computation - source term
void ComputeDensityError() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        double sourceTerm = REST_DENSITY - p.density;
        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->predictedVelocity).
                    dot(CubicSplineKernelGradient(p.position - neighbor->position));
            }
            else {
                sourceTerm -= TIME_STEP * neighbor->mass * p.predictedVelocity.dot(CubicSplineKernelGradient(p.position - neighbor->position));
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

        if (p.diagonal != 0) {
            p.pressure = std::max(0.0, OMEGA * p.sourceTerm / p.diagonal);
        }
        else {
            p.pressure = 0.0;
        }
    }
}

// IISPH compression convergence
void CompressionConvergence() {
    double densityError = REST_DENSITY;
    NB_ITERATIONS = 0;
    FIRST_ERR = 0.0;

    while ((DENSITY_ERR > ERR_THRESHOLD * REST_DENSITY || NB_ITERATIONS < 2) && NB_ITERATIONS < 100) {
        densityError = 0.0;

        // first loop: compute pressure acceleration
#pragma omp parallel for
        for (int i = 0; i < particles.size(); i++) {
            Particle& p = particles[i];
            if (p.isFluid == false) { continue; }

            Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
            for (auto neighbor : p.neighbors) {
                if (neighbor->isFluid) {
                    acceleration -= neighbor->mass * (p.pressure / pow(REST_DENSITY, 2) + neighbor->pressure / pow(REST_DENSITY, 2)) *
                        CubicSplineKernelGradient(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * GAMMA * neighbor->mass * p.pressure / pow(REST_DENSITY, 2) *
                        CubicSplineKernelGradient(p.position - neighbor->position);
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
                    divergenceVel += neighbor->mass * p.pressureAcceleration.dot(CubicSplineKernelGradient(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(TIME_STEP, 2);

            if (abs(p.diagonal) >= 1e-6) {
                p.pressure = std::max(0.0, p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal);
            }
            else {
                p.pressure = 0.0;
            }

            densityError += divergenceVel - p.sourceTerm;
        }

        densityError /= (PARTICLES_X * PARTICLES_Y * PARTICLES_Z);
        DENSITY_ERR = densityError;

        if (NB_ITERATIONS == 0) {
            FIRST_ERR = DENSITY_ERR;
        }

        NB_ITERATIONS++;
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
        double pressure = STIFFNESS * (p.density / REST_DENSITY - 1);
        p.pressure = std::max(0.0, pressure);
    }

    AVG_DENSITY = averageDensity / (PARTICLES_X * PARTICLES_Y);
}

void ComputeAcceleration2D() {
    for (auto& p : particles2D) {
        if (p.isFluid == false) { // Skip boundary particles
            continue;
        }

        // Gravity force
        Eigen::Vector2d acceleration = GRAVITY2D;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; } // Skip self
            const Eigen::Vector2d r = p.position - neighbor->position;
            const double rSquaredNorm = r.squaredNorm();
            const Eigen::Vector2d kernel = CubicSplineKernelGradient2D(r);

            // Pressure force
            acceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2))
                * kernel;

            // Viscosity force
            const Eigen::Vector2d v = p.velocity - neighbor->velocity;
            acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
                (rSquaredNorm + 0.01f * pow(SPACING, 2));
        }

        p.acceleration = acceleration;
    }
}

void Update2D() {
    for (auto& p : particles2D) {
		if (p.isFluid == false) { continue; }
		p.velocity += TIME_STEP * p.acceleration;
		p.position += TIME_STEP * p.velocity;
	}
}

// Euler integration
void UpdateParticles2D() {
    double maxVelocity = 0.0;

#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];

        if (p.isFluid == true) {
            p.velocity = p.predictedVelocity + TIME_STEP * p.pressureAcceleration;
        }
        p.position += TIME_STEP * p.velocity;
        maxVelocity = std::max(maxVelocity, p.velocity.norm());
    }

    TIME_STEP = (maxVelocity < std::numeric_limits<double>::epsilon()) ? MAX_TIME_STEP : std::min(MAX_TIME_STEP, 0.4 * SPACING / maxVelocity);
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

    AVG_DENSITY = avgDensity / (PARTICLES_X * PARTICLES_Y);

    /*if (0.4 * SPACING / maxVelocity <= MAX_TIME_STEP) {
        TIME_STEP = 0.4 * SPACING / maxVelocity;
    }*/
}

void ComputeSurface2D() {
    for (auto& p : particles2D) {
		
        Eigen::Vector2d normal = Eigen::Vector2d::Zero();
        for (auto neighbor : p.neighbors) {
            normal += neighbor->mass / neighbor->density * CubicSplineKernelGradient2D(p.position - neighbor->position);
		}

		p.normal = SPACING * normal;

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
        acceleration += 2 * VISCOSITY * neighbor->mass * v / neighbor->density * r.dot(kernel) /
            (r.squaredNorm() + 0.01f * pow(SPACING, 2));
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

        cohesion = -COHESION * p->mass * neighbor->mass * CohesionSpline2D(r) / r.norm() * r;
        curvature = -COHESION * p->mass * n;

        surfaceTension += 2 * REST_DENSITY / (p->density + neighbor->density) * (cohesion + curvature);
	}

	return surfaceTension;
}

// IISPH non-pressure acceleration computation and velocity prediction
void PredictVelocity2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) { continue; }

        if (SURFACE_TENSION) {
			p.acceleration = GRAVITY2D + ViscosityAcceleration2D(&p) + SurfaceTensionAcceleration2D(&p);
		}
        else {
			p.acceleration = GRAVITY2D + ViscosityAcceleration2D(&p);
		}

        p.predictedVelocity = p.velocity + TIME_STEP * p.acceleration;
    }
}

// IISPH current density error computation - source term
void ComputeDensityError2D() {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];
        if (p.isFluid == false) { continue; }

        double sourceTerm = REST_DENSITY - p.density;
        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                sourceTerm -= TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->predictedVelocity).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                sourceTerm -= TIME_STEP * neighbor->mass * p.predictedVelocity.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
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
                coefficient -= neighbor->mass / pow(REST_DENSITY, 2) * CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
            else {
                coefficient -= 2 * GAMMA * neighbor->mass / pow(REST_DENSITY, 2) * CubicSplineKernelGradient2D(p.position - neighbor->position);
            }
        }

        for (auto neighbor : p.neighbors) {
            if (neighbor->isFluid) {
                diagonal += neighbor->mass * (coefficient + p.mass / pow(REST_DENSITY, 2) *
                    CubicSplineKernelGradient2D(neighbor->position - p.position)).
                    dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
            else {
                diagonal += neighbor->mass * coefficient.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
            }
        }

        diagonal *= pow(TIME_STEP, 2);
        p.diagonal = diagonal;

        if (p.diagonal != 0) {
            p.pressure = std::max(0.0, OMEGA * p.sourceTerm / p.diagonal);
        }
        else {
            p.pressure = 0.0;
        }
    }
}

// IISPH compression convergence
void CompressionConvergence2D() {
    double densityError = REST_DENSITY;
    NB_ITERATIONS = 0;
    FIRST_ERR = 0.0;

    while ((DENSITY_ERR > ERR_THRESHOLD * REST_DENSITY || NB_ITERATIONS < 2) && NB_ITERATIONS < 1000) {
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
                    acceleration -= neighbor->mass * (p.pressure / pow(REST_DENSITY, 2) + neighbor->pressure / pow(REST_DENSITY, 2)) *
                        CubicSplineKernelGradient2D(p.position - neighbor->position);
                }
                else {
                    acceleration -= 2 * GAMMA * neighbor->mass * p.pressure / pow(REST_DENSITY, 2) *
                        CubicSplineKernelGradient2D(p.position - neighbor->position);
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
                    divergenceVel += neighbor->mass * p.pressureAcceleration.dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
                }
            }

            divergenceVel *= pow(TIME_STEP, 2);

            if (abs(p.diagonal) >= 1e-6) {
                p.pressure = std::max(0.0, p.pressure + OMEGA * (p.sourceTerm - divergenceVel) / p.diagonal);
            }
            else {
				p.pressure = 0.0;
			}
           
            densityError += divergenceVel - p.sourceTerm;
        }

        densityError /= (PARTICLES_X * PARTICLES_Y);
        DENSITY_ERR = densityError;

        if (NB_ITERATIONS == 0) {
            FIRST_ERR = DENSITY_ERR;
        }

        NB_ITERATIONS++;
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
			p.mass = REST_DENSITY * 0.75 / kernelSum;
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
            p.mass = REST_DENSITY * 0.8 / kernelSum;
        }
    }
}

void RotateBoundary2D() {
    double maxVelocity = 0.0;

    int width_max = WINDOW_WIDTH / SPACING / 5;
    int width_avg = WINDOW_WIDTH / SPACING / 10;
    int height_avg = WINDOW_HEIGHT / SPACING / 7;

#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        Particle2D& p = particles2D[i];

        if (p.isFluid == true) {
            p.velocity = p.predictedVelocity + TIME_STEP * p.pressureAcceleration;
            p.position += TIME_STEP * p.velocity;
        }
        else if (p.ID >= ROTATING_BOUNDARY_ID) {
            const Eigen::Vector2d center = Eigen::Vector2d((width_max + width_avg) * SPACING, height_avg * SPACING);

            double radius = (p.position - center).norm();
            double angularVelocity = 0.2;
            
            double angle = std::atan2(p.position.y() - center.y(), p.position.x() - center.x()) + angularVelocity * TIME_STEP;

            p.position = center + Eigen::Vector2d(radius * std::cos(angle), radius * std::sin(angle));
        }

        maxVelocity = std::max(maxVelocity, p.velocity.norm());
    }

    TIME_STEP = (maxVelocity < std::numeric_limits<double>::epsilon()) ? MAX_TIME_STEP : std::min(MAX_TIME_STEP, 0.4 * SPACING / maxVelocity);
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

        g.position += TIME_STEP * g.velocity;
	}
}
