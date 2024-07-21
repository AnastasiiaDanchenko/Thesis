#include "..\headers\Solvers.h"

Solver2D::Solver2D() {
	grid2D = Grid2D(parameters.spacing * 2);
}

void Solver2D::initBoundaries() {
    int width = (parameters.windowSize.width / 2) / parameters.spacing - 1;
    int hight = (parameters.windowSize.height / 2) / parameters.spacing - 1;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1) {
                Particle2D p;

                p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
                p.isFluid = false;
                p.ID = particles2D.size();

                particles2D.push_back(p);
            }
        }
    }
}

void Solver2D::initFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing);
            p.ID = particles2D.size();
            particles2D.push_back(p);
        }
    }
}

void Solver2D::neighborSearch() {
    grid2D.updateGrid();
    grid2D.neighborSearch(particles2D);
}

void Solver2D::computeDensityPressure() {
    double averageDensity = 0.0f;

    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }

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

void Solver2D::computeAcceleration() {
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }

        // Gravity force
        Eigen::Vector2d acceleration = parameters.gravity2D;

        for (auto neighbor : p.neighbors) {
            if (neighbor == &p) { continue; }
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

void Solver2D::updateParticles() {
    for (auto& p : particles2D) {
        if (p.isFluid == false) { continue; }
        p.velocity += parameters.timeStep * p.acceleration;
        p.position += parameters.timeStep * p.velocity;
    }
}

void Solver2D::computeDensity() {
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

void Solver2D::computeSurface() {
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

void Solver2D::predictVelocity() {
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

void Solver2D::computeDensityError() {
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

void Solver2D::computeLaplacian() {
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

void Solver2D::compressionConvergence() {
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

void Solver2D::advectParticles() {
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

void Solver2D::boundaryMassUpdate() {
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

Eigen::Vector2d rotatePlane(const Eigen::Vector2d& center, const Eigen::Vector2d& pos0, int angle) {
    double angleRad = angle * M_PI / 180;
    double x = center.x() + cos(angleRad) * (pos0.x() - center.x()) - sin(angleRad) * (pos0.y() - center.y());
    double y = center.y() + sin(angleRad) * (pos0.x() - center.x()) + cos(angleRad) * (pos0.y() - center.y());

    return Eigen::Vector2d(x, y);
}

void Solver2D::initRotatingBoundary() {
    int width_max = parameters.windowSize.width / parameters.spacing / 5;
    int width_avg = parameters.windowSize.width / parameters.spacing / 10;

    int height_avg = parameters.windowSize.height / parameters.spacing / 7;

    const Eigen::Vector2d center = Eigen::Vector2d((width_max + width_avg) * parameters.spacing, height_avg * parameters.spacing);
    const double angularVelocity = 0.05;

    for (int i = 0; i < width_max; i++) {
        Particle2D p0, p1, p2, p3;

        if (i == width_avg) {
            p0.position = Eigen::Vector2d((i + width_max + 2.5) * parameters.spacing, (height_avg + 1) * parameters.spacing);
            p0.isFluid = false;
            p0.ID = particles2D.size();

            particles2D.push_back(p0);

            p1.position = rotatePlane(center, p0.position, 45);
            p1.isFluid = false;
            p1.ID = particles2D.size();

            particles2D.push_back(p1);

            p2.position = rotatePlane(center, p0.position, 90);
            p2.isFluid = false;
            p2.ID = particles2D.size();

            particles2D.push_back(p2);

            p3.position = rotatePlane(center, p0.position, 135);
            p3.isFluid = false;
            p3.ID = particles2D.size();

            particles2D.push_back(p3);

            p0.position = Eigen::Vector2d((i + width_max - 2.5) * parameters.spacing, (height_avg - 1) * parameters.spacing);
        }
        else {
            p0.position = Eigen::Vector2d((i + width_max) * parameters.spacing, height_avg * parameters.spacing);
        }
        p0.isFluid = false;
        p0.ID = particles2D.size();

        particles2D.push_back(p0);

        p1.position = rotatePlane(center, p0.position, 45);
        p1.isFluid = false;
        p1.ID = particles2D.size();

        particles2D.push_back(p1);

        p2.position = rotatePlane(center, p0.position, 90);
        p2.isFluid = false;
        p2.ID = particles2D.size();

        particles2D.push_back(p2);

        p3.position = rotatePlane(center, p0.position, 135);
        p3.isFluid = false;
        p3.ID = particles2D.size();

        particles2D.push_back(p3);

        if (i == 0) {
            parameters.boundaryTestID = p0.ID;
        }
    }
}

void Solver2D::rotateBoundary() {
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

void Solver2D::initMovingBoundary() {
    int width = (parameters.windowSize.width / 2) / parameters.spacing - 1;
    int hight = (parameters.windowSize.height / 2) / parameters.spacing - 1;

    for (int i = 5; i < width - 5; i++) {
        for (int j = 5; j < hight - 5; j++) {
            if (i < 6 || i > width - 7 || j < 6 || j > hight - 7) {
                Particle2D p;

                p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
                p.isFluid = false;
                p.ID = particles2D.size();

                particles2D.push_back(p);
            }
        }
    }
}

void Solver2D::initMovingFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle2D p;
            p.position = Eigen::Vector2d((i + 7) * parameters.spacing, (j + 7) * parameters.spacing);
            p.ID = particles2D.size();
            particles2D.push_back(p);
        }
    }
}

void Solver2D::moveBoundary() {
    for (int i = -((parameters.windowSize.width / 2) / parameters.spacing - 1) / 2; i < 0; i++) {
        for (int j = 3; j < 6; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * parameters.spacing, j * parameters.spacing + parameters.windowSize.height / 10);
            p.isFluid = false;
            p.ID = particles2D.size();
            p.velocity = Eigen::Vector2d(20, 0.0);

            particles2D.push_back(p);
        }
    }
}

Solver::Solver() {
    grid = Grid(parameters.spacing * 2);
}

void Solver::initBoundaries() {
    int width = parameters.windowSize.width / parameters.spacing - 1;
    int hight = parameters.windowSize.height / parameters.spacing - 1;
    int depth = (parameters.windowSize.depth / parameters.spacing - 1) / 2;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            for (float k = 0; k < depth; k += 0.5) {
                if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1 || k < 0.5 || k > depth - 1) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, (k + 1) * parameters.spacing);
                    p.isFluid = false;
                    p.ID = particles.size();
                    parameters.boundaryTestID = p.ID;

                    particles.push_back(p);
                }
            }
        }
    }
}

void Solver::initFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            for (int k = 0; k < parameters.particlesPerDimension.z; k++) {
                Particle p;

                p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, (k + 2) * parameters.spacing);
                p.ID = particles.size();
                particles.push_back(p);
            }
        }
    }
}

void Solver::initGhostFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle p;

            p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, parameters.windowSize.depth / 4);
            p.ID = ghostParticles.size();
            ghostParticles.push_back(p);
        }
    }
}

void Solver::initGhostBoundary() {
    int width = parameters.windowSize.width / parameters.spacing - 1;
    int hight = parameters.windowSize.height / parameters.spacing - 1;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1) {
                Particle p;

                p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, parameters.windowSize.depth / 4);
                p.isFluid = false;
                p.ID = ghostParticles.size();
                ghostParticles.push_back(p);
            }
        }
    }
}

void Solver::initMovingBoundary() {
    for (int i = -20; i < 0; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = 0; k < (parameters.windowSize.depth / parameters.spacing - 1) / 2; k++) {
                Particle p;

                p.position = Eigen::Vector3d(i * parameters.spacing, j * parameters.spacing, (k + 1) * parameters.spacing);
                p.isFluid = false;
                p.ID = particles.size();
                p.velocity = Eigen::Vector3d(20, 0.0, 0.0);

                particles.push_back(p);
            }
        }
    }
}

void Solver::neighborSearch() {
    std::vector<Particle*> allParticlePtrs;
    for (auto& particle : particles) {
        allParticlePtrs.push_back(&particle);
    }
    for (auto& body : rigidBodies) {
        for (auto& particle : body.getOuterParticles()) {
            allParticlePtrs.push_back(&particle);
        }
    }

	grid.updateGrid(allParticlePtrs);
	grid.neighborSearch(allParticlePtrs);
}

void Solver::computeDensityPressure(){
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

void Solver::computeAcceleration() {
    for (auto& p : particles) {
        if (p.isFluid == false) { continue; }

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

void Solver::updateParticles() {
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
        if (p.position.x() < -parameters.windowSize.width || p.position.x() > parameters.windowSize.width * 2 ||
            p.position.y() < -parameters.windowSize.height || p.position.y() > parameters.windowSize.height * 2) {
            particles.erase(particles.begin() + i);
        }
    }

    parameters.timeStep = (maxVelocity < std::numeric_limits<double>::epsilon()) ? parameters.maxTimeStep :
        std::min(parameters.maxTimeStep, 0.4 * parameters.spacing / maxVelocity);
}

void Solver::computeDensity() {
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

void Solver::predictVelocity() {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        if (p.isFluid == false) { continue; }

        p.acceleration = parameters.gravity + ViscosityAcceleration(&p);
        p.predictedVelocity = p.velocity + parameters.timeStep * p.acceleration;
    }
}

void Solver::computeDensityError() {
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

void Solver::computeLaplacian() {
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

void Solver::compressionConvergence() {
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

void Solver::boundaryMassUpdate() {
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

void Solver::neighborSearchGhosts() {
	grid.neighborSearchGhosts();
}

void Solver::updateGhosts() {
#pragma omp parallel for
    for (int i = 0; i < ghostParticles.size(); i++) {
        if (!ghostParticles[i].isFluid) { continue; }
		Particle& g = ghostParticles[i];
		
        double density = 0, pressure = 0, weight = 0;
        Eigen::Vector3d velocity = Eigen::Vector3d::Zero();

        for (auto gn : g.neighbors) {
            // weighted sum of the density of the neighbors depending on the distance
            double distance = (g.position - gn->position).norm();
            weight += 1 / distance;

            density += 1 / distance * gn->density;
            pressure += 1 / distance * gn->pressure;
            velocity += 1 / distance * gn->velocity;
		}
        
		g.density = density / weight;
        g.pressure = pressure / weight;
		g.velocity.x() = velocity.x() / weight;
        g.velocity.y() = velocity.y() / weight;

        g.position += parameters.timeStep * g.velocity;
	}
}

void Solver::initRigidCube() {
    int depth = (parameters.windowSize.depth / parameters.spacing - 1) / 4;
	int width = depth, height = depth;

	std::vector<Particle> body, contour;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < depth; k++) {
				Particle p;

				p.position = Eigen::Vector3d(
                    (i + 2) * parameters.spacing,
                    (j + parameters.windowSize.depth / parameters.spacing - 2) * parameters.spacing,
                    (k + 2) * parameters.spacing
                );
				p.isFluid = false;
				body.push_back(p);

				if (i == 0 || i == width - 1 || j == 0 || j == height - 1 || k == 0 || k == depth - 1) {
					contour.push_back(p);
				}
			}
		}
	}

	RigidBody newBody(body);
    newBody.discardInnerParticles();
	newBody.setOuterParticles(contour);
    
	this->rigidBodies.push_back(newBody);
}
