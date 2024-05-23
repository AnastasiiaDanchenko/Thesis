#include "..\headers\ComputationsIISPH.h"

void CalculateDensity2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		float density = 0.0f;

		for (auto n : p.neighbors) {
			density += n->mass * CubicSplineKernel2D(p.position - n->position);
		}

		p.density = density;
	}
}

void PredictVelocityAdvection2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		Eigen::Vector2f acceleration = GRAVITY2D;

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

void ComputeDii() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		Eigen::Vector2f dii = Eigen::Vector2f::Zero();

		for (auto neighbor : p.neighbors) {
			dii -= neighbor->mass * CubicSplineKernelGradient2D(p.position - neighbor->position);
		}
		dii *= pow(TIME_STEP, 2) / pow(p.density, 2);
		p.dii = dii;
	}
}

void ComputeAii() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		float aii = 0.0f;

		for (auto neighbor : p.neighbors) {
			aii += p.mass * (p.dii + p.mass * pow(TIME_STEP, 2) / pow(p.density, 2) *
				CubicSplineKernelGradient2D(neighbor->position - p.position)).dot(
				CubicSplineKernelGradient2D(p.position - neighbor->position));
		}
		p.aii = aii;
	}
}

void PredictDensity2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		float density = p.density;

		for (auto neighbor : p.neighbors) {
			density += TIME_STEP * neighbor->mass * (p.predictedVelocity - neighbor->predictedVelocity).dot(
				CubicSplineKernelGradient2D(p.position - neighbor->position));
		}
		p.predictedDensity = density;
	}
}

void SolvePressure2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		p.pressure *= 0.5;
	}

	NB_ITERATIONS = 0;
	bool conversion = true;

	while (conversion || NB_ITERATIONS < 2) {
#pragma omp parallel for
		for (int i = 0; i < particles2D.size(); i++) {

			if (particles2D[i].isFluid == false) {
				continue;
			}

			Particle2D& p = particles2D[i];
			Eigen::Vector2f ci = Eigen::Vector2f::Zero();

			for (auto neighbor : p.neighbors) {
				ci -= neighbor->pressure / pow(neighbor->density, 2) *
					CubicSplineKernelGradient2D(p.position - neighbor->position);
			}
			ci *= p.mass * pow(TIME_STEP, 2);
			p.ci = ci;
		}

		AVG_DENSITY = 0.0f;
#pragma omp parallel for 
		for (int i = 0; i < particles2D.size(); i++) {

			if (particles2D[i].isFluid == false) {
				continue;
			}

			Particle2D& p = particles2D[i];
			float newPressure = 0.0f;

			if (p.aii != 0) {
				p.predictedPressure = REST_DENSITY - p.predictedDensity;

				for (auto neighbor : p.neighbors) {
					if (!neighbor->isFluid) {
						p.predictedPressure -= neighbor->mass * p.ci.dot(
							CubicSplineKernelGradient2D(p.position - neighbor->position));
					}
					else {
						p.predictedPressure -= neighbor->mass * (p.ci - neighbor->dii * neighbor->pressure - neighbor->ci - 
							p.mass * pow(TIME_STEP, 2) / pow(p.density, 2) * CubicSplineKernelGradient2D(neighbor->position - p.position) *
							p.pressure).dot(CubicSplineKernelGradient2D(p.position - neighbor->position));
					}
				}
				newPressure = std::max(0.0f, (1 - OMEGA) * p.pressure + OMEGA / p.aii * p.predictedPressure);
			}

			p.pressure = newPressure;
			AVG_DENSITY += p.density;
		}

		NB_ITERATIONS++;
		AVG_DENSITY /= PARTICLES_X * PARTICLES_Y;
		DENSITY_ERR = abs(AVG_DENSITY - REST_DENSITY);

		if (DENSITY_ERR <= ERR_THRESHOLD * REST_DENSITY) {
			conversion = false;
		}

		if (NB_ITERATIONS == 100) { // Maximum nb of iterations
			conversion = false;
		}
	}
}

void ComputePressureAcceleration2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		Eigen::Vector2f pressureAcceleration = Eigen::Vector2f::Zero();

		for (auto neighbor : p.neighbors) {
			pressureAcceleration -= neighbor->mass * (p.pressure / pow(p.density, 2) + neighbor->pressure / pow(neighbor->density, 2)) *
				CubicSplineKernelGradient2D(p.position - neighbor->position);
		}
		p.pressureAcceleration = pressureAcceleration;
	}
}

void AdvectParticles2D() {
#pragma omp parallel for
	for (int i = 0; i < particles2D.size(); i++) {

		if (particles2D[i].isFluid == false) {
			continue;
		}

		Particle2D& p = particles2D[i];
		p.velocity = p.predictedVelocity + TIME_STEP * p.pressureAcceleration;
		p.position += TIME_STEP * p.velocity;
	}
}
