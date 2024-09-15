#pragma once
#include "Particle.h"
#include <map>

class RigidBody {
private:
	std::vector<Particle> innerParticles;
	std::vector<Particle> outerParticles;
	Eigen::Vector3d positionCM; // X
	Eigen::Vector3d velocityCM; // V
	Eigen::Vector3d angularVelocity; // omega
	Eigen::Vector3d force; // F
	Eigen::Vector3d torque; // tau
	Eigen::Vector3d angularMomentum; // L
	Eigen::Matrix3d rotationMatrix; // R
	Eigen::Matrix3d invInertiaTensor; // ^ -1
	Eigen::Matrix3d invInitialInertiaTensor; // ^ -1
	double mass;
	double density;
	
public:
	RigidBody();
	RigidBody(std::vector<Particle> particles, double density);
	RigidBody(std::vector<Particle> innerParticles, std::vector<Particle> outerParticles, double density);
	RigidBody(std::vector<Particle> particles, Eigen::Vector3d position, Eigen::Vector3d velocity, 
		Eigen::Matrix3d rotation, Eigen::Vector3d angularMomentum, Eigen::Matrix3d inertiaTensor, double mass);

	void initializeRigidBody(std::vector<Particle> particles, double density);
	void initializeRigidBody(std::vector<Particle> inner, std::vector<Particle> outer, double bodyDensity);
	void discardInnerParticles();

	void computeParticleQuantities();
	void updateBodyQuantities();
	void updateParticles();

	void setOuterParticles(std::vector<Particle> particles) { this->outerParticles = particles; }
	void setConstantVelocity(Eigen::Vector3d velocity) { this->velocityCM = velocity; }

	std::vector<Particle>& getOuterParticles() { return this->outerParticles; }
	Eigen::Vector3d& getPositionCM() { return this->positionCM; }
	double getMass() { return this->mass; }
	Eigen::Matrix3d getInvInertiaTensor() { return this->invInertiaTensor; }
	Eigen::Vector3d getVelocityCM() { return this->velocityCM; }
	Eigen::Vector3d getAngularVelocity() { return this->angularVelocity; }
};
