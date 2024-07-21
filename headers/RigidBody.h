#pragma once
#include "Particle.h"

class RigidBody {
private:
	std::vector<Particle> innerParticles;
	std::vector<Particle> outerParticles;
	Eigen::Vector3d positionCM;
	Eigen::Vector3d velocityCM;
	Eigen::Vector3d angularVelocity;
	Eigen::Vector3d force;
	Eigen::Vector3d torque;
	Eigen::Vector3d angularMomentum;
	Eigen::Matrix3d rotationMatrix;
	Eigen::Matrix3d inertiaTensor; // ^ -1
	Eigen::Quaterniond orientation;
	double mass;
	double density;
	
public:
	RigidBody();
	RigidBody(std::vector<Particle> particles);
	RigidBody(std::vector<Particle> particles, Eigen::Vector3d position, Eigen::Vector3d velocity, 
		Eigen::Matrix3d rotation, Eigen::Vector3d angularMomentum, Eigen::Matrix3d inertiaTensor, double mass, 
		double density);
	void discardInnerParticles();
	void computeParticleQuantities();
	void updateBodyQuantities();
	void updateParticles();

	void setOuterParticles(std::vector<Particle> particles) { outerParticles = particles; }

	std::vector<Particle>& getOuterParticles() { return outerParticles; }
};
