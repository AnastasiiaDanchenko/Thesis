#include "..\headers\RigidBody.h"

RigidBody::RigidBody() {
}

RigidBody::RigidBody(std::vector<Particle> particles) {
	this->innerParticles = particles;
	
	this->mass = innerParticles.size() * parameters.restDensity * pow(parameters.spacing, 3);
	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Zero();
	this->invInitialInertiaTensor = Eigen::Matrix3d::Zero();
	this->rotationMatrix = Eigen::Matrix3d::Identity();

	for (auto& p : innerParticles) {
		this->positionCM += p.position * p.mass;
	}
	this->positionCM /= mass;

	for (auto& p : innerParticles) {
		Eigen::Vector3d r = p.position - positionCM;
		this->invInertiaTensor -= p.mass * r * r.transpose();
	}
	this->invInertiaTensor = this->invInertiaTensor.inverse();
	this->invInitialInertiaTensor = this->invInertiaTensor;
}

RigidBody::RigidBody(std::vector<Particle> particles, Eigen::Vector3d position, Eigen::Vector3d velocity,
	Eigen::Matrix3d rotation, Eigen::Vector3d angularMomentum, Eigen::Matrix3d inertiaTensor, double mass) {
	this->outerParticles = particles;
	this->positionCM = position;
	this->velocityCM = velocity;
	this->rotationMatrix = rotation;
	this->angularMomentum = angularMomentum;
	this->mass = mass;
	this->invInertiaTensor = inertiaTensor;
}

void RigidBody::computeParticleQuantities() {
	this->force = this->mass * parameters.gravity;
	this->torque = Eigen::Vector3d::Zero();

	for (auto& p : outerParticles) {
		this->force += p.mass * p.acceleration;
		this->torque += (p.position - this->positionCM).cross(p.mass * p.acceleration);
	}
}

void RigidBody::updateBodyQuantities() {
	this->linearMomentum += this->force * parameters.timeStep;
	this->velocityCM = this->linearMomentum / this->mass;
	this->prevPositionCM = this->positionCM;
	this->positionCM += this->velocityCM * parameters.timeStep;

	this->invInertiaTensor = this->rotationMatrix * this->invInitialInertiaTensor * this->rotationMatrix.transpose();
	this->angularMomentum += this->torque * parameters.timeStep;
	this->angularVelocity = this->invInertiaTensor * this->angularMomentum;

	Eigen::Quaterniond deltaRotation(0, angularVelocity.x(), angularVelocity.y(), angularVelocity.z());
	deltaRotation *= this->orientation;
	deltaRotation.coeffs() *= 0.5 * parameters.timeStep;
	this->orientation.coeffs() += deltaRotation.coeffs();
	this->orientation.normalize();
	this->rotationMatrix = this->orientation.toRotationMatrix();
	//std::cout << "Rotation matrix: " << this->rotationMatrix << std::endl;
}

void RigidBody::updateParticles() {
	for (auto& p : outerParticles) {
		p.position = rotationMatrix * (p.position - prevPositionCM) + positionCM;
		p.velocity = angularVelocity.cross(p.position - prevPositionCM) + velocityCM;
	}
}

void RigidBody::discardInnerParticles() {
	innerParticles.clear();
}
