#include "..\headers\RigidBody.h"

RigidBody::RigidBody() {
}

RigidBody::RigidBody(std::vector<Particle> particles, double density) {
	this->innerParticles = particles;
	this->density = density;
	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
	this->invInitialInertiaTensor = Eigen::Matrix3d::Zero();
	this->rotationMatrix = Eigen::Matrix3d::Identity();
	this->orientation = Eigen::Quaterniond::Identity();

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

RigidBody::RigidBody(std::vector<Particle> inner, std::vector<Particle> outer, double density) {
	this->innerParticles = inner;
	this->outerParticles = outer;
	this->density = density;

	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
	this->invInitialInertiaTensor = Eigen::Matrix3d::Zero();
	this->rotationMatrix = Eigen::Matrix3d::Identity();
	this->orientation = Eigen::Quaterniond::Identity();

	for (auto& p : innerParticles) {
		this->positionCM += p.position * p.mass;
	}
	this->positionCM /= mass;

	for (auto& p : innerParticles) {
		p.relativePosition = p.position - positionCM;
		this->invInertiaTensor -= p.mass * p.relativePosition * p.relativePosition.transpose();
	}
	for (auto& p : outerParticles) {
		p.relativePosition = p.position - positionCM;
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
	// position first, using old velocity
	this->prevPositionCM = this->positionCM;
	this->positionCM += this->velocityCM * parameters.timeStep;

	// velocity 
	this->linearMomentum += this->force * parameters.timeStep;
	this->velocityCM = this->linearMomentum / this->mass;
	
	// orientation A
	Eigen::Quaterniond deltaRotation(0, angularVelocity.x(), angularVelocity.y(), angularVelocity.z());
	deltaRotation *= this->orientation;
	deltaRotation.coeffs() *= 0.5 * parameters.timeStep;
	this->orientation.coeffs() += deltaRotation.coeffs();
	this->orientation.normalize();
	//std::cout << "Orientation: " << orientation.coeffs() << std::endl;
	this->rotationMatrix = this->orientation.toRotationMatrix();

	// angular momentum L
	this->angularMomentum += this->torque * parameters.timeStep;

	// inertia tensor I
	this->invInertiaTensor = this->rotationMatrix * this->invInitialInertiaTensor * this->rotationMatrix.transpose();
	
	// angular velocity w
	this->angularVelocity = this->invInertiaTensor * this->angularMomentum;
	//std::cout << "Angular velocity: " << angularVelocity.transpose() << std::endl;

	
}

void RigidBody::updateParticles() {
	for (auto& p : outerParticles) {
		p.position = rotationMatrix * p.relativePosition + positionCM;
		p.velocity = angularVelocity.cross(p.relativePosition) + velocityCM;
	}
}

void RigidBody::discardInnerParticles() {
	innerParticles.clear();
}
