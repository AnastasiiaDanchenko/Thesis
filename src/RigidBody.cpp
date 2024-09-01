#include "..\headers\RigidBody.h"

Eigen::Matrix3d vectorToMatrix(Eigen::Vector3d& vector) {
	Eigen::Matrix3d matrix;
	matrix << 0, -vector.z(), vector.y(),
		vector.z(), 0, -vector.x(),
		-vector.y(), vector.x(), 0;

	return matrix;
}

RigidBody::RigidBody() {
}

void RigidBody::initializeRigidBody(std::vector<Particle> particles, double bodyDensity) {
	this->innerParticles = particles;
	this->density = bodyDensity;
	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
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

RigidBody::RigidBody(std::vector<Particle> particles, double bodyDensity) {
	this->innerParticles = particles;
	this->density = bodyDensity;
	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
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

void RigidBody::initializeRigidBody(std::vector<Particle> inner, std::vector<Particle> outer, double bodyDensity) {
	this->innerParticles = inner;
	this->outerParticles = outer;
	this->density = bodyDensity;

	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
	this->invInitialInertiaTensor = Eigen::Matrix3d::Zero();
	this->rotationMatrix = Eigen::Matrix3d::Identity();

	for (auto& p : innerParticles) {
		this->positionCM += p.position * p.mass;
	}
	this->positionCM /= innerParticles.size() * parameters.restDensity * pow(parameters.spacing, 3);

	for (auto& p : innerParticles) {
		p.relativePosition = p.position - positionCM;
		this->invInertiaTensor -= p.mass * vectorToMatrix(p.relativePosition) * vectorToMatrix(p.relativePosition);
	}
	for (auto& p : outerParticles) {
		p.relativePosition = p.position - positionCM;
	}
	this->invInertiaTensor = this->invInertiaTensor.inverse();
	this->invInitialInertiaTensor = this->invInertiaTensor;
}

RigidBody::RigidBody(std::vector<Particle> inner, std::vector<Particle> outer, double bodyDensity) {
	this->innerParticles = inner;
	this->outerParticles = outer;
	this->density = bodyDensity;

	this->mass = innerParticles.size() * this->density * pow(parameters.spacing, 3);

	this->velocityCM = Eigen::Vector3d::Zero();
	this->positionCM = Eigen::Vector3d::Zero();
	this->angularVelocity = Eigen::Vector3d::Zero();
	this->invInertiaTensor = Eigen::Matrix3d::Identity();
	this->invInitialInertiaTensor = Eigen::Matrix3d::Zero();
	this->rotationMatrix = Eigen::Matrix3d::Identity();

	for (auto& p : innerParticles) {
		this->positionCM += p.position * p.mass;
	}
	this->positionCM /= innerParticles.size() * parameters.restDensity * pow(parameters.spacing, 3);

	for (auto& p : innerParticles) {
		p.relativePosition = p.position - positionCM;
		this->invInertiaTensor -= p.mass * vectorToMatrix(p.relativePosition) * vectorToMatrix(p.relativePosition);
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
		this->torque += (rotationMatrix * p.relativePosition).cross(p.mass * p.acceleration);
	}
}

void RigidBody::updateBodyQuantities() {
	// position first, using old velocity
	this->positionCM += this->velocityCM * parameters.timeStep;

	// velocity 
	this->velocityCM += parameters.timeStep * this->force / this->mass;
	
	// orientation A
	this->rotationMatrix += parameters.timeStep * vectorToMatrix(this->angularVelocity) * this->rotationMatrix;

	// angular momentum L
	this->angularMomentum += this->torque * parameters.timeStep;

	// inertia tensor I
	this->invInertiaTensor = this->rotationMatrix * this->invInitialInertiaTensor * this->rotationMatrix.transpose();
	
	// angular velocity w
	this->angularVelocity = this->invInertiaTensor * this->angularMomentum;

	
}

void RigidBody::updateParticles() {
	for (auto& p : outerParticles) {
		p.position = rotationMatrix * p.relativePosition + positionCM;
		p.velocity = angularVelocity.cross(rotationMatrix * p.relativePosition) + velocityCM;
		p.acceleration = Eigen::Vector3d::Zero();
	}
}

void RigidBody::discardInnerParticles() {
	innerParticles.clear();
}
