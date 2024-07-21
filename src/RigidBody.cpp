#include "..\headers\RigidBody.h"

RigidBody::RigidBody() {
}

RigidBody::RigidBody(std::vector<Particle> particles) {
	innerParticles = particles;
	
	mass = 0;
	positionCM = Eigen::Vector3d::Zero();
	inertiaTensor = Eigen::Matrix3d::Zero();

	for (auto& p : innerParticles) {
		mass += p.mass;
		positionCM += p.position * p.mass;
	}
	positionCM /= mass;

	for (auto& p : innerParticles) {
		Eigen::Vector3d r = p.position - positionCM;
		inertiaTensor += p.mass * (r.squaredNorm() * Eigen::Matrix3d::Identity() - r * r.transpose());
		//inertiaTensor -= p.mass * (p.position - positionCM).transpose() * (p.position - positionCM).transpose();
	}
}

RigidBody::RigidBody(std::vector<Particle> particles, Eigen::Vector3d position, Eigen::Vector3d velocity,
	Eigen::Matrix3d rotation, Eigen::Vector3d angularMomentum, Eigen::Matrix3d inertiaTensor, double mass, 
	double density) {
	outerParticles = particles;
	positionCM = position;
	velocityCM = velocity;
	rotationMatrix = rotation;
	this->angularMomentum = angularMomentum;
	this->mass = mass;
	this->inertiaTensor = inertiaTensor;
	this->density = density;
}

void RigidBody::computeParticleQuantities() {
	Eigen::Vector3d force = Eigen::Vector3d::Zero(), torque = Eigen::Vector3d::Zero();
	/*for (auto& p : outerParticles) {
		torque += (p.position - positionCM).cross(p.acceleration / p.mass);
		force += p.acceleration / p.mass;
	}

	this->torque = torque;
	this->force = force;*/

	for (auto& p : outerParticles) {
		auto acceleration = parameters.gravity;
		//Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
		for (auto n : p.neighbors) {
			acceleration += n->acceleration + n->pressureAcceleration;// -parameters.gravity;
		}
		force += acceleration;
		torque += (p.position - positionCM).cross(acceleration);
	}
	this->force = force;
	this->torque = torque;
}

void RigidBody::updateBodyQuantities() {
	positionCM += velocityCM * parameters.timeStep;
	velocityCM += force * parameters.timeStep;

	Eigen::Quaterniond deltaRotation;
	Eigen::Vector3d deltaTheta = angularVelocity * parameters.timeStep;
	double angle = deltaTheta.norm();
	if (angle > 1e-12) {
		Eigen::Vector3d axis = deltaTheta.normalized();
		deltaRotation = Eigen::AngleAxisd(angle, axis);
	}
	else {
		deltaRotation = Eigen::Quaterniond(1, deltaTheta.x() / 2, deltaTheta.y() / 2, deltaTheta.z() / 2);
	}
	orientation = (deltaRotation * orientation).normalized();
	rotationMatrix = orientation.toRotationMatrix();
	//rotationMatrix += angularVelocity * rotationMatrix * parameters.timeStep;

	angularMomentum += torque * parameters.timeStep;
	inertiaTensor = rotationMatrix * inertiaTensor * rotationMatrix.transpose();
	angularVelocity = inertiaTensor * angularMomentum;
}

void RigidBody::updateParticles() {
	for (auto& p : outerParticles) {
		auto positionPrev = p.position;
		auto velocityPrev = p.velocity;
		p.position = rotationMatrix * (p.position - positionCM) + positionCM;
		p.velocity = angularVelocity.cross(p.position - positionCM) + velocityCM;

		if (p.position != positionPrev) {
			std::cout << "Position changed" << std::endl;
		}
		if (p.velocity != velocityPrev) {
			std::cout << "Velocity changed" << std::endl;
		}
	}
}

void RigidBody::discardInnerParticles() {
	innerParticles.clear();
}
