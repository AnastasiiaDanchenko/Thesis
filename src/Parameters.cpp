#include "..\headers\Parameters.h"

int WINDOW_WIDTH;
int WINDOW_HEIGHT;
int SCENE_DEPTH;
int PARTICLE_NEIGHBORS;

int PARTICLES_PER_DIMENSION;
float SPACING;
int CELL_SIZE;

float SUPPORT;
float REST_DENSITY;
float TIME_STEP;
float STIFFNESS;
float VISCOSITY;

Eigen::Vector3f GRAVITY = Eigen::Vector3f(0.0f, -9.8f, 0.0f);

std::string NS_METHOD;
std::string SIMULATION;

bool VISUALIZATION;

void readParameters() {
	std::ifstream file(fileName);

	if (!file.is_open()) {
		std::cout << "Error opening file " << fileName << std::endl;
		std::exit(1);
	}

	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string parameterName;
		std::string parameterValue;
		std::getline(iss, parameterName, ',');
		std::getline(iss, parameterValue, ',');
		if (parameterName == "window_width") {
			WINDOW_WIDTH = std::stoi(parameterValue);
		}
		else if (parameterName == "window_hight") {
			WINDOW_HEIGHT = std::stoi(parameterValue);
		}
		else if (parameterName == "depth") {
			SCENE_DEPTH = std::stoi(parameterValue);
		}
		else if (parameterName == "spacing") {
			SPACING = std::stod(parameterValue);
			CELL_SIZE = 2 * SPACING;
		}
		else if (parameterName == "support") {
			SUPPORT = SPACING * std::stod(parameterValue);
		}
		else if (parameterName == "density") {
			REST_DENSITY = std::stod(parameterValue);
		}
		else if (parameterName == "particles") {
			PARTICLES_PER_DIMENSION = std::stoi(parameterValue);
		}
		else if (parameterName == "timestep") {
			TIME_STEP = std::stod(parameterValue);
		}
		else if (parameterName == "stiffness") {
			STIFFNESS = std::stod(parameterValue);
		}
		else if (parameterName == "viscosity") {
			VISCOSITY = std::stod(parameterValue);
		}
		else if (parameterName == "neighbors") {
			PARTICLE_NEIGHBORS = std::stoi(parameterValue);
		}
	}
}
