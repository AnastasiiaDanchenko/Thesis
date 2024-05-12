#include "..\headers\Parameters.h"

int WINDOW_WIDTH = 800;
int WINDOW_HEIGHT = 800;
int SCENE_DEPTH = 600;
int PARTICLE_NEIGHBORS = 0;

int PARTICLES_PER_DIMENSION = 10;
int PARTICLES_X = 10;
int PARTICLES_Y = 10;
int PARTICLES_Z = 10;
float SPACING = 10.0f;
int CELL_SIZE;

float SUPPORT = 2.1 * SPACING;
float REST_DENSITY = 1000.0f;
float TIME_STEP = 0.01f;
float STIFFNESS = 100000.0f;
float VISCOSITY = 0.1f;

Eigen::Vector3f GRAVITY = Eigen::Vector3f(0.0f, -9.8f, 0.0f);
Eigen::Vector2f GRAVITY2D = Eigen::Vector2f(0.0f, -9.8f);

float GAMMA = 0.7f;
float OMEGA = 0.5f;
float AVG_DENSITY = 0.0f;
float DENSITY_ERR = 0.0f;
float ERR_THRESHOLD = 0.01f;

std::string NS_METHOD;
std::string SIMULATION;

bool VISUALIZATION = true;
int DIMENSIONS = 3;

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
		else if (parameterName == "particlesX") {
			PARTICLES_X = std::stoi(parameterValue);
		}
		else if (parameterName == "particlesY") {
			PARTICLES_Y = std::stoi(parameterValue);
		}
		else if (parameterName == "particlesZ") {
			PARTICLES_Z = std::stoi(parameterValue);
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
		else if (parameterName == "gamma") {
			GAMMA = std::stod(parameterValue);
		}
		else if (parameterName == "omega") {
			OMEGA = std::stod(parameterValue);
		}
		else if (parameterName == "dimensions") {
			DIMENSIONS = std::stoi(parameterValue);
		}
		else if (parameterName == "error") {
			ERR_THRESHOLD = std::stod(parameterValue);
		}
	}
}
