#include "..\headers\Parameters.h"

int WINDOW_WIDTH = 800;
int WINDOW_HEIGHT = 800;
int SCENE_DEPTH = 600;
int PARTICLE_NEIGHBORS = 0;

int PARTICLES_PER_DIMENSION = 10;
int PARTICLES_X = 10;
int PARTICLES_Y = 10;
int PARTICLES_Z = 10;
double SPACING = 10.0;
double CELL_SIZE;

double SUPPORT = 2.0 * SPACING;
double REST_DENSITY = 1000.0f;
double TIME_STEP = 0.01f;
double STIFFNESS = 100000.0f;
double VISCOSITY = 0.01;
double COHESION = 0.00001f;
double MAX_TIME_STEP = 0.05f;
int ITERATIONS_COUNT = 0;

Eigen::Vector3d GRAVITY = Eigen::Vector3d(0.0f, -9.8f, 0.0f);
Eigen::Vector2d GRAVITY2D = Eigen::Vector2d(0.0, -9.8);

double GAMMA = 0.7f;
double OMEGA = 0.5f;
double AVG_DENSITY = 0.0f;
double DENSITY_ERR = 0.0f;
double FIRST_ERR = 0.0f;
double ERR_THRESHOLD = 0.001f;
int NB_ITERATIONS = 0;

std::string NS_METHOD;
std::string SIMULATION;

bool VISUALIZATION = true;
bool SURFACE_TENSION = true;
int DIMENSIONS = 3;
bool FIRST_STEP = true;
double FIRST_STEP_CORRECTION = 1.0;

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
		else if (parameterName == "max_timestep") {
			MAX_TIME_STEP = std::stod(parameterValue);
			MAX_TIME_STEP = 0.005 * SPACING;
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
		else if (parameterName == "error%") {
			ERR_THRESHOLD = std::stod(parameterValue) * 0.01; // Convert to percentage
		}
		else if (parameterName == "surface_tension") {
			if (parameterValue == "0") {
				SURFACE_TENSION = false;
			}
			else {
				SURFACE_TENSION = true;
			}
		}
		else if (parameterName == "cohesion") {
			COHESION = std::stod(parameterValue);
		}
	}

	if (DIMENSIONS == 2) {
		MAX_TIME_STEP = 0.005 * SPACING;
	}
	else if (DIMENSIONS == 3) {
		MAX_TIME_STEP = 0.0025 * SPACING;
	}
}
