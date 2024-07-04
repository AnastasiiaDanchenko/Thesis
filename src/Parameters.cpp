#include "..\headers\Parameters.h"

using json = nlohmann::json;

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

    json j;
    file >> j;

    if (j.contains("dimensions")) DIMENSIONS = j["dimensions"].get<int>();
    if (j.contains("support")) SUPPORT = SPACING * j["support"].get<double>();
    if (j.contains("density")) REST_DENSITY = j["density"].get<double>();
    if (j.contains("visualize neighbors")) PARTICLE_NEIGHBORS = j["visualize neighbors"].get<int>();

    std::string dimension_key = (DIMENSIONS == 2) ? "2D" : "3D";

    if (j.contains(dimension_key)) {
        json dim_data = j[dimension_key];

        if (dim_data.contains("window")) {
            auto window = dim_data["window"];
            WINDOW_WIDTH = window[0].get<int>();
            WINDOW_HEIGHT = window[1].get<int>();
            if (DIMENSIONS == 3 && window.size() > 2) {
                SCENE_DEPTH = window[2].get<int>();
            }
        }

        if (dim_data.contains("spacing")) {
            SPACING = dim_data["spacing"].get<double>();
            CELL_SIZE = 2 * SPACING;
        }

        if (dim_data.contains("particles nb")) {
            auto particles_nb = dim_data["particles nb"];
            PARTICLES_X = particles_nb[0].get<int>();
            PARTICLES_Y = particles_nb[1].get<int>();
            if (DIMENSIONS == 3 && particles_nb.size() > 2) {
                PARTICLES_Z = particles_nb[2].get<int>();
            }
        }

        if (dim_data.contains("timestep")) TIME_STEP = dim_data["timestep"].get<double>();
        if (dim_data.contains("stiffness")) STIFFNESS = dim_data["stiffness"].get<double>();
        if (dim_data.contains("viscosity")) VISCOSITY = dim_data["viscosity"].get<double>();
        if (dim_data.contains("gamma")) GAMMA = dim_data["gamma"].get<double>();
        if (dim_data.contains("omega")) OMEGA = dim_data["omega"].get<double>();
        if (dim_data.contains("error, %")) ERR_THRESHOLD = dim_data["error, %"].get<double>() * 0.01;

        if (dim_data.contains("surface_tension")) {
            json surface_tension = dim_data["surface_tension"];
            if (surface_tension.contains("enabled")) SURFACE_TENSION = surface_tension["enabled"].get<bool>();
            if (surface_tension.contains("cohesion")) COHESION = surface_tension["cohesion"].get<double>();
        }
    }

    if (DIMENSIONS == 2) {
        MAX_TIME_STEP = 0.005 * SPACING;
    }
    else if (DIMENSIONS == 3) {
        MAX_TIME_STEP = 0.0025 * SPACING;
    }
}
