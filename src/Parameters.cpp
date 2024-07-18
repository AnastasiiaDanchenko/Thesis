#include "..\headers\Parameters.h"

using json = nlohmann::json;

Parameters::Parameters() : 
    gravity(0.0, -9.8, 0.0),
    gravity2D(0.0, -9.8) {}

Parameters parameters;

void Parameters::readParameters() {
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cout << "Error opening file " << fileName << std::endl;
        std::exit(1);
    }

    json j;
    file >> j;

    if (j.contains("dimensions")) dimensions = j["dimensions"].get<int>();
    if (j.contains("support")) support = j["support"].get<double>();
    if (j.contains("density")) restDensity = j["density"].get<double>();
    if (j.contains("visualize neighbors")) visualizeNeighbors = j["visualize neighbors"].get<int>();

    std::string dimension_key = (dimensions == 2) ? "2D" : "3D";

    if (j.contains(dimension_key)) {
        json dim_data = j[dimension_key];

        if (dim_data.contains("window")) {
            auto window = dim_data["window"];
            windowSize.width = window[0].get<int>();
            windowSize.height = window[1].get<int>();
            if (dimensions == 3 && window.size() > 2) {
                windowSize.depth = window[2].get<int>();
            }
        }

        if (dim_data.contains("spacing")) {
            spacing = dim_data["spacing"].get<double>();
            support *= spacing;
        }

        if (dim_data.contains("particles nb")) {
            auto particles_nb = dim_data["particles nb"];
            particlesPerDimension.x = particles_nb[0].get<int>();
            particlesPerDimension.y = particles_nb[1].get<int>();
            if (dimensions == 3 && particles_nb.size() > 2) {
                particlesPerDimension.z = particles_nb[2].get<int>();
            }
        }

        if (dim_data.contains("timestep")) timeStep = dim_data["timestep"].get<double>();
        if (dim_data.contains("stiffness")) stiffness = dim_data["stiffness"].get<double>();
        if (dim_data.contains("viscosity")) viscosity = dim_data["viscosity"].get<double>();
        if (dim_data.contains("gamma")) gamma = dim_data["gamma"].get<double>();
        if (dim_data.contains("omega")) omega = dim_data["omega"].get<double>();
        if (dim_data.contains("error, %")) errThreshold = dim_data["error, %"].get<double>() * 0.01;

        if (dim_data.contains("surface_tension")) {
            json surface_tension = dim_data["surface_tension"];
            if (surface_tension.contains("enabled")) surfaceTension = surface_tension["enabled"].get<bool>();
            if (surface_tension.contains("cohesion")) cohesion = surface_tension["cohesion"].get<double>();
        }
    }

    if (dimensions == 2) {
        maxTimeStep = 0.005 * spacing;
    }
    else if (dimensions == 3) {
        maxTimeStep = 0.0025 * spacing;
    }
}
