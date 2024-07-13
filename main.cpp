#include "headers/Visualization.h"

// Close the window when pressing ESC
void keyCallback_old(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
}

void SEvsPPE() {
    std::vector<std::string> METHODS = { "SE", "PPE" };
    int steps = 100;

    for (auto method : METHODS) {
        Initialization2D();

        auto start = std::chrono::high_resolution_clock::now();
        if (method == "SE") {
            TIME_STEP = 0.000001;
            for (int i = 0; i < steps; i++) {
                Simulation2D();
            }
            //SaveToDisk2D();
        }
        else if (method == "PPE") {
            TIME_STEP = 0.05;
            for (int i = 0; i < steps; i++) {
                SimulationIISPH2D();
            }
            //SaveToDisk2D();
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        std::cout << "Method: " << method
                << "\nNumber of particles: " << PARTICLES_X * PARTICLES_Y
                << "\nTime step: " << TIME_STEP
                << "\nAverage density by the end of simulation: " << AVG_DENSITY << " kg/m^3"
                << "\nExecution time total: " << duration.count() << "s" << std::endl;
        std::cout << "Execution time per step: " << duration.count() / steps << "s" << std::endl << std::endl;

        particles2D.clear();
        grid2D.clear();
    }
}

int main() {
    // Read parameters from the input file
    readParameters();

    if (DIMENSIONS == 2){
        Initialization2D();
        Visualize2D();
        //SEvsPPE();
    }
    else if (DIMENSIONS == 3) {
        Initialization();
        //Visualize();
        VisualizeGhosts();
	}

    return 0;
}
