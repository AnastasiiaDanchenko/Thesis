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
            parameters.timeStep = 0.000001;
            for (int i = 0; i < steps; i++) {
                //Simulation2D();
            }
            //SaveToDisk2D();
        }
        else if (method == "PPE") {
            parameters.timeStep = 0.05;
            for (int i = 0; i < steps; i++) {
                //SimulationIISPH2D(grid2D);
            }
            //SaveToDisk2D();
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        std::cout << "Method: " << method
			<< "\nNumber of particles: " << parameters.particlesPerDimension.x * parameters.particlesPerDimension.y
            << "\nTime step: " << parameters.timeStep
            << "\nAverage density by the end of simulation: " << parameters.avgDensity << " kg/m^3"
            << "\nExecution time total: " << duration.count() << "s" << std::endl;
        std::cout << "Execution time per step: " << duration.count() / steps << "s" << std::endl << std::endl;

        particles2D.clear();
		/*grid2D.updateGrid();*/
    }
}

int main() {
    // Read parameters from the input file
    parameters.readParameters();

    if (parameters.dimensions == 2){
        Initialization2D();
        Visualize2D();
        //SEvsPPE();
    }
    else if (parameters.dimensions == 3) {
        Initialization();
        Visualize();
        //VisualizeGhosts();
	}

    return 0;
}
