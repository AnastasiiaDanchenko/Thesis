#include "headers/Visualization.h"

// Close the window when pressing ESC
void keyCallback_old(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
}

void SEvsPPE(Solver2D& solver) {
    std::vector<std::string> METHODS = { "PPE" };
    int steps = 50;

    for (auto method : METHODS) {
        Initialization2D(solver);
        std::cout << "Method: " << method
            << "\nNumber of particles: " << parameters.particlesPerDimension.x * parameters.particlesPerDimension.y << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        if (method == "SE") {
            parameters.timeStep = 0.000005;
            for (int i = 0; i < steps; i++) {
                Simulation2D(solver);
            }
        }
        else if (method == "PPE") {
            parameters.timeStep = 0.05;
            for (int i = 0; i < steps; i++) {
                SimulationIISPH2D(solver);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        std::cout << "\nTime step: " << parameters.timeStep
            << "\nAverage density by the end of simulation: " << parameters.avgDensity << " kg/m^3"
            << "\nExecution time per step: " << duration.count() / steps << "s" << std::endl << std::endl;

        particles2D.clear();
    }
}

void AnalyzeIISPH3D(Solver& solver) {
    int steps = 50;
    std::cout << "Number of particles: " << parameters.particlesPerDimension.x * parameters.particlesPerDimension.y *
        parameters.particlesPerDimension.z << std::endl;

    auto startTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationPressure = std::chrono::duration<double>::zero();
    for (int i = 0; i < steps; i++) {
        solver.neighborSearch();
        solver.boundaryMassUpdate();
        solver.computeDensity();
        solver.predictVelocity();
        solver.computeSourceTerm();
        solver.computeDiagonalElement();

        auto startPressure = std::chrono::high_resolution_clock::now();
        solver.compressionConvergence();
        auto endPressure = std::chrono::high_resolution_clock::now();
        durationPressure += endPressure - startPressure;

        solver.updateParticles();
    }
    auto endTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationTotal = endTotal - startTotal;

    std::cout << "\nTime step: " << parameters.timeStep
        << "\nAverage density by the end of simulation: " << parameters.avgDensity << " kg/m^3"
        << "\nNumber of pressure solver iterations: " << parameters.iterationsCount
        << "\nExecution time total: " << durationTotal.count() << "s" << std::endl;
    std::cout << "Execution time for pressure convergence: " << durationPressure.count() / steps << "s" << std::endl << std::endl;
    std::cout << "Execution time per step: " << durationTotal.count() / steps << "s" << std::endl;

    particles.clear();
}

void AnalyzeRelaxation(Solver& solver) {
    std::cout << "Number of particles: " << parameters.particlesPerDimension.x * parameters.particlesPerDimension.y *
        parameters.particlesPerDimension.z << std::endl;

    while (parameters.avgDensity <= parameters.restDensity) {
        SimulationIISPH(solver);
    }

    int steps = 50;
    std::vector<float> omegas = { 0.05, 0.1, 0.2, 0.5, 0.7 };
    for (float omega : omegas) {
        parameters.omega = omega;
        auto startTotal = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> durationPressure = std::chrono::duration<double>::zero();
        for (int i = 0; i < steps; i++) {
			solver.neighborSearch();
			solver.boundaryMassUpdate();
			solver.computeDensity();
			solver.predictVelocity();
			solver.computeSourceTerm();
			solver.computeDiagonalElement();

			auto startPressure = std::chrono::high_resolution_clock::now();
			solver.compressionConvergence();
			auto endPressure = std::chrono::high_resolution_clock::now();
			durationPressure += endPressure - startPressure;

			solver.updateParticles();
		}
        auto endTotal = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> durationTotal = endTotal - startTotal;

        std::cout << "\nTime step: " << parameters.timeStep
            << "\nRelaxation factor: " << parameters.omega 
			<< "\nAverage density by the end of simulation: " << parameters.avgDensity << " kg/m^3"
			<< "\nNumber of pressure solver iterations: " << parameters.iterationsCount
			<< "\nExecution time total: " << durationTotal.count() << "s"
            << "\nExecution time for pressure convergence: " << durationPressure.count() / steps << "s"
            << "\nExecution time per step: " << durationTotal.count() / steps << "s" << std::endl << std::endl;
    }

}

void ExportSequencePLY(Solver& solver) {
    for (int exp = 0; exp < 30; exp++){
		ExportPLY(solver, exp);
		//Visualize(solver);

        const int totalSteps = 100;

        for (int i = 0; i < totalSteps; i++) {
            SimulationIISPH(solver);

            float progress = (i + 1) / static_cast<float>(totalSteps);
            int barWidth = 50;
            int pos = static_cast<int>(progress * barWidth);

            std::cout << "\r[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
            std::flush(std::cout);
        }
        std::cout << std::endl;
    }
}

int main() {
    // Read parameters from the input file
    parameters.readParameters();

    if (parameters.dimensions == 2){
		Solver2D solver;
        Initialization2D(solver);
        Visualize2D(solver);
    }
    else if (parameters.dimensions == 3) {
		Solver solver;
        Initialization(solver);
        //Visualize(solver);
        ExportSequencePLY(solver);
	}

    return 0;
}
