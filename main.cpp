#include "headers/Visualization.h"

// Close the window when pressing ESC
void keyCallback_old(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
}

int main() {
    // Read parameters from the input file
    readParameters();

    if (DIMENSIONS == 2){
        Initialization2D();
        //Visualize2D();
        SaveToDisk2D();
    }
    else if (DIMENSIONS == 3) {
        Initialization();
        Visualize();
	}

    return 0;
}
