#include "../headers/Visualization.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../headers/stb_image_write.h"
#include <filesystem>

const int IMGUI_WINDOW_WIDTH = 300;

bool isSimulationRunning = false;
int simulationType = 0;

typedef enum {
    POSITION_ATTRIBUTE = 0,
    COLOR_ATTRIBUTE,
    TRANCPARENCY_ATTRIBUTE,
    COUNT_ATTRIBUTES
} Attributes;

typedef struct {
    float x, y, z;
    float r, g, b;
    float t;
} Vertex;

typedef enum {
    POSITION_ATTRIBUTE2D = 0,
    COLOR_ATTRIBUTE2D,
    COUNT_ATTRIBUTES2D
} Attributes2D;

typedef struct {
    float x, y;
    float r, g, b;
} Vertex2D;

GLuint vao = 0;
GLuint vbo = 0;

#define VERTICES_CAPACITY 1000000
Vertex vertices[VERTICES_CAPACITY];
Vertex2D vertices2D[VERTICES_CAPACITY];
size_t verticesCount = 0;

// Close the window when pressing ESC
void keyCallbackVBO(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
}

double magnitude(Eigen::Vector3d vec) {
    return sqrt(vec.x() * vec.x() + vec.y() * vec.y() + vec.z() * vec.z());
}

double magnitude2D(Eigen::Vector2d vec) {
    return sqrt(vec.x() * vec.x() + vec.y() * vec.y());
}

double mapColor(double value, double inMin, double inMax, double outMin, double outMax) {
    return outMin + (outMax - outMin) * ((value - inMin) / (inMax - inMin));
}

void HSVtoRGB(double* r, double* g, double* b, double h, double s, double v) {
	int i;
	double f, p, q, t;

    if (s == 0) {
		*r = *g = *b = v; // Achromatic color (gray)
		return;
	}

	h /= 60;
	i = floor(h);
	f = h - i;
	p = v * (1 - s);
	q = v * (1 - s * f);
	t = v * (1 - s * (1 - f));

    switch (i) {
	case 0:
		*r = v;
		*g = t;
		*b = p;
		break;
	case 1:
		*r = q;
		*g = v;
		*b = p;
		break;
	case 2:
		*r = p;
		*g = v;
		*b = t;
		break;
	case 3:
		*r = p;
		*g = q;
		*b = v;
		break;
	case 4:
		*r = t;
		*g = p;
		*b = v;
		break;
	default:
		*r = v;
		*g = p;
		*b = q;
		break;
	}
}

void initBuffers() {
    // Create the vertex array object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create the vertex buffer object
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);

    // Initialize the vertex attributes
    glEnableVertexAttribArray(POSITION_ATTRIBUTE);
    glVertexAttribPointer(POSITION_ATTRIBUTE, 3, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (void*)offsetof(Vertex, x));
    glVertexAttribDivisor(POSITION_ATTRIBUTE, 1);

    glEnableVertexAttribArray(COLOR_ATTRIBUTE);
    glVertexAttribPointer(COLOR_ATTRIBUTE, 3, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (void*)offsetof(Vertex, r));
    glVertexAttribDivisor(COLOR_ATTRIBUTE, 1);

    glEnableVertexAttribArray(TRANCPARENCY_ATTRIBUTE);
    glVertexAttribPointer(TRANCPARENCY_ATTRIBUTE, 1, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (void*)offsetof(Vertex, t));
    glVertexAttribDivisor(TRANCPARENCY_ATTRIBUTE, 1);
}

void syncBuffers() {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices[0]) * verticesCount, vertices);
}

void pushVertex(Eigen::Vector3d position, double r, double g, double b, float t) {
    if (verticesCount >= VERTICES_CAPACITY) {
        std::cerr << "Vertex buffer overflow!" << std::endl;
        exit(EXIT_FAILURE);
    }

    vertices[verticesCount].x = position.x();
    vertices[verticesCount].y = position.y();
    vertices[verticesCount].z = position.z();
    vertices[verticesCount].r = r;
    vertices[verticesCount].g = g;
    vertices[verticesCount].b = b;
    vertices[verticesCount].t = t;

    verticesCount++;
}

void clearBuffers() {
    verticesCount = 0;
}

void Visualize(Solver& solver) {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set the OpenGL version to 3.3 to use modern shaders
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(
        IMGUI_WINDOW_WIDTH + parameters.windowSize.width, parameters.windowSize.height, 
        "SPH solver in 3D", nullptr, nullptr
    );

    if (!window) {
        std::cerr << "Failed to create GLFW window!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set the background color
    glClearColor(.05f, .05f, .05f, 1.0f);

    // Copile and link the shaders
    ShaderProgramSource source = ParseShader("Shaders/vbo.vert", "Shaders/vbo.frag");
    unsigned int shader = CreateShader(source.VertexSource, source.FragmentSource);

    glUseProgram(shader);

    GLint lightPosLoc = glGetUniformLocation(shader, "lightPosition");
    GLint lightColorLoc = glGetUniformLocation(shader, "lightColor");

    glUniform3f(lightPosLoc, 0.0f, parameters.windowSize.height, parameters.windowSize.depth);
    glUniform3f(lightColorLoc, 1.0f, 1.0f, 1.0f);

    // Shader inputs were here
    Camera camera(parameters.windowSize.width, parameters.windowSize.height, parameters.windowSize.depth);

    // Initialize the buffers
    initBuffers();

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // boundary box
    Eigen::Vector3d minBound = Eigen::Vector3d(parameters.spacing, parameters.spacing, parameters.spacing);
    Eigen::Vector3d maxBound = Eigen::Vector3d(
        parameters.windowSize.width     - parameters.spacing / 2, 
        parameters.windowSize.height    - parameters.spacing / 2,
        parameters.windowSize.depth / 2 - parameters.spacing / 2
    );

    Grid grid(parameters.spacing * 2);

    // event loop
    while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0) {
        if (isSimulationRunning) {
            SimulationIISPH(solver, parameters.simulationType);
        }

        clearBuffers();

        // Render ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glViewport(IMGUI_WINDOW_WIDTH, 0, parameters.windowSize.width, parameters.windowSize.height);

        for (auto p : particles) {
            if (p.position.z() <= parameters.slicingPlane) {
                if (p.isFluid) {
                    double speed = magnitude(p.velocity);
                    double hue = mapColor(speed, 0.0f, 100.0f, 240.0f, 0.0f);
                    double r, g, b;
                    HSVtoRGB(&r, &g, &b, hue, 1.0f, 1.0f);

                    pushVertex(p.position, r, g, b, 1.0f);
                }
                else {
                    double hue = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 0.0, 30.0);
                    double saturation = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 0.0,
                        1.0);
                    double value = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 1.0,
                        0.6);
                    double r, g, b;
                    HSVtoRGB(&r, &g, &b, hue, saturation, value);

                    pushVertex(p.position, r, g, b, 0.3f);
                }
            }
		}

        if (parameters.simulationType != 0) {
            for (auto& body : solver.getRigidBodies()) {
                pushVertex(body.getPositionCM(), 0.0f, 1.0f, 0.0f, 1.0f);
                for (auto& p : body.getOuterParticles()) {
					if (p.position.z() <= parameters.slicingPlane) {
						pushVertex(p.position, 1.0f, 0.0f, 0.0f, 1.0f);

                        /*double hue = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.rigidBody.density, 
                            0.0, 30.0);
                        double saturation = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * 
                            parameters.rigidBody.density, 0.0, 1.0);
                        double value = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.rigidBody.density,
                            1.0, 0.6);
                        double r, g, b;
                        HSVtoRGB(&r, &g, &b, hue, saturation, value);

                        pushVertex(p.position, r, g, b, 1.0f);*/
					}
                }
            }
        }

        syncBuffers();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        camera.Inputs(window);
        camera.Matrix(parameters.windowSize.width, parameters.windowSize.height, parameters.windowSize.depth, shader);

        glPointSize(parameters.spacing / 1.25);
        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_POINTS, 0, 1, verticesCount);
        glDisable(GL_DEPTH_TEST);

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height - 100));

        ImGui::Begin("Simulation Parameters");
        ImGui::Text("Number of fluid particles: %d", 
            parameters.particlesPerDimension.x * 
            parameters.particlesPerDimension.y * 
            parameters.particlesPerDimension.z
        );

        float tempTimeStep = static_cast<float>(parameters.timeStep);
        ImGui::Text("\Time step:");
        if (ImGui::SliderFloat("##TimeStep", &tempTimeStep, 0.0001, parameters.maxTimeStep)) {
            parameters.timeStep = static_cast<double>(tempTimeStep);
        }

        float tempGamma = static_cast<float>(parameters.gamma);
        ImGui::Text("\Gamma:");
        if (ImGui::SliderFloat("##Gamma", &tempGamma, 0.01, 1.0)) {
            parameters.gamma = static_cast<double>(tempGamma);
        }

        float tempOmega = static_cast<float>(parameters.omega);
        ImGui::Text("\Omega:");
        if (ImGui::SliderFloat("##Omega", &tempOmega, 0.01, 1.0)) {
            parameters.omega = static_cast<double>(tempOmega);
        }

        ImGui::Text("\Rest Density: %f", parameters.restDensity);
        ImGui::Text("\Average Density: %f", parameters.avgDensity);
        ImGui::Text("\Density Error, %%: %f", parameters.densityErr / parameters.restDensity * 100);
        ImGui::Text("\Density Error before PPE, %%: %f", parameters.firstErr / parameters.restDensity * 100);
        ImGui::Text("\Number of iterations (l): %d", parameters.nbIterations);
        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(0, parameters.windowSize.height - 150));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height));

        ImGui::Begin("Simulation Controls");
        if (ImGui::Button("Start/Pause")) { isSimulationRunning = !isSimulationRunning; }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            particles.clear();
            Initialization(solver, parameters.simulationType);
        }
        if (ImGui::Button("Move forward one time step")) { SimulationIISPH(solver, parameters.simulationType); }
        if (ImGui::Button("Start moving boundary")) { solver.initMovingBoundary(); }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteProgram(shader);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);

    glfwTerminate();
}

void initBuffers2D() {
    // Create the vertex array object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create the vertex buffer object
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices2D), vertices2D, GL_DYNAMIC_DRAW);

    // Initialize the vertex attributes
    glEnableVertexAttribArray(POSITION_ATTRIBUTE2D);
    glVertexAttribPointer(POSITION_ATTRIBUTE2D, 2, GL_FLOAT, GL_FALSE, sizeof(vertices2D[0]), (void*)offsetof(Vertex2D, x));
    glVertexAttribDivisor(POSITION_ATTRIBUTE2D, 1);

    glEnableVertexAttribArray(COLOR_ATTRIBUTE2D);
    glVertexAttribPointer(COLOR_ATTRIBUTE2D, 3, GL_FLOAT, GL_FALSE, sizeof(vertices2D[0]), (void*)offsetof(Vertex2D, r));
    glVertexAttribDivisor(COLOR_ATTRIBUTE2D, 1);
}

void syncBuffers2D() {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices2D[0]) * verticesCount, vertices2D);
}

// 2D visualization
void pushVertex2D(float x, float y, float r, float g, float b) {
    if (verticesCount >= VERTICES_CAPACITY) {
        std::cerr << "Vertex buffer overflow!" << std::endl;
        exit(EXIT_FAILURE);
    }

    vertices2D[verticesCount].x = x;
    vertices2D[verticesCount].y = y;
    vertices2D[verticesCount].r = r;
    vertices2D[verticesCount].g = g;
    vertices2D[verticesCount].b = b;

    verticesCount++;
}

void saveImage(const char* filepath, GLFWwindow* w) {
    int width, height;
    glfwGetFramebufferSize(w, &width, &height);
    GLsizei nrChannels = 3;
    GLsizei stride = nrChannels * width;
    stride += (stride % 4) ? (4 - stride % 4) : 0;
    GLsizei bufferSize = stride * height;
    std::vector<char> buffer(bufferSize);
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath, width, height, nrChannels, buffer.data(), stride);
}

void Visualize2D(Solver2D& solver) {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set the OpenGL version to 3.3 to use modern shaders
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(
        IMGUI_WINDOW_WIDTH + parameters.windowSize.width, 
        parameters.windowSize.height, "SPH solver in 2D", nullptr, nullptr
    );
    if (!window) {
        std::cerr << "Failed to create GLFW window!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set the background color
    glClearColor(.05f, .05f, .05f, 1.0f);

    // Copile and link the shaders
    ShaderProgramSource source = ParseShader("Shaders/vbo2D.vert", "Shaders/vbo2D.frag");
    unsigned int shader = CreateShader(source.VertexSource, source.FragmentSource);

    glUseProgram(shader);
    GLuint resolutionUniform = glGetUniformLocation(shader, "resolution");
    glUniform2f(resolutionUniform, parameters.windowSize.width, parameters.windowSize.height);

    // Initialize the buffers
    initBuffers2D();

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    //create save directory
    /*int count = 0;
    std::string path = "output/" + std::to_string(parameters.particlesPerDimension.x * parameters.particlesPerDimension.y) + "_" + std::to_string(ERR_THRESHOLD);
    std::filesystem::create_directories(path);*/

    Grid2D grid2D(parameters.spacing * 2);

    // event loop
    while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0) {
        if (isSimulationRunning) {
            if (simulationType == 2) RotatingBoundaryIISPH2D(solver);
			else SimulationIISPH2D(solver);
        }

        clearBuffers();

        // Render ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glViewport(IMGUI_WINDOW_WIDTH, 0, parameters.windowSize.width, parameters.windowSize.height);

        for (auto& p : particles2D) {
            if (p.isFluid) {
                /*double speed = magnitude2D(p.velocity);
                double hue = mapColor(speed, 0.0f, 100.0f, 240.0f, 0.0f);
                double r, g, b;
                HSVtoRGB(&r, &g, &b, hue, 1.0f, 1.0f);

                pushVertex2D(p.position.x(), p.position.y(), r, g, b);*/

                bool isNeighbor = false;
                for (auto& n : p.neighbors) {
                    if (n == &particles2D[parameters.visualizeNeighbors]) { isNeighbor = true; break; }
                }
                if (p.ID == parameters.visualizeNeighbors) { pushVertex2D(p.position.x(), p.position.y(), 1.0f, 1.0f, 0.0f); }
                else if (isNeighbor) { pushVertex2D(p.position.x(), p.position.y(), 0.0f, 1.0f, 0.0f); }
                else { pushVertex2D(p.position.x(), p.position.y(), 0.2f, 0.5f, 1.0f); }
            }
            else {
                double hue = mapColor(p.mass, 0.0, parameters.spacing * parameters.spacing * 
                    parameters.restDensity, 0.0, 30.0);
                double saturation = mapColor(p.mass, 0.0, parameters.spacing * parameters.spacing * 
                    parameters.restDensity, 0.0, 1.0);
                double value = mapColor(p.mass, 0.0, parameters.spacing * parameters.spacing * 
                    parameters.restDensity, 1.0, 0.6);
                double r, g, b;
                HSVtoRGB(&r, &g, &b, hue, saturation, value);

                pushVertex2D(p.position.x(), p.position.y(), r, g, b);
            }
        }

        syncBuffers2D();

        glClear(GL_COLOR_BUFFER_BIT);

        glPointSize(parameters.spacing * 2); // Set the point size
        glDrawArraysInstanced(GL_POINTS, 0, 1, verticesCount);

        //set fixed position for imgui window
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height - 100));

        ImGui::Begin("Simulation Parameters");
        ImGui::Text("Number of fluid particles: %d", 
            parameters.particlesPerDimension.x * parameters.particlesPerDimension.y
        );

        float tempTimeStep = static_cast<float>(parameters.timeStep);
        ImGui::Text("\Time step:");
        if (ImGui::SliderFloat("##TimeStep", &tempTimeStep, 0.0001, parameters.maxTimeStep)) {
            parameters.timeStep = static_cast<double>(tempTimeStep);
        }

        float tempGamma = static_cast<float>(parameters.gamma);
        ImGui::Text("\Gamma:");
        if (ImGui::SliderFloat("##Gamma", &tempGamma, 0.01, 1.0)) {
            parameters.gamma = static_cast<double>(tempGamma);
        }

        float tempOmega = static_cast<float>(parameters.omega);
        ImGui::Text("\Omega:");
        if (ImGui::SliderFloat("##Omega", &tempOmega, 0.01, 1.0)) {
            parameters.omega = static_cast<double>(tempOmega);
        }

        ImGui::Text("\Rest Density: %f", parameters.restDensity);
        ImGui::Text("\Average Density: %f", parameters.avgDensity);
        ImGui::Text("\Density Error, %%: %f", parameters.densityErr / parameters.restDensity * 100);
        ImGui::Text("\Density Error before PPE, %%: %f", parameters.firstErr / parameters.restDensity * 100);
        ImGui::Text("\Number of iterations (l): %d", parameters.nbIterations);
        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(0, parameters.windowSize.height - 150));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height));

        ImGui::Begin("Simulation Controls");
        if (ImGui::Button("Start/Pause")) { isSimulationRunning = !isSimulationRunning; }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) { 
            particles2D.clear(); 
            if (simulationType == 0) Initialization2D(solver);
            else if (simulationType == 1) MovingBoundaryInitialization(solver);
            else if (simulationType == 2) RotatingBoundaryInitialization(solver);
        }
        if (ImGui::Button("Surface Tension: ON/OFF")) {
            parameters.surfaceTension = !parameters.surfaceTension;
            particles2D.clear(); Initialization2D(solver);
        }
        if (ImGui::Button("Move forward one time step")) { 
            if (simulationType == 2) RotatingBoundaryIISPH2D(solver);
			else SimulationIISPH2D(solver); 
        }
        if (ImGui::Button("Change simulation type")) {
            simulationType = (simulationType + 1) % 3;
			particles2D.clear();
            if (simulationType == 0) {
                parameters.particlesPerDimension.x *= 2;
                Initialization2D(solver);
            }
			else if (simulationType == 1) MovingBoundaryInitialization(solver);
            else if (simulationType == 2) {
                parameters.particlesPerDimension.x /= 2;
                RotatingBoundaryInitialization(solver);
            }
        }
        if (simulationType == 1) {
            if (ImGui::Button("Start moving boundary")) { solver.moveBoundary(); }
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();

        /*std::stringstream filenameStream;
        filenameStream << path << "/" << count << ".png";
        std::string filename = filenameStream.str();
        saveImage(filename.c_str(), window);

        count++;*/
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteProgram(shader);
    glfwTerminate();
}

void VisualizeGhosts(Solver& solver) {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set the OpenGL version to 3.3 to use modern shaders
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(IMGUI_WINDOW_WIDTH + parameters.windowSize.width, parameters.windowSize.height, "SPH solver in 3D : projection onto plane", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set the background color
    glClearColor(.05f, .05f, .05f, 1.0f);

    // Copile and link the shaders
    ShaderProgramSource source = ParseShader("Shaders/vbo2D.vert", "Shaders/vbo2D.frag");
    unsigned int shader = CreateShader(source.VertexSource, source.FragmentSource);

    glUseProgram(shader);
    GLuint resolutionUniform = glGetUniformLocation(shader, "resolution");
    glUniform2f(resolutionUniform, parameters.windowSize.width, parameters.windowSize.height);

    // Initialize the buffers
    initBuffers2D();

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // event loop
    while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0) {
        if (isSimulationRunning) {
            SimulationIISPH(solver, parameters.simulationType);
        }

        clearBuffers();

        // Render ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glViewport(IMGUI_WINDOW_WIDTH, 0, parameters.windowSize.width, parameters.windowSize.height);

        for (auto& p : ghostParticles) {
            if (p.isFluid) {
                double speed = sqrt(p.velocity.x() * p.velocity.x() + p.velocity.y() * p.velocity.y());
                double hue = mapColor(speed, 0.0f, 100.0f, 240.0f, 0.0f);
                double r, g, b;
                HSVtoRGB(&r, &g, &b, hue, 1.0f, 1.0f);

                pushVertex2D(p.position.x() / 2, p.position.y() / 2, r, g, b);
            }
            else {
                double hue = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 0.0, 30.0);
                double saturation = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 0.0, 
                    1.0);
                double value = mapColor(p.mass, 0.0, pow(parameters.spacing, 3) * parameters.restDensity, 1.0, 0.6);
                double r, g, b;
                HSVtoRGB(&r, &g, &b, hue, saturation, value);

                pushVertex2D(p.position.x() / 2, p.position.y() / 2, r, g, b);
            }
        }

        syncBuffers2D();

        glClear(GL_COLOR_BUFFER_BIT);

        glPointSize(parameters.spacing); // Set the point size
        glDrawArraysInstanced(GL_POINTS, 0, 1, verticesCount);

        //set fixed position for imgui window
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height - 100));

        ImGui::Begin("Simulation Parameters");
        ImGui::Text("Number of fluid particles: %d", parameters.particlesPerDimension.x * 
            parameters.particlesPerDimension.y * parameters.particlesPerDimension.z
        );

        float tempTimeStep = static_cast<float>(parameters.timeStep);
        ImGui::Text("\Time step:");
        if (ImGui::SliderFloat("##TimeStep", &tempTimeStep, 0.0001, parameters.maxTimeStep)) {
            parameters.timeStep = static_cast<double>(tempTimeStep);
        }

        float tempGamma = static_cast<float>(parameters.gamma);
        ImGui::Text("\Gamma:");
        if (ImGui::SliderFloat("##Gamma", &tempGamma, 0.01, 1.0)) {
            parameters.gamma = static_cast<double>(tempGamma);
        }

        float tempOmega = static_cast<float>(parameters.omega);
        ImGui::Text("\Omega:");
        if (ImGui::SliderFloat("##Omega", &tempOmega, 0.01, 1.0)) {
            parameters.omega = static_cast<double>(tempOmega);
        }

        ImGui::Text("\Rest Density: %f", parameters.restDensity);
        ImGui::Text("\Average Density: %f", parameters.avgDensity);
        ImGui::Text("\Density Error, %%: %f", parameters.densityErr / parameters.restDensity * 100);
        ImGui::Text("\Density Error before PPE, %%: %f", parameters.firstErr / parameters.restDensity * 100);
        ImGui::Text("\Number of iterations (l): %d", parameters.nbIterations);
        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(0, parameters.windowSize.height - 150));
        ImGui::SetNextWindowSize(ImVec2(IMGUI_WINDOW_WIDTH, parameters.windowSize.height));

        ImGui::Begin("Simulation Controls");
        if (ImGui::Button("Start/Pause")) { isSimulationRunning = !isSimulationRunning; }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            particles.clear();
			ghostParticles.clear();
            Initialization(solver, parameters.simulationType);
        }
        if (ImGui::Button("Move forward one time step")) {
            SimulationIISPH(solver, parameters.simulationType);
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteProgram(shader);
    glfwTerminate();
}

void ExportPLY(Solver& solver) {
    int nbFrames = 100;
    Grid grid(parameters.spacing * 2);

    std::cout << "Exporting simulation frames to .ply files..." << std::endl;

    int nbRigidParticles = 0;
    for (auto& body : solver.getRigidBodies()) {
        nbRigidParticles += body.getOuterParticles().size();
	}

    for (int i = 0; i < nbFrames; i++) {
		SimulationIISPH(solver, parameters.simulationType);

        double scaleFactor = 0.01;
		
        std::string fileName = "output/blender_frames/" + std::to_string(i) + ".ply";
        std::ofstream file(fileName);

        if (!file.is_open()) {
			std::cerr << "Failed to open .ply file for writing." << std::endl;
			exit(EXIT_FAILURE);
		}

        file << "ply" << std::endl;
		file << "format ascii 1.0" << std::endl;
		file << "element vertex " << //particles.size() + nbrigidParticles << std::endl;
            parameters.particlesPerDimension.x * 
            parameters.particlesPerDimension.y * 
			parameters.particlesPerDimension.z + nbRigidParticles << std::endl;
		file << "property float x" << std::endl;
		file << "property float y" << std::endl;
		file << "property float z" << std::endl;
		file << "end_header" << std::endl;

		for (auto& p : particles) {
            if (p.position.z() < 0) {
                std::cout << "Negative z-coordinate, " << p.ID << std::endl;
            }

            if (p.isFluid) {
                file << p.position.x() * scaleFactor << " " << 
                        p.position.z() * scaleFactor << " " << 
                        p.position.y() * scaleFactor << std::endl;
            }
		}

        for (auto& body : solver.getRigidBodies()) {
			for (auto& p : body.getOuterParticles()) {
				file << p.position.x() * scaleFactor << " " << 
						p.position.z() * scaleFactor << " " << 
						p.position.y() * scaleFactor << std::endl;
			}
		}

		file.close();

        int barWidth = 50;
        float progress = static_cast<float>(i + 1) / nbFrames;
        std::cout << "[";
        int pos = static_cast<int>(barWidth * progress);

        for (int j = 0; j < barWidth; ++j) {
            if (j < pos) std::cout << "=";
            else if (j == pos) std::cout << ">";
            else std::cout << " ";
        }

        std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100.0) << " %\r";
        std::cout.flush();
	}

    std::cout << "\nExported " << nbFrames << " frames to .ply files." << std::endl;
}
