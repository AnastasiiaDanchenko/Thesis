#include "../headers/Visualization.h"

const int IMGUI_WINDOW_WIDTH = 300;

bool isSimulationRunning = true;

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

GLuint vao = 0;
GLuint vbo = 0;

#define VERTICES_CAPACITY 1000000
Vertex vertices[VERTICES_CAPACITY];
size_t verticesCount = 0;

// Close the window when pressing ESC
void keyCallbackVBO(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
}

float magnitude(Eigen::Vector3f vec) {
    return sqrt(vec.x() * vec.x() + vec.y() * vec.y() + vec.z() * vec.z());
}

float mapColor(float value, float inMin, float inMax, float outMin, float outMax) {
    return outMin + (outMax - outMin) * ((value - inMin) / (inMax - inMin));
}

void HSVtoRGB(float* r, float* g, float* b, float h, float s, float v) {
	int i;
	float f, p, q, t;

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

void pushVertex(Eigen::Vector3f position, float r, float g, float b, float t) {
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

void Visualize() {
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
    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "SPH solver in 3D", nullptr, nullptr);
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

    // Shader inputs were here
    Camera camera(WINDOW_WIDTH, WINDOW_HEIGHT, SCENE_DEPTH);

    // Initialize the buffers
    initBuffers();

    // boundary box
    Eigen::Vector3f minBound = Eigen::Vector3f(SPACING * 3, SPACING * 3, SPACING * 3);
    Eigen::Vector3f maxBound = Eigen::Vector3f(WINDOW_WIDTH - SPACING * 3, WINDOW_HEIGHT - SPACING * 3, SCENE_DEPTH - SPACING * 3);

    // event loop
    while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0) {
        //Simulation();
        SimulationIISPH();

        clearBuffers();

        glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

        // Draw the particles
        for (auto p : particles) {
            if (p.isFluid) {
                //colorize particles based on their speed
                float speed = magnitude(p.velocity);
                float hue = mapColor(speed, 0.0f, 100.0f, 240.0f, 0.0f);
                float r, g, b;
                HSVtoRGB(&r, &g, &b, hue, 1.0f, 1.0f);

                pushVertex(p.position, r, g, b, 1.0f);
            }
            else {
                if ((p.position.x() == minBound.x() || p.position.x() == maxBound.x()) && 
                    (p.position.y() == minBound.y() || p.position.y() == maxBound.y()) &&
                    p.position.z() >= minBound.z() && p.position.z() <= maxBound.z() ||
                    (p.position.x() == minBound.x() || p.position.x() == maxBound.x()) &&
                    (p.position.z() == minBound.z() || p.position.z() == maxBound.z()) &&
                    p.position.y() >= minBound.y() && p.position.y() <= maxBound.y() || 
                    (p.position.z() == minBound.z() || p.position.z() == maxBound.z()) &&
                    (p.position.y() == minBound.y() || p.position.y() == maxBound.y()) &&
                    p.position.x() >= minBound.x() && p.position.x() <= maxBound.x()) {
					pushVertex(p.position, 1.0f, 1.0f, 1.0f, 1.0f);
				}
			}
		}

        syncBuffers();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        camera.Inputs(window);
        camera.Matrix(WINDOW_WIDTH, WINDOW_HEIGHT, SCENE_DEPTH, shader);

        glPointSize(15.0f); // Set the point size
        glEnable(GL_DEPTH_TEST); // Enable depth test
        glDrawArraysInstanced(GL_POINTS, 0, 1, verticesCount);
        glDisable(GL_DEPTH_TEST); // Disable depth test

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteProgram(shader);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);

    glfwTerminate();
}
