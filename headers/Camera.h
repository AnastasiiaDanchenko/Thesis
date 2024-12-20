#ifndef CAMERA_CLASS_H
#define CAMERA_CLASS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/vector_angle.hpp>

#include "CompileShaders.h"

class Camera {
public:
	glm::vec3 position;
	glm::vec3 orientation;
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

	bool firstClick = true;

	int width, height;

	float speed = 0.1f;
	float sensitivity = 50.0f;

	Camera(int width, int height, int depth);

	void Matrix(float width, float height, float depth, unsigned int shader);
	void Inputs(GLFWwindow* window);
};

#endif