#include "../headers/Camera.h"

Camera::Camera(int width, int height, int depth) {
    this->width = width;
    this->height = height;
    this->position = glm::vec3(width / 5, height / 5, -depth);
    this->orientation = glm::vec3(width, height, 0.0f);
}

void Camera::Matrix(float width, float height, float depth, unsigned int shader) {
    glm::mat4 projection = glm::perspective(glm::radians(25.0f), width / height, 0.1f, 5 * depth);
    projection = glm::translate(projection, glm::vec3(-width * 1.15, -height / 3, -4 * depth));
    GLuint projectionUniform = glGetUniformLocation(shader, "projection");
    glUniformMatrix4fv(projectionUniform, 1, GL_FALSE, glm::value_ptr(projection));

    glm::mat4 view = glm::lookAt(this->position, this->position + this->orientation, this->up);
    GLuint viewUniform = glGetUniformLocation(shader, "view");
    glUniformMatrix4fv(viewUniform, 1, GL_FALSE, glm::value_ptr(view));

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, -1.0f, 0.0f));
    model = glm::rotate(model, glm::radians(50.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    GLuint modelUniform = glGetUniformLocation(shader, "model");
    glUniformMatrix4fv(modelUniform, 1, GL_FALSE, glm::value_ptr(model));
}

void Camera::Inputs(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
        position += orientation * speed * 0.01f;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) {
        position -= orientation * speed * 0.01f;
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        position -= glm::normalize(glm::cross(orientation, up)) * speed * 5.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        position += glm::normalize(glm::cross(orientation, up)) * speed * 5.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        position += up * speed * 5.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        position -= up * speed * 5.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
        speed = 0.5f;
    }
    else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_RELEASE) {
        speed = 0.1f;
    }

    // Mouse motion
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

        if (firstClick) {
            glfwSetCursorPos(window, width / 2, height / 2);
            firstClick = false;
        }

        double mouseX, mouseY;
        glfwGetCursorPos(window, &mouseX, &mouseY);

        float rotationX = sensitivity * (float)(mouseX - width / 2) / width;
        float rotationY = sensitivity * (float)(mouseY - height / 2) / height;

        orientation = glm::rotate(orientation, glm::radians(-rotationX), glm::normalize(glm::cross(orientation, up)));
        orientation = glm::rotate(orientation, glm::radians(-rotationY), up);

        glfwSetCursorPos(window, width / 2, height / 2);
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        firstClick = true;
    }
}