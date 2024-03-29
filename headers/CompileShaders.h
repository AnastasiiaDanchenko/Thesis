#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <stdlib.h>
#include <string.h>

struct ShaderProgramSource {
    std::string VertexSource;
    std::string FragmentSource;
};

static ShaderProgramSource ParseShader(const std::string& vertexShaderPath, const std::string& fragmentShaderPath) {
    ShaderProgramSource source;
    std::ifstream vertexFile(vertexShaderPath);
    std::ifstream fragmentFile(fragmentShaderPath);

    std::stringstream vertexStream, fragmentStream;

    if (!vertexFile.is_open() || !fragmentFile.is_open()) {
        std::cerr << "Failed to open shader files!" << std::endl;
        return source;
    }

    vertexStream << vertexFile.rdbuf();
    fragmentStream << fragmentFile.rdbuf();

    source.VertexSource = vertexStream.str();
    source.FragmentSource = fragmentStream.str();

    return source;
}

static unsigned int CompileShader(unsigned int type, const std::string& source) {
    unsigned int id = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(id, 1, &src, nullptr);
    glCompileShader(id);

    int result;
    glGetShaderiv(id, GL_COMPILE_STATUS, &result);

    if (result == GL_FALSE) {
        int length;
        glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
        char* message = (char*)alloca(length * sizeof(char));
        glGetShaderInfoLog(id, length, &length, message);
        std::cerr << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << " shader!" << std::endl;
        std::cerr << message << std::endl;
        glDeleteShader(id);
        return 0;
    }

    return id;
}

// Function to compile and link a shader program
static unsigned int CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
    unsigned int program = glCreateProgram();
    unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
    unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);

    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glValidateProgram(program);

    glDeleteShader(vs);
    glDeleteShader(fs);

    return program;
}
