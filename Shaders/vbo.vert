#version 330 core

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 color;
layout(location = 2) in double transparency;

out vec3 Color;
out double Transparency;

void main() {
	gl_Position = projection * view * model * vec4(position, 1.0);

	Color = color;
	Transparency = transparency;
}