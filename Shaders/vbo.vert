#version 330 core

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 color;
layout(location = 2) in float transparency;

out vec3 Color;
out float Transparency;

void main() {
	gl_Position = projection * view * model * vec4(position, 1.0);

	Color = color;
	Transparency = transparency;
}