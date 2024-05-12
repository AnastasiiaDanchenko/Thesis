#version 330 core

float PARTICLE_SIZE = 0.5f; // customized particle size in pixels

uniform vec2 resolution;

layout(location = 0) in vec2 position;
layout(location = 1) in vec3 color;

out vec3 Color;

vec2 screen_to_ndc(vec2 pos) {
	return (pos - resolution / 4.0) / (resolution / 4.0);
}

void main() {
	vec2 uv = vec2(
		float(gl_VertexID & 1),
		float((gl_VertexID >> 1) & 1)
	);
	gl_Position = vec4(screen_to_ndc(position + uv * PARTICLE_SIZE), 0.0, 1.0);

	Color = color;
}