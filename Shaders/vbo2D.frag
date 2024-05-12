#version 330 core

in vec3 Color;

void main() {
	float distance = length(gl_PointCoord - vec2(0.5, 0.5));
    if (distance > 0.5) {
        discard;
    }
	gl_FragColor = vec4(Color, 1.0);
}