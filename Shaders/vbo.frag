#version 330 core

in vec3 Color;
in float Transparency;

void main() {
	float distance = length(gl_PointCoord - vec2(0.5, 0.5));
    if (distance > 0.5) {
        discard; // Discard the fragments outside the circle
    }

	gl_FragColor = vec4(Color, Transparency);
}