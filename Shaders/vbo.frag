#version 330 core

in vec3 Color;
in float Transparency;

uniform vec3 viewPos;
uniform vec3 lightPosition;
uniform vec3 lightColor;

void main() {
    vec2 coord = gl_PointCoord - vec2(0.5, 0.5);
    float distance = length(coord);

    if (distance > 0.5) {
        discard; // Discard the fragments outside the circle
    }

    vec3 normal = vec3(coord, sqrt(0.25 - distance * distance) * 2.0);
    normal = normalize(normal);

    // Lighting calculations
    vec3 lightDir = normalize(lightPosition - vec3(gl_FragCoord.xy, 0.0));
    float diff = max(dot(normal, lightDir), 0.0);

    vec3 viewDir = normalize(viewPos - vec3(gl_FragCoord.xy, 0.0));
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);

    vec3 ambient = 0.1 * lightColor;
    vec3 diffuse = diff * lightColor;
    vec3 specular = spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * Color;

	gl_FragColor = vec4(result, Transparency);
}