#version 430 core
//layout (location = 0) in vec3 aPos;
//layout (location = 1) in vec2 aTexCoords;

//out vec2 TexCoords;
out float rr;

// Input
layout(std140, binding = 0) buffer readonlybuffer
{
    vec4 lhs[];
};

uniform mat4 projectionMatrix;
uniform int numparticles;

void main()
{
    //TexCoords = aTexCoords;

    vec3 aPos = lhs[gl_InstanceID].xyz;
    //vec3 aPos = vec3(gl_InstanceID%256,gl_InstanceID/256,0);
    rr =  float(gl_InstanceID)/numparticles;

    gl_Position = projectionMatrix*vec4(aPos, 1.0);
    gl_PointSize = 2.0f;
}
