#version 330 core
out vec4 FragColor;

in vec3 aNorm;

void main()
{
    FragColor = vec4(aNorm,1.0);
    // Add visibility
    //FragColor.a = 1.0f;
}
