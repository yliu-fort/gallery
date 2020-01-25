#version 330 core
out vec4 FragColor;

in GS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec3 Color;
    //vec2 TexCoords;
} fs_in;

void main()
{
    FragColor = vec4(fs_in.Normal,1.0);
    // Add visibility
    //FragColor.a = 1.0f;
}
