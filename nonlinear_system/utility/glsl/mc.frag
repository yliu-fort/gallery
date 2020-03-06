#version 330 core
out vec4 FragColor;

in GS_OUT {
    vec3 FragPos;
    vec3 Normal;
    //vec3 Color;
    vec3 TexCoords;
} fs_in;

uniform sampler3D volumeTex;
//uniform sampler2D floorTexture;
uniform vec3 viewPos;
vec3 lightPos;
bool blinn;
bool trueColor;

void main()
{
    // Will switch to uniform and integrate to gui interface
    trueColor = false;
    blinn = true;
    lightPos = vec3(0.5f,0.5f,3.0f);

    // fetch data
    vec3 color = texture(volumeTex, fs_in.TexCoords).rgb;

    // if rgb color
    if(trueColor)
    {
        FragColor = vec4(color, 1.0f);
        return;
    }
    // ambient
    vec3 ambient = 0.55 * color;
    // diffuse
    vec3 lightDir = normalize(lightPos - fs_in.FragPos); // directional light
    vec3 normal = normalize(fs_in.Normal);
    float diff = max(dot(lightDir, normal), 0.0);
    vec3 diffuse = 0.3 * diff * color;
    // specular
    vec3 viewDir = normalize(viewPos - fs_in.FragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = 0.0;
    if(blinn)
    {
        vec3 halfwayDir = normalize(lightDir + viewDir);
        spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);
    }
    else
    {
        vec3 reflectDir = reflect(-lightDir, normal);
        spec = pow(max(dot(viewDir, reflectDir), 0.0), 8.0);
    }
    vec3 specular = vec3(0.3) * spec; // assuming bright white light color
    FragColor = vec4(ambient + diffuse + specular, 1.0);
}
