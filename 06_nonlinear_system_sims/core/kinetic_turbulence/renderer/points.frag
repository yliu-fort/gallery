#version 430 core
out vec4 FragColor;

//in vec2 TexCoords;
in float rr;

void main()
{

    //FragColor = vec4(rr,1.0f - rr,4.0f*rr*(1.0f-rr),1.0f);
    FragColor = vec4(0.0f,0.0f,0.0f,1.0f);
}
