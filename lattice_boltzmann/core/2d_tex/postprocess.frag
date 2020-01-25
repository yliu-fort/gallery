#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D renderTarget;
//uniform sampler2D occlusion;
uniform int renderType;
uniform sampler1D colorMap;

layout (std140) uniform LatticeConstants
{
    mat3 wi;
    mat3 cx;
    mat3 cy;
    float csSqrInv;
    float csSqr;
    float Reg;
    float u_max;
    float tau;
    int NX,NY;
};

uniform float slidebar;

ivec2 relocUV(ivec2 v, int nx, int ny)
{
    return ivec2(v.x - ((v.x + NX)/NX - 1)*NX, v.y - ((v.y + NY)/NY - 1)*NY);
}

void main()
{
    vec4 val = vec4(texture(renderTarget, TexCoords));
    vec3 result = vec3(0.0);
    //float shadow = texture(occlusion, TexCoords).w;

    // Vorticity
    if (renderType == 0)
    {
        ivec2 cord[4];
        cord[0] = ivec2(TexCoords.x*NX, TexCoords.y*NY) - ivec2(cx[0][1], cy[0][1]);cord[0] = relocUV(cord[0],NX, NY);
        cord[1] = ivec2(TexCoords.x*NX, TexCoords.y*NY) - ivec2(cx[0][2], cy[0][2]);cord[1] = relocUV(cord[1],NX, NY);
        cord[2] = ivec2(TexCoords.x*NX, TexCoords.y*NY) - ivec2(cx[1][0], cy[1][0]);cord[2] = relocUV(cord[2],NX, NY);
        cord[3] = ivec2(TexCoords.x*NX, TexCoords.y*NY) - ivec2(cx[1][1], cy[1][1]);cord[3] = relocUV(cord[3],NX, NY);

        float v1 = texelFetch(renderTarget, cord[0],0).y;
        float u1 = texelFetch(renderTarget, cord[1],0).x;
        float v2 = texelFetch(renderTarget, cord[2],0).y;
        float u2 = texelFetch(renderTarget, cord[3],0).x;
        float vort = (((v1-v2)-(u1-u2))/slidebar+1)/2;

        result = texture(colorMap, vort).rgb;

        //result = (1-shadow)*result + shadow*vec3(0.2);
    }
    // Umag
    if (renderType == 1)
    {
        vec2 uv = (val.xy)/0.127;
        float index = sqrt(uv.x*uv.x + uv.y*uv.y);
        result = texture(colorMap, index).rgb;

        //result = (1-shadow)*result + shadow*vec3(0.2);
    }
    if (renderType == 2)
    {
        //float index = 1.0 - exp(-abs(val.z - 1.0) * 1200.0*slidebar);
        //result = 1.0 - vec3(index);
        //result = (1-shadow)*result + shadow*vec3(0.2);
        //result = vec3(min(val.z,1-shadow));
    }

    FragColor = vec4(result, 1.0);
}
