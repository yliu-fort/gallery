#version 330 core
layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

//out vec3 aNorm;
out GS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec3 TexCoords;
} gs_out;

in VS_OUT {
    ivec3 gridPos;
    vec3 normal;
} gs_in[];

uniform mat4 projectionMatrix;

// Marching cube
uniform sampler3D volumeTex;
//uniform isampler1D triTex;
//uniform isampler1D numVertsTex;
//uniform float isoValue;
uniform float depth;
uniform ivec3 gridSize;
uniform vec3 voxelSize;

// sample volume data set at a point
float sampleVolume(ivec3 p, ivec3 gridSize)
{
    p.x = min(p.x, gridSize.x - 1);
    p.y = min(p.y, gridSize.y - 1);
    p.z = min(p.z, gridSize.z - 1);
    return texelFetch(volumeTex, p, 0).a;
}

void render_quad(vec3 p);
void render_point(vec3 color, vec4 position);
void main() {
    render_quad(gl_in[0].gl_Position.xyz);
}

// calculate triangle normal
vec3 calcNormal(vec3 v0, vec3 v1, vec3 v2)
{
    vec3 edge0 = v1 - v0;
    vec3 edge1 = v2 - v0;
    // note - it's faster to perform normalization in vertex shader rather than here
    return cross(edge0, edge1);
}

void render_quad(vec3 p)
{
    //aNorm = gs_in[0].gridPos*voxelSize; // gs_in[0] since there's only one input vertex
    vec3 v[4];
    v[0] = p + vec3(          0,           0, depth);
    v[1] = p + vec3(voxelSize.x,           0, depth);
    v[2] = p + vec3(          0, voxelSize.y, depth);
    v[3] = p + vec3(voxelSize.x, voxelSize.y, depth);

    // calculate triangle surface normal
    gs_out.Normal = gs_in[0].normal;
    //vec3 n = normalize(calcNormal(v[0], v[1], v[2]));

    // 1:bottom-left
    gs_out.FragPos = v[0]/(voxelSize/voxelSize.y);;
    //gs_out.Normal = n;
    gs_out.TexCoords = v[0] + 0.5f*vec3(voxelSize.xy, 0.0f);
    gl_Position = projectionMatrix*vec4(v[0], 1.0);
    EmitVertex();

    // 2:bottom-right
    gs_out.FragPos = v[1]/(voxelSize/voxelSize.y);;
    //gs_out.Normal = n;
    gs_out.TexCoords = v[1] + 0.5f*vec3(voxelSize.xy, 0.0f);
    gl_Position = projectionMatrix*vec4(v[1], 1.0);
    EmitVertex();

    // 3:top-left
    gs_out.FragPos = v[2]/(voxelSize/voxelSize.y);;
    //gs_out.Normal = n;
    gs_out.TexCoords = v[2] + 0.5f*vec3(voxelSize.xy, 0.0f);
    gl_Position = projectionMatrix*vec4(v[2], 1.0);
    EmitVertex();

    // 4:top-right
    gs_out.FragPos = v[3]/(voxelSize/voxelSize.y);;
    //gs_out.Normal = n;
    gs_out.TexCoords = v[3] + 0.5f*vec3(voxelSize.xy, 0.0f);
    gl_Position = projectionMatrix*vec4(v[3], 1.0);
    EmitVertex();

    // For wireframe
    //gl_Position = position + vec4(voxelSize*vec3(0, 0, 0.0),0.0); // 1:bottom-left
    //gl_Position = projectionMatrix*gl_Position;
    //EmitVertex();
    EndPrimitive();
}

void render_point(vec3 color, vec4 position)
{
    //aNorm = color; // gs_in[0] since there's only one input vertex
    gl_Position = projectionMatrix*position; // 1:bottom-left
    EmitVertex();
    EndPrimitive();
}
