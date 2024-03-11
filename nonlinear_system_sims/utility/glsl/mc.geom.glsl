#version 330 core
layout (points) in;
layout (triangle_strip, max_vertices = 15) out;

out GS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec3 TexCoords;
} gs_out;

in VS_OUT {
    ivec3 gridPos;
} gs_in[];

uniform mat4 projectionMatrix;

// Marching cube
uniform sampler3D volumeTex;
uniform isampler1D triTex;
uniform isampler1D numVertsTex;
uniform float isoValue;
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

// classify voxel based on number of vertices it will generate
// one thread per voxel
int classifyVoxel()
{
    ivec3 gridPos = gs_in[0].gridPos;

    // read field values at neighbouring grid vertices
    float field[8];
    field[0] = sampleVolume(gridPos, gridSize);
    field[1] = sampleVolume(gridPos + ivec3(1, 0, 0), gridSize);
    field[2] = sampleVolume(gridPos + ivec3(1, 1, 0), gridSize);
    field[3] = sampleVolume(gridPos + ivec3(0, 1, 0), gridSize);
    field[4] = sampleVolume(gridPos + ivec3(0, 0, 1), gridSize);
    field[5] = sampleVolume(gridPos + ivec3(1, 0, 1), gridSize);
    field[6] = sampleVolume(gridPos + ivec3(1, 1, 1), gridSize);
    field[7] = sampleVolume(gridPos + ivec3(0, 1, 1), gridSize);

    // calculate flag indicating if each vertex is inside or outside isosurface
    int cubeindex;
    cubeindex =  int(field[0] < isoValue);
    cubeindex += int(field[1] < isoValue)*2;
    cubeindex += int(field[2] < isoValue)*4;
    cubeindex += int(field[3] < isoValue)*8;
    cubeindex += int(field[4] < isoValue)*16;
    cubeindex += int(field[5] < isoValue)*32;
    cubeindex += int(field[6] < isoValue)*64;
    cubeindex += int(field[7] < isoValue)*128;

    // read number of vertices from texture
    int numVerts = texelFetch(numVertsTex, cubeindex, 0).a;
    return numVerts;
}

// compute interpolated vertex along an edge
vec3 vertexInterp(float isolevel, vec3 p0, vec3 p1, float f0, float f1)
{
    float t = (isolevel - f0) / (f1 - f0);
    return mix(p0, p1, t);
}

// compute interpolated vertex position and normal along an edge
void vertexInterp2(float isolevel, vec3 p0, vec3 p1, vec4 f0, vec4 f1, inout vec3 p, inout vec3 n)
{
    float t = (isolevel - f0.w) / (f1.w - f0.w);
    p = mix(p0, p1, t);
    n.x = mix(f0.x, f1.x, t);
    n.y = mix(f0.y, f1.y, t);
    n.z = mix(f0.z, f1.z, t);
    //    n = normalize(n);
}

// calculate triangle normal
vec3 calcNormal(vec3 v0, vec3 v1, vec3 v2)
{
    vec3 edge0 = v1 - v0;
    vec3 edge1 = v2 - v0;
    // note - it's faster to perform normalization in vertex shader rather than here
    return cross(edge0, edge1);
}

// version that calculates flat surface normal for each triangle
void generateTriangles(ivec3 gridPos, vec3 p)
{

    // calculate cell vertex positions
    vec3 v[8];
    v[0] = p;
    v[1] = p + vec3(voxelSize.x, 0, 0);
    v[2] = p + vec3(voxelSize.x, voxelSize.y, 0);
    v[3] = p + vec3(0, voxelSize.y, 0);
    v[4] = p + vec3(0, 0, voxelSize.z);
    v[5] = p + vec3(voxelSize.x, 0, voxelSize.z);
    v[6] = p + vec3(voxelSize.x, voxelSize.y, voxelSize.z);
    v[7] = p + vec3(0, voxelSize.y, voxelSize.z);

    float field[8];
    field[0] = sampleVolume(gridPos, gridSize);
    field[1] = sampleVolume(gridPos + ivec3(1, 0, 0), gridSize);
    field[2] = sampleVolume(gridPos + ivec3(1, 1, 0), gridSize);
    field[3] = sampleVolume(gridPos + ivec3(0, 1, 0), gridSize);
    field[4] = sampleVolume(gridPos + ivec3(0, 0, 1), gridSize);
    field[5] = sampleVolume(gridPos + ivec3(1, 0, 1), gridSize);
    field[6] = sampleVolume(gridPos + ivec3(1, 1, 1), gridSize);
    field[7] = sampleVolume(gridPos + ivec3(0, 1, 1), gridSize);

    // recalculate flag
    int cubeindex;
    cubeindex =  int(field[0] < isoValue);
    cubeindex += int(field[1] < isoValue)*2;
    cubeindex += int(field[2] < isoValue)*4;
    cubeindex += int(field[3] < isoValue)*8;
    cubeindex += int(field[4] < isoValue)*16;
    cubeindex += int(field[5] < isoValue)*32;
    cubeindex += int(field[6] < isoValue)*64;
    cubeindex += int(field[7] < isoValue)*128;

    // output triangle vertices
    if((cubeindex == 0) || (cubeindex == 255)) return;

    int numVerts = texelFetch(numVertsTex, cubeindex, 0).r;

    // find the vertices where the surface intersects the cube
    vec3 vertlist[12];

    vertlist[ 0] = vertexInterp(isoValue, v[0], v[1], field[0], field[1]);
    vertlist[ 1] = vertexInterp(isoValue, v[1], v[2], field[1], field[2]);
    vertlist[ 2] = vertexInterp(isoValue, v[2], v[3], field[2], field[3]);
    vertlist[ 3] = vertexInterp(isoValue, v[3], v[0], field[3], field[0]);

    vertlist[ 4] = vertexInterp(isoValue, v[4], v[5], field[4], field[5]);
    vertlist[ 5] = vertexInterp(isoValue, v[5], v[6], field[5], field[6]);
    vertlist[ 6] = vertexInterp(isoValue, v[6], v[7], field[6], field[7]);
    vertlist[ 7] = vertexInterp(isoValue, v[7], v[4], field[7], field[4]);

    vertlist[ 8] = vertexInterp(isoValue, v[0], v[4], field[0], field[4]);
    vertlist[ 9] = vertexInterp(isoValue, v[1], v[5], field[1], field[5]);
    vertlist[10] = vertexInterp(isoValue, v[2], v[6], field[2], field[6]);
    vertlist[11] = vertexInterp(isoValue, v[3], v[7], field[3], field[7]);


    for (int i=0; i<numVerts; i+=3)
    {
        //int index = numVertsScanned[voxel] + i;

        vec3 v[3];
        int edge;
        edge = texelFetch(triTex, (cubeindex*16) + i    , 0).r;
        v[0] = vertlist[edge];

        edge = texelFetch(triTex, (cubeindex*16) + i + 1, 0).r;
        v[1] = vertlist[edge];

        edge = texelFetch(triTex, (cubeindex*16) + i + 2, 0).r;
        v[2] = vertlist[edge];

        // calculate triangle surface normal
        vec3 n = normalize(calcNormal(v[0], v[1], v[2]));

        // Draw triangle in order 2->1->0 for face culling
        //aNorm = n;
        gs_out.FragPos = v[2]/(voxelSize/voxelSize.y); // restore the scaled world coordinate -> correct lighting behaviour
        gs_out.Normal = n;
        gs_out.TexCoords = v[2] + 0.5f*voxelSize;
        gl_Position = projectionMatrix*vec4(v[2], 1.0); // tri vertex 0
        EmitVertex();

        //aNorm = n;
        gs_out.FragPos = v[1]/(voxelSize/voxelSize.y);
        gs_out.Normal = n;
        gs_out.TexCoords = v[1] + 0.5f*voxelSize;
        gl_Position = projectionMatrix*vec4(v[1], 1.0); // tri vertex 1
        EmitVertex();

        //aNorm = n;
        gs_out.FragPos = v[0]/(voxelSize/voxelSize.y);
        gs_out.Normal = n;
        gs_out.TexCoords = v[0] + 0.5f*voxelSize;
        gl_Position = projectionMatrix*vec4(v[0], 1.0); // tri vertex 2
        EmitVertex();
        EndPrimitive();
    }
}

void render_quad(vec4 position);
void render_point(vec3 color, vec4 position);
void main() {
    //render_quad(gl_in[0].gl_Position);
    //render_point(vec3(1.0)*classifyVoxel()/15.0f, gl_in[0].gl_Position);
    //render_point(vec3(1.0)*sampleVolume(gs_in[0].gridPos, gridSize), gl_in[0].gl_Position); // Debug volume texture
    //render_point(vec3(1.0)*texelFetch(numVertsTex, 1, 0).r, gl_in[0].gl_Position); // Debug numverts texture
    //gs_out.Color = gs_in[0].gridPos*voxelSize/4.0;
    //gs_out.Color = vec3(0.7);
    generateTriangles(gs_in[0].gridPos, gl_in[0].gl_Position.xyz);
}

void render_quad(vec4 position)
{
    gs_out.TexCoords = gs_in[0].gridPos*voxelSize; // gs_in[0] since there's only one input vertex
    gl_Position = position + vec4(voxelSize*vec3(0, 0, 0.0),0.0); // 1:bottom-left
    gl_Position = projectionMatrix*gl_Position;
    EmitVertex();
    gl_Position = position + vec4(voxelSize*vec3(1, 0, 0.0),0.0); // 2:bottom-right
    gl_Position = projectionMatrix*gl_Position;
    EmitVertex();
    gl_Position = position + vec4(voxelSize*vec3(1, 1, 0.0),0.0); // 4:top-right
    gl_Position = projectionMatrix*gl_Position;
    EmitVertex();
    gl_Position = position + vec4(voxelSize*vec3(0, 1, 0.0),0.0); // 3:top-left
    gl_Position = projectionMatrix*gl_Position;
    EmitVertex();
    gl_Position = position + vec4(voxelSize*vec3(0, 0, 0.0),0.0); // 1:bottom-left
    gl_Position = projectionMatrix*gl_Position;
    EmitVertex();
    EndPrimitive();
}

void render_point(vec3 color, vec4 position)
{
    gs_out.TexCoords = color; // gs_in[0] since there's only one input vertex
    gl_Position = projectionMatrix*position; // 1:bottom-left
    EmitVertex();
    EndPrimitive();
}
