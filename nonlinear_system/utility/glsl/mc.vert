#version 330 core
out VS_OUT {
    ivec3 gridPos;
} vs_out;

ivec3 fetch3D(int ind, ivec3 gridSize) {

    ivec3 gridPos;
    // works for arbitrary size
    int tileXY = gridSize.x*gridSize.y;
    int nElem = tileXY*gridSize.z;
    gridPos.x = (ind%tileXY)%gridSize.x;
    gridPos.y = (ind%tileXY)/gridSize.x;
    gridPos.z = (ind%nElem)/tileXY;
    return gridPos;

}

/*ivec3 fetch3Dpow2(int ind, ivec3 gridSize) {

    ivec3 gridPos;
    // only works for pow of 2 size
    gridPos.x = ind & gridSize.x;
    gridPos.y = (ind >> log2(gridSize.x)) & gridSize.y;
    gridPos.z = (ind >> (log2(gridSize.x) + log2(gridSize.y))) & gridSize.z;
    return gridPos;

}*/

uniform ivec3 gridSize;
uniform vec3 voxelSize;

void main()
{
    vs_out.gridPos = fetch3D(gl_InstanceID, gridSize);

    gl_Position = vec4(vs_out.gridPos*voxelSize, 1.0); // normalize to vec3(0.0 - 1.0) for texture lookup
}
