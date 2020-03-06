#version 330 core
out VS_OUT {
    ivec3 gridPos;
    vec3 normal;
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
uniform float depth; // normalized, [0, 1)
uniform ivec3 gridSize;
uniform vec3 voxelSize;

void main()
{
    // This only valid for Z dir, need to implement something to rotate the plane...
    vs_out.gridPos = fetch3D(gl_InstanceID, gridSize);
    vs_out.normal = vec3(0,0,1);

    gl_Position = vec4(vs_out.gridPos*voxelSize, 1.0);
}
