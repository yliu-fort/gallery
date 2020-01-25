#version 430
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image3D img_output;
layout(std140, binding = 1) buffer ParticleDistribution
{
    double f[];
};

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

uniform ivec3 gridSize;
uniform float time;


void main() {
  // get index in global work group i.e x,y position
  ivec3 gridPos = ivec3(gl_GlobalInvocationID.xyz);
  int linearIndex = gridPos.x + gridPos.y*gridSize.x + gridPos.z*gridSize.x*gridSize.y;
  // 1. Load from image
  //float alpha = imageLoad(img_output, gridPos).a;
  // 2. Load from array
  // Minimal load size -> 4*float
  //float alpha = f[5*linearIndex + 0].x;// fetch mat4 f[0].x
  //float alpha = f[5*linearIndex + 0].y;// should fail
  double alpha = f[linearIndex];// fetch mat4 f[1].x
  //float alpha = f[4*linearIndex + 2];// fetch mat4 f[2].x
  //float alpha = f[4*linearIndex + 3];// fetch mat4 f[3].x
  //float alpha = f[linearIndex];// fetch float
  //
  // interesting stuff happens here later
  //
  vec4 color = vec4(vec3(gridPos)/vec3(gridSize),float(alpha));

  // output to a specific pixel in the image
  imageStore(img_output, gridPos, color);
  //f[linearIndex] = 0.5f*(f[linearIndex] + sin(0.1f*time)+1);
}
