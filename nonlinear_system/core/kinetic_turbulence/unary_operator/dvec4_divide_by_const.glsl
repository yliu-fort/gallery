#version 430
// Kernel
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

// Input
layout(std140, binding = 1) buffer writeonlybuffer
{
    vec4 lhs[];
};

// Uniforms
uniform vec4 rhs;


layout (std140, binding = 0) uniform meshParameters
{
    ivec4 dim;
};

// Global variables
ivec3 gridPos;
int linearIndex;

int ind()
{
    return (gridPos.x + gridPos.y*dim.x + gridPos.z*dim.x*dim.y);
}

void main()
{

    // get index in global work group i.e x,y position
    gridPos = ivec3(gl_GlobalInvocationID.xyz);
    linearIndex = ind();
    if(linearIndex >= dim.w) { return; }

    // Actual compute subroutine
    lhs[linearIndex] /= rhs;

}


