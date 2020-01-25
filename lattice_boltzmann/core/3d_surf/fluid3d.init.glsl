#version 430
layout(local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

// Definition
struct DiscreteLattice{
    vec4 v[7];
};

// Input
layout(std140, binding = 0) buffer uvwrRead
{
    vec4 _uvwr[];
};

layout(std140, binding = 1) buffer occlRead
{
    vec4 _occl[];
};
layout(std140, binding = 2) buffer uvwrWrite
{
    vec4 uvwr[];
};
layout(std140, binding = 3) buffer occlWrite
{
    vec4 occl[];
};
// f // fread: 11 12 13 14 15 16 17 fwrite: 4 5 6 7 8 9 10
layout(std140, binding = 4) buffer fFieldWrite0
{
    mat2x4 _fw0[];
};
layout(std140, binding = 5) buffer fFieldWrite1
{
    mat2x4 _fw1[];
};
layout(std140, binding = 6) buffer fFieldWrite2
{
    mat2x4 _fw2[];
};
layout(std140, binding = 7) buffer fFieldWrite3
{
    vec4 _fw3[];
};

layout(std140, binding = 8) buffer fFieldRead0
{
    mat2x4 _fr0[];
};
layout(std140, binding = 9) buffer fFieldRead1
{
    mat2x4 _fr1[];
};
layout(std140, binding = 10) buffer fFieldRead2
{
    mat2x4 _fr2[];
};
layout(std140, binding = 11) buffer fFieldRead3
{
    vec4 _fr3[];
};

layout(rgba32f, binding = 0) uniform image3D vis;


// Uniform
layout (std140, binding = 0) uniform LatticeConstants
{
    DiscreteLattice wi;
    DiscreteLattice cx;
    DiscreteLattice cy;
    DiscreteLattice cz;
    int NX,NY,NZ,nElem;
    float csSqrInv;
    float csSqr;
    float Reg;
    float u_max;
    float tau;

};

layout (std140, binding = 1) uniform RuntimeParameters
{
    ivec4 ntasks;
    int time;
};

// Functions
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

int ind(ivec3 gridPos)
{
    return (gridPos.x + gridPos.y*NX + gridPos.z*NX*NY);
}

ivec3 relocUV(ivec3 v)
{
    return ivec3(v.x - ((v.x + NX)/NX - 1)*NX,
                 v.y - ((v.y + NY)/NY - 1)*NY,
                 v.z - ((v.z + NZ)/NZ - 1)*NZ);
}

float draw_sphere(vec3 center, float radius, ivec3 gridPos)
{
    float occlusion = 0.0f;
    float dist = distance(vec3(gridPos.xyz), center);
    if(dist < radius ) occlusion = 1.0f;
    return occlusion;
}

float draw_cylinderZ(vec2 center, float radius, ivec3 gridPos)
{
    float occlusion = 0.0f;
    float dist = distance(vec2(gridPos.xy), center);
    if(dist < radius ) occlusion = 1.0f;
    return occlusion;
}

float clipZ(float a, float z, ivec3 gridPos)
{
    if((z > 0 && gridPos.z < z) || (z < 0 && gridPos.z > -z)) return a;
    return 0;
}

float draw_airfoil2d(vec2 center, float chord, float aoa, ivec3 gridPos)
{
    float occlusion = 0.0f;
    float dist = distance(gridPos.xy, center);
    float x1 = float(gridPos.x) - center.x; // 0 < dx < NX/12.0
    float y1 = float(gridPos.y) - center.y; // 0 < dy < f(x)
    float dx = x1*cos(aoa) - y1*sin(aoa);
    float dy = x1*sin(aoa) + y1*cos(aoa);
    if((dx >= 0.0f) && (dx <= chord))
    {
        float x = dx/chord; // chord
        float t = 0.15f; //  t gives the last two digits in the NACA 4-digit denomination divided by 100
        float surf = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*pow(x,3.0f) - 0.1015*pow(x,4.0f));
        if(abs(dy/chord) < surf) occlusion = 1.0;
    }
    return occlusion;
}

DiscreteLattice init_equilibrium(vec4 init, float s)
{
    float u=init.x, v=init.y, w=init.z, rho=init.w;
    DiscreteLattice workMatrix;
    for(int i = 0; i < 7; i++)
    {
        vec4 cidotu = cx.v[i]*u*s + cy.v[i]*v*s + cz.v[i]*w*s;
        workMatrix.v[i] = wi.v[i]*rho*(1.0f+3.0f*cidotu+4.5f*cidotu*cidotu-1.5f*(u*u+v*v+w*w));
    }

    return workMatrix;
}

vec4 taylor_green2d(float t, float rho0, ivec3 gridPos)
{
    float pi = 3.141592654;
    float nu = csSqr*(tau-0.5);
    float kx = pi/(16.0);
    float ky = pi/(16.0);
    float td = 1.0/(nu*(kx*kx+ky*ky));
    float X = gridPos.x;
    float Y = gridPos.y;
    float ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    float uy =  u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
    float P = -0.25*rho0*u_max*u_max*( (ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y) )*exp(-2.0*t/td);
    float rho = rho0+3.0*P;

    return vec4(ux, uy, 0.0, rho);
}

vec4 zeros()
{
    return vec4(0.0,0.0,0.0,1.0);
}

vec4 perturb3d(vec3 center, ivec3 gridPos)
{
    float d2 = distance(gridPos.xyz, center);
    float rho = 1.0 + 1.0*exp(-d2/10.0);
    return vec4(0.0,0.0,0.0,rho);
}

float gaussian(vec2 center, float sigma,vec2 gridPos)
{
    float d2 = distance(center, gridPos);
    return exp(-d2*d2/sigma/sigma);
}
float circle(vec2 center, float radius,vec2 gridPos)
{
    float d2 = distance(center, gridPos);
    return (d2<radius?1.0:0.0);
}

vec4 parabolic(float uc, ivec3 gridPos)
{
    return vec4(uc*4.0f*(gridPos.y-1.0f)/(float(NY)-1.0f)*(1.0f-(gridPos.y-1.0f)/(float(NY)-1.0f)),0.0f,0.0f,1.0f);
}

float sinusoidal(float uc, float amp, float k, ivec3 gridPos)
{
    return uc*(1.0f-amp+amp*sin(2.0f*3.14159f*gridPos.y/float(NY/k)));
}

vec4 setupBoundary(inout vec4 init, ivec3 gridPos)
{
    float occlusion = 0.0;
    {
        // Boundaries
        if(gridPos.x >(NX-2)) occlusion = 1.0;// Right
        if(gridPos.x <1) occlusion = 1.0;// Left
        if(gridPos.y >(NY-2)) occlusion = 1.0;// Top
        if(gridPos.y <1) occlusion = 1.0;// Bottom
        if(gridPos.z >(NZ-2)) occlusion = 1.0;// Back
        if(gridPos.z <1) occlusion = 1.0;// Front

        // vertical plate
        /*
        if(gridPos.x <min(NY/2.0,NX/8.0)+1 && gridPos.x >min(NY/2.0,NX/8.0)
                && gridPos.y <(1.0*NY/8.0) && gridPos.y >(0.5*NY/8.0)) occlusion = 1.0;
        if(gridPos.x <min(NY/2.0,NX/8.0)+1 && gridPos.x >min(NY/2.0,NX/8.0)
                && gridPos.y >(7.0*NY/8.0) && gridPos.y <(7.5*NY/8.0)) occlusion = 1.0;
        if(gridPos.x <min(NY/2.0,NX/8.0)+1 && gridPos.x >min(NY/2.0,NX/8.0)
                && gridPos.y >(2.5*NY/8.0) && gridPos.y <(3.0*NY/8.0)) occlusion = 1.0;
        if(gridPos.x <min(NY/2.0,NX/8.0)+1 && gridPos.x >min(NY/2.0,NX/8.0)
                && gridPos.y >(5.0*NY/8.0) && gridPos.y <(5.5*NY/8.0)) occlusion = 1.0;
        */
        // vertical grid: to trigger flow separation, grid width must >= 4
        //if(gridPos.x <min(16.0,NX/16.0)+4 && gridPos.x >min(16.0,NX/16.0)
        //        && (int(gridPos.y+6)%41 <8) && (int(gridPos.z+6)%41 <8)) occlusion = 1.0;
        //if(gridPos.x <min(16.0,NX/16.0)+42 && gridPos.x >min(16.0,NX/16.0)+40
        //        && (int(gridPos.y)%23 <4) && (int(gridPos.z)%23 <4)) occlusion = 1.0;
        //if(gridPos.x <min(16.0,NX/16.0)+82 && gridPos.x >min(16.0,NX/16.0)+80
        //        && ((int(gridPos.y)%13 <2) && (int(gridPos.z)%13 <2))) occlusion = 1.0;

        // vertical assembled grid
        //if(gridPos.x <min(16.0,NX/16.0)+14 && gridPos.x >min(16.0,NX/16.0)+10
        //        && ((int(gridPos.y+6)%41 <4) || (int(gridPos.z+6)%41 <4))) occlusion = 1.0;
        //if(gridPos.x <min(16.0,NX/16.0)+52 && gridPos.x >min(16.0,NX/16.0)+50
        //        && ((int(gridPos.y)%23 <2) || (int(gridPos.z)%23 <2))) occlusion = 1.0;
        //if(gridPos.x <min(16.0,NX/16.0)+92 && gridPos.x >min(16.0,NX/16.0)+90
        //        && ((int(gridPos.y)%13 <1) || (int(gridPos.z)%13 <1))) occlusion = 1.0;

        // horizontal plate
        //if(gridPos.y <16 && gridPos.y >0
        //       /*&& gridPos.x <(3.0*NX/4.0)*/ && gridPos.x >(NX/8.0)) occlusion = 1.0;

        //Sphere
        //occlusion += draw_sphere(vec3(NX/4.0,NY/2.0,NZ/2.0), 0.2f*NY, gridPos);

        //Cylinder-Zdir
        occlusion += clipZ(clipZ(
                               draw_cylinderZ(vec2(NX/2.0,0.1f*NY), 0.2f*NY, gridPos),
                               0.7*NZ, gridPos),-0.05*NZ, gridPos);

        //NACA0015
        //occlusion += clipZ(clipZ(
        //                draw_airfoil2d(vec2(NX/8.0, NY/2.0+0.05*NZ), 0.5f*NY, radians(15.0),  gridPos),
        //                0.7*NZ, gridPos),-0.05*NZ, gridPos);
    }

    float robin_constant = 1.0f; // 1.0 - dirchlet, 0.0 - neumann
    {
        //if(gridPos.x >(NX-2)) robin_constant = 0.0f; // Right
    }
    // Boundary velocity
    if(occlusion > 0.0) {init.x = 0.0f;init.y = 0.0f;init.z = 0.0;}
    {

        if(gridPos.z <1) {init = vec4(0,0,0,1);} // Front
        if(gridPos.z >(NZ-2)) {init = vec4(0,0,0,1);} // Back
        if(gridPos.y <1) {init = vec4(0,0,0,1);} // Bottom
        if(gridPos.x <1) {init = vec4(0,0,0,1);} // Left
        if(gridPos.x >(NX-2)) {init = vec4(0,0,0,1);} // Right

        // For corner stability, interface between dirchlet boundary assigned with
        // different velocities, larger value must be assigned to corner cell.

        if(gridPos.y >(NY-2)) {init = vec4(0,0,0,1);} // Top

    }

    // Free-surface
    float mass_fraction = 0.001f;

    {
        // Rect region
        if(gridPos.x <(NX/4.0) &&
           gridPos.y <(NY/1.2) &&
           gridPos.z <(NZ/1.2)) mass_fraction = 1.0f;// Left

        // No mass fraction on bc (for clean visualization)
        if(occlusion > 0.0) mass_fraction = 0.001f;// Left
    }

    return vec4(robin_constant, mass_fraction, 1.0f, occlusion);
}

void main() {

    // get index in global work group i.e x,y position
    ivec3 gridPos = ivec3(gl_GlobalInvocationID.xyz) + ntasks.xyz*ntasks.w;
    if(gridPos.x >= NX || gridPos.y >= NY || gridPos.z >= NZ) return;
    int linearIndex = gridPos.x + gridPos.y*NX + gridPos.z*NX*NY;

    // Initial vel and rho
    vec4 initField =zeros();
    //initField =perturb3d(vec3(NX/4.0,	NY/4.0,NZ/4.0),gridPos);
    //vec4 initField =taylor_green2d(0,1.0,gridPos);

    //Boundary
    vec4 initBC = setupBoundary(initField, gridPos);

    // Equilibrium: relaxation forbiddened for intial condition sensitive simulation
    float relaxation = 1.0f/max(256,max(NZ, max(NX,NY)));
    //relaxation = 1.0f;
    DiscreteLattice f = init_equilibrium(initField, relaxation);

    // output
    uvwr[linearIndex] = initField;
    occl[linearIndex] = initBC;

    for(int j = 0; j < 4; j++){int i = 0; _fw0[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 1; _fw0[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 2; _fw1[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 3; _fw1[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 4; _fw2[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 5; _fw2[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j]; }
    for(int j = 0; j < 4; j++){int i = 6; _fw3[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))]   [j] = f.v[i][j]; }

    // Post-processing
    vec4 color = vec4(vec3(sqrt(initField.x*initField.x+
                                initField.y*initField.y+
                                initField.z*initField.z)/u_max),initBC.w);
    //vec4 color = vec4(f.v[0].yz, initBC.w);
    imageStore(vis, gridPos, color);
    //vis[linearIndex] = color;

}
