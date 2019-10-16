#version 330 core
layout (location = 0) out vec4 uvrf;
layout (location = 1) out vec4 f1_4;
layout (location = 2) out vec4 f5_8;
layout (location = 3) out vec4 occl;
// this
//in vec2 TexCoords;

uniform int NX;
uniform int NY;
uniform mat3 wi;
uniform mat3 cx;
uniform mat3 cy;
uniform float csSqrInv;
uniform float csSqr;
uniform float Reg;
uniform float u_max;
uniform float tau;
/*
layout (std140) uniform Lattice_Constants
{
    mat3 wi;
    mat3 cx;
    mat3 cy;
    float csSqrInv;
    float csSqr;
    float Reg;
    float u_max;
    float tau;
};*/

void debug();

float draw_sphere(vec2 center, float radius)
{
    float occlusion = 0.0f;
    float dist = distance(gl_FragCoord.xy, center);
    if(dist < radius*NY ) occlusion = 1.0f;
    return occlusion;
}

float draw_airfoil(vec2 center, float chord, float aoa)
{
    float occlusion = 0.0f;
    float dist = distance(gl_FragCoord.xy, center);
    float x1 = float(gl_FragCoord.x) - center.x; // 0 < dx < NX/12.0
    float y1 = float(gl_FragCoord.y) - center.y; // 0 < dy < f(x)
    float dx = x1*cos(aoa) - y1*sin(aoa);
    float dy = x1*sin(aoa) + y1*cos(aoa);
    if((dx >= 0.0f) && (dx <= chord))
    {
        float x = dx/chord; // chord
        float t = 0.15f; //  t gives the last two digits in the NACA 4-digit denomination divided by 100
        float surf = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*pow(x,3f) - 0.1015*pow(x,4f));
        if(abs(dy/chord) < surf) occlusion = 1.0;
    }
    return occlusion;
}

mat3 mulMat3(mat3 A, mat3 B)
{
    return mat3(A[0].x*B[0].x,A[0].y*B[0].y,A[0].z*B[0].z,
            A[1].x*B[1].x,A[1].y*B[1].y,A[1].z*B[1].z,
            A[2].x*B[2].x,A[2].y*B[2].y,A[2].z*B[2].z);
}

mat3 init_equilibrium(float u, float v, float rho)
{
    mat3 cidotu = cx*u + cy*v;
    mat3 workMatrix = rho*(1.0+3.0*cidotu+4.5*mulMat3(cidotu, cidotu)-1.5*(u*u+v*v));

    return mulMat3(workMatrix, wi);
}

vec3 taylor_green(float t, int NX, int NY, float rho0, float tau, float u_max, float csSqr)
{
    float pi = 3.141592654;
    float nu = csSqr*(tau-0.5);
    float kx = pi/(16.0);
    float ky = pi/(16.0);
    float td = 1.0/(nu*(kx*kx+ky*ky));
    float X = gl_FragCoord.x;
    float Y = gl_FragCoord.y;
    float ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    float uy =  u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);

    float P = -0.25*rho0*u_max*u_max*( (ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y) )*exp(-2.0*t/td);

    float rho = rho0+3.0*P;

    return vec3(ux, uy, rho);

}

vec3 zeros()
{
    return vec3(0.0,0.0,1.0);
}

vec3 perturb()
{
    float d2 = distance(gl_FragCoord.xy, vec2(NX/2.0, NY/2.0));
    float rho = 1.0 + 1.0*exp(-d2/10.0);
    return vec3(0.0,0.0,rho);
}

vec3 parabolic(float uc)
{
    return vec3(uc*4.0f*(gl_FragCoord.y-1.0f)/(float(NY)-1.0f)*(1.0f-(gl_FragCoord.y-1.0f)/(float(NY)-1.0f)),0.0f,1.0f);
}

vec4 setupBoundary(inout vec3 init)
{
    float occlusion = 0.0;
    {
        // Top
        if(gl_FragCoord.y >(NY-1)) occlusion = 1.0;

        // Bottom
        if(gl_FragCoord.y <1) occlusion = 1.0;

        // Right
        if(gl_FragCoord.x >(NX-1)) occlusion = 1.0;

        // Left
        if(gl_FragCoord.x <1) occlusion = 1.0;

        // plate
        //if(gl_FragCoord.x <(NX/16.0+1) && gl_FragCoord.x >(NX/16.0)
        //        && gl_FragCoord.y <(4.5*NY/8.0) && gl_FragCoord.y >(3.5*NY/8.0)) occlusion = 1.0;

        //Sphere
        occlusion += draw_sphere(vec2(NX/6.0,NY/2.0), 0.1f);

        //NACA0015
        //occlusion += draw_airfoil(vec2(NX/6.0, NY/2.0), NX/6.0, radians(5.0));

    }

    float robin_constant = 1.0f; // 1.0 - dirchlet, -1.0 - neumann
    {
        // Right
        if(gl_FragCoord.x >(NX-1)) robin_constant = 0.0f;

    }
    // Boundary velocity
    if(occlusion > 0.0) {init.x = 0.0f;init.y = 0.0f;}
    {
        // Right
        if(gl_FragCoord.x >(NX-1)) {init.x =u_max; init.y = 0.0;}
        // when using CBC outlet, make sure inlet

        // Left
        if(gl_FragCoord.x <1) {init.x =u_max; init.y = 0.0;}

        // Top
        if(gl_FragCoord.y >(NY-1)) {init.x =u_max; init.y = 0.0;}

        // Bottom
        if(gl_FragCoord.y <1) {init.x =u_max; init.y = 0.0;}

    }

    return vec4(robin_constant, 0.0f, 0.0f, occlusion);
}

void main()
{
    // textureSize(sampler2D, P) to retrieve NX, NY
    //vec3 initField =taylor_green(0, NX, NY, 1.0, tau, u_max, csSqr);
    vec3 initField =zeros();
    //vec3 initField =parabolic(u_max);

    initField.x += u_max;
    //initField += perturb();

    // Boundary
    occl = setupBoundary(initField);

    // compute equilibrium
    mat3 f = init_equilibrium(initField.x,initField.y,initField.z);

    // Output
    uvrf = vec4(initField, f[0][0]);
    f1_4 = vec4(f[0][1],f[0][2],f[1][0],f[1][1]);
    f5_8 = vec4(f[1][2],f[2][0],f[2][1],f[2][2]);


}

