#version 440

// Kernel
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

// Input
layout(std140, binding = 0) buffer readonlybuffer
{
    vec4 lhs[];
};

layout(std140, binding = 1) buffer readonlybuffer2
{
    vec4 lhs_2[];
};

// Uniforms
layout (std140, binding = 0) uniform meshParameters
{
    ivec4 dim;
};

uniform float t0;
uniform float t_end;
uniform float tol;

#define NK (16)

struct Param
{
    vec4 a_xt[NK];
    vec4 a_yt[NK];
    vec4 a_zt[NK];
    vec4 k_xt[NK];
    vec4 k_yt[NK];
    vec4 k_zt[NK];
    vec4 om[NK];
};

layout (std140, binding = 1) uniform kineticParameters
{
    Param param;
};

uniform float tau_inv;
uniform vec4 g;

// Global variables
ivec3 gridPos;
int linearIndex;

int ind()
{
    return (gridPos.x + gridPos.y*dim.x + gridPos.z*dim.x*dim.y);
}

// Compute subroutine
// Evaluate velocity
vec4 f(vec4 p, float t)
{
    vec4 dpdt = vec4(0,0,0,0);
    vec4 uMat = vec4(0);
    vec4 vMat = vec4(0);
    vec4 wMat = vec4(0);

    // A large L-vortex appeared near origin?
    // Add some constant in position to avoid that...
    vec3 _p = p.xyz + vec3(32.0f,32.0f,32.0f);

    // Unrolled loops
    for(int i = 0; i < NK; i++)
    {
        vec4 KS,dKS;

        KS = _p.x*param.k_xt[i] + _p.y*param.k_yt[i] + _p.z*param.k_zt[i] + param.om[i]*t;
        dKS = cos( KS ) - sin( KS );

        uMat += param.a_xt[i]*dKS;
        vMat += param.a_yt[i]*dKS;
        wMat += param.a_zt[i]*dKS;

    }

    dpdt.x += dot(uMat,vec4(1));
    dpdt.y += dot(vMat,vec4(1));
    dpdt.z += dot(wMat,vec4(1));

    return dpdt;
}

void main()
{

    // get index in global work group i.e x,y position
    gridPos = ivec3(gl_GlobalInvocationID.xyz);
    linearIndex = ind();
    if(linearIndex >= dim.w) { return; }

    // Actual compute subroutine
    float t = t0;
    float dt0 = t_end - t0;
    float dt = t_end - t0;
    dt = dt > 0 ? dt:1;

    vec4 a1,a2,a3,a4,a5,a6,dpdt,ds;
    vec4 p = lhs[linearIndex];

    float s;
    dpdt = f(p,t);

    while(t < t_end)
    {
        // 4/5th Runge-kutta adaptive
        for(int inner_cycle = 0; inner_cycle < 16; inner_cycle++)
        {
            a1 = dt*f(p,t);
            a2 = dt*f(p+0.25f*a1,t+0.25f*dt);
            a3 = dt*f(p+0.09375f*a1+0.28125f*a2,t+0.375f*dt);
            a4 = dt*f(p+0.8793809740555303f*a1+3.277196176604461f*a2+3.3208921256258535f*a3,t+0.9230769230769231f*dt);
            a5 = dt*f(p+2.0324074074074074f*a1-8.0f              *a2+7.1734892787524360f*a3-0.20589668615984405f*a4,t+dt);
            a6 = dt*f(p-0.2962962962962963f*a1+2.0f              *a2-1.3816764132553607f*a3+0.4529727095516569f *a4-0.275f*a5,t+0.5f*dt);

            // Residual estimation
            ds = 0.002777777777778f*a1-0.029941520467836f*a3-0.029591504049595f*a4-0.02f*a5+0.036363636363636f*a6;

            ds = abs(ds);
            ds.x = ds.x > ds.y?ds.x:ds.y;
            ds.x = ds.x > ds.z?ds.x:ds.z;
            ds.x = ds.x > ds.w?ds.x:ds.w;
            ds.x += 1e-6;

            s = 0.84f*pow(tol/ds[0],0.25f);

            if (s < 1) // h is too large
            {
                dt = s*dt;
            }
            else
            {
                if (s < 2) // h is appropriate
                {
                    if(t + dt > t_end)
                    {
                        dt = t_end - t;
                    }
                    break;
                }
                else // s >= 2, h is too small
                {
                    dt = min(s*dt,dt0);
                }
            }
        }

        // Compute optimal dpdt
        // Caution: must store dt*f or you will have numerical accuracy issue (truncation)
        {
            a1 = dt*f(p,t);
            a2 = dt*f(p+0.25f*a1,t+0.25f*dt);
            a3 = dt*f(p+0.09375f*a1+0.28125f*a2,t+0.375f*dt);
            a4 = dt*f(p+0.8793809740555303f*a1+3.277196176604461f*a2+3.3208921256258535f*a3,t+0.9230769230769231f*dt);
            a5 = dt*f(p+2.0324074074074074f*a1-8.0f*a2+7.173489278752436f *a3-0.20589668615984405f*a4,t+dt);
            a6 = dt*f(p-0.2962962962962963f*a1+2.0f*a2-1.3816764132553607f*a3+0.4529727095516569f *a4-0.275f*a5,t+0.5f*dt);

            dpdt = 0.11574074074074074f*a1 + 0.5489278752436647f*a3 + 0.5357229943916118f*a4 - 0.2f*a5;
        }

        t = t + dt;
        p = p + dpdt;
    }

    lhs[linearIndex] = p;
    lhs_2[linearIndex] = dpdt/dt; // velocity
}


