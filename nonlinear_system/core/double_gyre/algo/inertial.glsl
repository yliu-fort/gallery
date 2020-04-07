#version 430
// Kernel
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

// Input
layout(std140, binding = 0) buffer readonlybuffer
{
    vec4 lhs[];
};

layout(std140, binding = 0) buffer readonlybuffer2
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

struct Param
{
    vec4 a_st[16];
    vec4 a_ct[16];
    vec4 k_st[16];
    vec4 k_ct[16];
    vec4 om[16];
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
vec4 p_vel;
float t;
vec4 f(vec4 p, float tc)
{
    vec4 dpdt = vec4(0,0,0,0);
    vec4 uMat = vec4(0);
    vec4 vMat = vec4(0);

    // A large L-vortex appeared near origin?
    // Add some constant in position to avoid that...
    vec2 _p = p.xy + vec2(32.0f,32.0f);

    for(int j = 0; j < 16; j++)
    {
            vec4 KS = _p.x*param.k_st[j] + _p.y*param.k_ct[j] + param.om[j]*tc;
            vec4 dKS = cos( KS )-sin( KS );

            uMat += param.a_ct[j]*dKS;
            vMat -= param.a_st[j]*dKS;
    }

    dpdt.x += uMat[0];dpdt.x += uMat[1];dpdt.x += uMat[2];dpdt.x += uMat[3];
    dpdt.y += vMat[0];dpdt.y += vMat[1];dpdt.y += vMat[2];dpdt.y += vMat[3];

    // Inertial correction
    vec4 d2pdt2 = tau_inv*(dpdt - p_vel) + g;

    return (tc-t)*d2pdt2;
}

void main()
{

    // get index in global work group i.e x,y position
    gridPos = ivec3(gl_GlobalInvocationID.xyz);
    linearIndex = ind();
    if(linearIndex >= dim.w) { return; }

    // Actual compute subroutine
    t = t0;
    float dt0 = t_end - t0;
    float dt = t_end - t0;

    vec4 a1,a2,a3,a4,a5,a6,dpdt,ds;
    vec4 p = lhs[linearIndex];
    p_vel = lhs_2[linearIndex];

    float s;

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
        {
            a1 = f(p,t);
            a2 = f(p+0.25f*a1,t+0.25f*dt);
            a3 = f(p+0.09375f*a1+0.28125f*a2,t+0.375f*dt);
            a4 = f(p+0.8793809740555303f*a1+3.277196176604461f*a2+3.3208921256258535f*a3,t+0.9230769230769231f*dt);
            a5 = f(p+2.0324074074074074f*a1-8.0f*a2+7.173489278752436f *a3-0.20589668615984405f*a4,t+dt);
            a6 = f(p-0.2962962962962963f*a1+2.0f*a2-1.3816764132553607f*a3+0.4529727095516569f *a4-0.275f*a5,t+0.5f*dt);

            dpdt = 0.11574074074074074f*a1 + 0.5489278752436647f*a3 + 0.5357229943916118f*a4 - 0.2f*a5;
        }

        // Euler method timestepping may leads unstability
        t = t + dt;
        p_vel += dpdt;
        p += dt*p_vel;

    }

    lhs[linearIndex] = vec4(p.x,p.y,0,1);
    lhs_2[linearIndex] = p_vel;
}


