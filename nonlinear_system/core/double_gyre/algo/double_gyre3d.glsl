#version 430
// Kernel
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

// Input
layout(std140, binding = 0) buffer readonlybuffer
{
    vec4 lhs[];
};

// Uniforms
layout (std140, binding = 0) uniform meshParameters
{
    ivec4 dim;
};

uniform float t0;
uniform float t_end;
uniform float tol;

// Global variables
ivec3 gridPos;
int linearIndex;

int ind()
{
    return (gridPos.x + gridPos.y*dim.x + gridPos.z*dim.x*dim.y);
}

// Compute subroutine
// Current: Double-gyre
vec4 f(vec4 p, float t)
{
    vec3 _p =mod(p.xyz,vec3(2.0f,1.0f,1.0f));
    vec3 dpdt;
    float pi = 3.141592654f;
    float eps1 = 0.25f, A = 0.1f, w = 2.0f*pi/10.0f;

    float a = eps1*sin(w*(t)), b = 1.0f - 2.0f*a;
    float fx = a*_p.x*_p.x + b*_p.x;
    float fp = 2*a*_p.x + b;

    dpdt.x = -pi*A*sin(pi*fx)*cos(pi*_p.y);
    dpdt.y =  pi*A*cos(pi*fx)*sin(pi*_p.y)*fp;
    dpdt.z =  pi*A*0.2f*_p.z*(1.0f-_p.z)*(_p.z - eps1*sin(4.0f*pi*w*t)-0.5f);

    return vec4(dpdt,0.0f);
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

    vec4 a1,a2,a3,a4,a5,a6,dpdt,ds;
    vec4 p = lhs[linearIndex];

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
}


