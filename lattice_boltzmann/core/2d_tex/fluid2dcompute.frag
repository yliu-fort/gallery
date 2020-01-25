#version 430 core
layout (location = 0) out vec4 uvrf;
layout (location = 1) out vec4 f1_4;
layout (location = 2) out vec4 f5_8;
//layout (location = 3) out vec4 occl;
//in vec2 TexCoords;

layout(binding = 0) uniform sampler2D uvrf_old;
layout(binding = 1) uniform sampler2D f1_4_old;
layout(binding = 2) uniform sampler2D f5_8_old;
//layout(binding = 3) uniform sampler2D occl_old;

layout(std140, binding = 0) buffer occlRead
{
    vec4 _occl[];
};
layout(std140, binding = 1) buffer occlWrite
{
    vec4 occl[];
};

uniform int nstep;

layout (std140, binding = 0) uniform LatticeConstants
{
    mat3 wi;
    mat3 cx;
    mat3 cy;
    float csSqrInv;
    float csSqr;
    float Reg;
    float u_max;
    float tau;
    int NX,NY;
};

int ind()
{
    return int(gl_FragCoord.x-0.5f) + NX*int(gl_FragCoord.y-0.5f);
}

float sumMat3(mat3 f);
float weightedSumMat3(mat3 w, mat3 f);
mat3 mulMat3(mat3 A, mat3 B);

// lattice boltzmann method
void streaming(inout mat3 f, inout float r, inout float u, inout float v);
mat3 fetchLattice(inout float r, inout float u, inout float v, ivec2 dir);
void fetchLatticev(inout float r, inout float u, inout float v, ivec2 dir);
mat3 compute_equilibrium(const float rho, const float u, const float v);
mat3 collision(const mat3 f, const mat3 feq, const float alpha, const float beta);

// entropy constraint
float compute_g(const mat3 f,const mat3 feq,const float a,const float b);
float compute_gradg(const mat3 f,const mat3 feq,const float a,const float b);
float constrain_entropy(const mat3 f,const mat3 feq,const float b,const float aold);

// boundary
vec4 setupBoundary(inout float uw,inout float vw,inout float rhow);

ivec2 relocUV(ivec2 v, int nx, int ny)
{
    return ivec2(v.x - ((v.x + NX)/NX - 1)*NX, v.y - ((v.y + NY)/NY - 1)*NY);
}
void bounce_back(inout mat3 f)
{
    f = mat3(f[0][0], f[1][0], f[1][1],
            f[0][1], f[0][2], f[2][1],
            f[2][2], f[1][2], f[2][0]);
}
float one_sided_diff(float x1, float x2, float x3, float sign)
{
    return sign*0.5f*(-3.0f*x1+4.0f*x2-x3);
}
float central_diff(float x1, float x2)
{
    return 0.5f*(x1-x2);
}

uniform float tau_gui;

void main()
{
    float tau = tau_gui + 0.5f;
    // Streaming
    mat3 f;
    float rho,u,v;
    streaming(f,rho,u,v);

    // Fetch boundary info
    //vec4 boundary = texelFetch(occl_old, ivec2(gl_FragCoord.x, gl_FragCoord.y), 0);
    vec4 boundary = _occl[ind()];
    float alpha_old = 2.0f*boundary.z; // will move to a new attachment later
    //boundary = setupBoundary(u,v,rho); // if dynamic

    // no-slip boundary
    float occlusion = boundary.w;
    float robin_constant = boundary.x;
    // bi = [1 4 5 2 3 8 9 6 7] depend on ci
    if(occlusion > 0.0f) // boundary node
    {
        // Initial equilibrium relaxation
        float scale= min(1.0f, float(nstep)/max(NX,NY));

        // Dirchlet
        if(robin_constant > 0.0f)
        {
            mat3 cidotu = cx*u*scale + cy*v*scale;
            f -= mulMat3(2.0f*1.0f*csSqrInv*cidotu, wi);
            bounce_back(f);

        }
        // Neumann
        if(robin_constant == 0.0f)
        {

            // interior node
            mat3 fi;
            float rhoi,ui,vi;


            // Here is the tricky part: using outlet boundary at top & bottom may lead to reverse flow
            // which will cause instabilities. using value@ivec2(0,2) improves stability.
            ivec2 offset = ivec2(0,0);
            if(gl_FragCoord.y >(NY-1))  offset = ivec2( 0, 1); // Top
            if(gl_FragCoord.y <1)       offset = ivec2( 0,-1); // Bottom
            if(gl_FragCoord.x <1)       offset = ivec2(-1, 0);// Left
            if(gl_FragCoord.x >(NX-1))  offset = ivec2( 1, 0); // Right
            fi = fetchLattice(rhoi,ui,vi,offset); // Right

            // Improved stability when reverse flow may appear
            /*if(gl_FragCoord.y <1)       {vi = min(0.0, vi); }// Bottom
            if(gl_FragCoord.y >(NY-1))  {vi = max(0.0, vi); } // Top
            if(gl_FragCoord.x <1)       {ui = min(0.0, ui); }// Left
            if(gl_FragCoord.x >(NX-1))  {ui = max(0.0, ui); } // Right*/

            // RBC
            //f = compute_equilibrium(rhoi, ui ,vi);

            // CBC
            mat3 X = mat3(u,rho,0,csSqr/rho,u,0,0,0,u);
            mat3 Y = mat3(v,0,rho,0,v,0,csSqr/rho,0,v);

            float cs = sqrt(csSqr);
            float P1 = 0.5f*csSqrInv;
            float P2 = 0.5f/rho/cs;
            //mat3 Pxinv= mat3(P1,0,P1,
            //                 -P2,0,P2,
            //                 0,1,0);

            // interior* node
            float rhoi_1,ui_1,vi_1;
            fetchLatticev(rhoi_1,ui_1,vi_1,ivec2(2,0));

            // y+ node
            float rho_yp,u_yp,v_yp;
            fetchLatticev(rho_yp,u_yp,v_yp,ivec2(0,-1));
            // y- node
            float rho_ym,u_ym,v_ym;
            fetchLatticev(rho_ym,u_ym,v_ym,ivec2(0,1));

            // gradients
            float drdx = one_sided_diff(rho, rhoi, rhoi_1, -1.0f);
            float dudx = one_sided_diff(u, ui, ui_1, -1.0f);
            float dvdx = one_sided_diff(v, vi, vi_1, -1.0f);
            float drdy = central_diff(rho_yp, rho_ym);
            float dudy = central_diff(u_yp,u_ym);
            float dvdy = central_diff(v_yp,v_ym);

            // for right outlet, L1 = 0
            vec3 Lx = vec3((u - cs)*(csSqr*drdx - cs*rho*dudx),u*dvdx,(u + cs)*(csSqr*drdx + cs*rho*dudx));
            Lx.x = 0.0;

            // construct vectors
            vec3 m = vec3(rho,u,v);

            // Apply CBC
            float gamma = 0.0f;
            m.x += -( P1*Lx.x+P1*Lx.z) - gamma*(v*drdy + rho*dvdy);
            m.y += -(-P2*Lx.x+P2*Lx.z) - gamma*(v*dudy);
            m.z += -(Lx.y)             - gamma*(csSqr/rho*drdy + v*dvdy);

            rho = m.x;
            u = m.y;
            v = m.z;

            f = compute_equilibrium(rho,u*scale,v*scale);

        }
    }
    else // interior node -> collision
    {
        // Compute rho, u and v
        rho = sumMat3(f);
        u = weightedSumMat3(cx, f)/rho;
        v = weightedSumMat3(cy, f)/rho;

        // External force
        //float Fx = 8f*csSqr*(tau-0.5f)*u_max/pow(float(NY-2), 2f);
        float Fx = 0;
        float Fy = 0;

        Fx = weightedSumMat3(mulMat3(wi, cx)*Fx*csSqrInv, cx);
        Fy = weightedSumMat3(mulMat3(wi, cy)*Fy*csSqrInv, cy);

        u += 0.5f*Fx/rho;
        v += 0.5f*Fy/rho;

        float ueq = u + (tau-0.5f)*Fx/rho;
        float veq = v + (tau-0.5f)*Fy/rho;

        // Collide
        mat3 feq = compute_equilibrium(rho, ueq, veq);
        float beta = 0.5f/(tau);

        //float alpha_old = 2.0f*boundary.z; // will move to a new attachment later
        float alpha = constrain_entropy(f, feq, beta, alpha_old);
        f = collision(f, feq, alpha, beta);

        boundary.z = min(1.0f,0.5f*alpha+0.005f); // update alpha value with relaxation as initial guess for next time step
    }

    // Output
    uvrf = vec4(u,v,rho,f[0][0]);
    f1_4 = vec4(f[0][1],f[0][2],f[1][0],f[1][1]);
    f5_8 = vec4(f[1][2],f[2][0],f[2][1],f[2][2]);
    occl[ind()] = boundary;

    // if dynamic boundaries, setup new boundary
    //vec4 dynbc = setupBoundary(u,v,rho);

}

void streaming(inout mat3 f, inout float r, inout float u, inout float v)
{
    // Streaming

    ivec2 cord[9];
    cord[0] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[0][0], cy[0][0]);cord[0] = relocUV(cord[0],NX, NY);

    cord[1] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[0][1], cy[0][1]);cord[1] = relocUV(cord[1],NX, NY);
    cord[2] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[0][2], cy[0][2]);cord[2] = relocUV(cord[2],NX, NY);
    cord[3] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[1][0], cy[1][0]);cord[3] = relocUV(cord[3],NX, NY);
    cord[4] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[1][1], cy[1][1]);cord[4] = relocUV(cord[4],NX, NY);

    cord[5] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[1][2], cy[1][2]);cord[5] = relocUV(cord[5],NX, NY);
    cord[6] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[2][0], cy[2][0]);cord[6] = relocUV(cord[6],NX, NY);
    cord[7] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[2][1], cy[2][1]);cord[7] = relocUV(cord[7],NX, NY);
    cord[8] = ivec2(gl_FragCoord.x, gl_FragCoord.y) - ivec2(cx[2][2], cy[2][2]);cord[8] = relocUV(cord[8],NX, NY);

    f[0][1] = texelFetch(f1_4_old, cord[1], 0).x; // f1
    f[0][2] = texelFetch(f1_4_old, cord[2], 0).y; // f2
    f[1][0] = texelFetch(f1_4_old, cord[3], 0).z; // f3
    f[1][1] = texelFetch(f1_4_old, cord[4], 0).w; // f4

    f[1][2] = texelFetch(f5_8_old, cord[5], 0).x; // f5
    f[2][0] = texelFetch(f5_8_old, cord[6], 0).y; // f6
    f[2][1] = texelFetch(f5_8_old, cord[7], 0).z; // f7
    f[2][2] = texelFetch(f5_8_old, cord[8], 0).w; // f8


    vec4 a = texelFetch(uvrf_old, cord[0], 0);
    r = a.z;
    u = a.x;
    v = a.y;
    f[0][0] = a.w;

}

mat3 fetchLattice(inout float r, inout float u, inout float v, ivec2 dir)
{
    mat3 f;
    // Streaming
    ivec2 cord = ivec2(gl_FragCoord.x, gl_FragCoord.y) - dir;cord = relocUV(cord,NX, NY);

    vec4 a = texelFetch(uvrf_old, cord, 0);
    vec4 b = texelFetch(f1_4_old, cord, 0); // f1
    vec4 c = texelFetch(f5_8_old, cord, 0); // f5

    f[0][0] = a.w; // f0

    f[0][1] = b.x; // f1
    f[0][2] = b.y; // f2
    f[1][0] = b.z; // f3
    f[1][1] = b.w; // f4

    f[1][2] = c.x; // f5
    f[2][0] = c.y; // f6
    f[2][1] = c.z; // f7
    f[2][2] = c.w; // f8

    // Compute rho, u and v
    r = a.z;
    u = a.x;
    v = a.y;

    return f;
}

void fetchLatticev(inout float r, inout float u, inout float v, ivec2 dir)
{

    // Streaming
    ivec2 cord = ivec2(gl_FragCoord.x, gl_FragCoord.y) - dir;cord = relocUV(cord,NX, NY);

    vec4 a = texelFetch(uvrf_old, cord, 0);

    // Compute rho, u and v
    r = a.z;
    u = a.x;
    v = a.y;
}

mat3 compute_equilibrium(const float rho, const float u, const float v)
{
    mat3 cidotu = cx*u + cy*v;
    mat3 feq = rho*(1.0+3.0*cidotu+4.5*mulMat3(cidotu, cidotu)-1.5*(u*u+v*v));
    return mulMat3(feq, wi);
}

mat3 collision(const mat3 f, const mat3 feq, const float alpha, const float beta)
{
    mat3 df = f - feq;
    float ab = alpha*beta;
    return (f - ab*df);
}

float sumMat3(mat3 f)
{
    return (  f[0].x + f[0].y + f[0].z
            + f[1].x + f[1].y + f[1].z
            + f[2].x + f[2].y + f[2].z );
}

float weightedSumMat3(mat3 w, mat3 f)
{
    return (  w[0].x*f[0].x + w[0].y*f[0].y + w[0].z*f[0].z
            + w[1].x*f[1].x + w[1].y*f[1].y + w[1].z*f[1].z
            + w[2].x*f[2].x + w[2].y*f[2].y + w[2].z*f[2].z);
}

mat3 mulMat3(mat3 A, mat3 B)
{
    return mat3(A[0].x*B[0].x,A[0].y*B[0].y,A[0].z*B[0].z,
            A[1].x*B[1].x,A[1].y*B[1].y,A[1].z*B[1].z,
            A[2].x*B[2].x,A[2].y*B[2].y,A[2].z*B[2].z);
}

float compute_g(
                  const mat3 f,
                  const mat3 feq,
                  const float a,const float b) {

        float result = 0.0f;
        mat3 c = collision(f,feq,a,b);

        for(int i = 0;i < 3; i++)
        {   for(int j = 0;j < 3; j++)
            {

                if(c[i][j] < 0.0f) {c[i][j] = 0.00001f;}
                result += c[i][j]*log(c[i][j]/wi[i][j]) - f[i][j]*log(f[i][j]/wi[i][j]);
            }
        }
        return result;
}

float compute_gradg(
                  const mat3 f,
                  const mat3 feq,
                  const float a,const float b)
{

        float result = 0.0f;
        mat3 c = collision(f,feq,a,b);

        for(int i = 0;i < 3; i++)
        {   for(int j = 0;j < 3; j++)
            {
                if(c[i][j] < 0.0f) {c[i][j] = 0.00001f;}
                result += -b*(f[i][j] - feq[i][j])*(log(c[i][j]/wi[i][j]) + 1.0f);
            }
        }
        return result;
}

float constrain_entropy(
        const mat3 f,
        const mat3 feq,
        const float b,
        const float alpha_old)
{
        // calculate deviation
        float amin=0.2f;
        float amax=alpha_old;
        float maxDeviation = 0.0f;
        for(int i = 0;i < 3; i++)
        {   for(int j = 0;j < 3; j++)
            {
                float deviation = abs(f[i][j]-feq[i][j])/feq[i][j];
                if(deviation > maxDeviation)
                    maxDeviation = deviation;
            }
        }

        // if deviation is too large
        //float stableDeviation = 0.2;
        if(maxDeviation < 0.1f) return amax;

        // compute G value
        float Gmin = compute_g(f,feq,amin,b);
        float Gmax = compute_g(f,feq,amax,b);
        float gradGmin = compute_gradg(f,feq,amin,b);
        float gradGmax = compute_gradg(f,feq,amax,b);
        if((Gmin*Gmax) > 0.0f) return amax;
        if(Gmin > 0.0f){ float tmp = amin;amin = amax;amax = tmp; }

        float a = 0.5f*(amin + amax);
        float da = abs(amax - amin);
        float a_o = a;
        //float da_o = da;
        float G = compute_g(f,feq,a,b);
        float gradG = compute_gradg(f,feq,a,b);

        int maxIter = 20;
        float tolerance = 0.0001f;
        for(int it = 0; it < maxIter; it++)
        {
            if( ( ((a-amax)*gradG-G)*((a-amin)*gradG-G) >= 0.0f )
            ||  ( abs(a_o*gradG-G-1.0f) > 1.0f ) )
            {
                // bisection
                //da_o = da;
                da = 0.5f*(amax - amin);
                a = amin-amax;
                if(amin == a) return a;
            }else
            {
                //da_o = da;
                da = G/gradG;
                a_o = a;
                a -= da;
                if(a_o == a) return a;
            }
            if(abs(da) < tolerance) return a;

            G = compute_g(f,feq,a,b);
            gradG = compute_gradg(f,feq,a,b);
            if(G < 0.0f) {amin = a;}
            else {amax = a;}
        }

        return amax;

}

vec4 setupBoundary(inout float uw,inout float vw,inout float rhow)
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


        // sphere grid
        /*float interval = NX/4.0;
        for(int n = 0; n < 3; n++)
        {
            float dist = distance(gl_FragCoord.xy, vec2(n*interval + NX/4.0+NX/11.0, NY/8.0));
            if(dist < 64.0*(NY/1024.0)) occlusion = 1.0;

            dist = distance(gl_FragCoord.xy, vec2(n*interval + NX/4.0+NX/11.0, 7.0*NY/8.0));
            if(dist < 64.0*(NY/1024.0)) occlusion = 1.0;
        }


        float dist = distance(gl_FragCoord.xy, vec2((NY)/2.0, NY/2.0));
        float dir = dot(gl_FragCoord.xy- vec2(NY/2.0, NY/2.0), vec2(-1.0,0.0));
        if((dist > (NY/2.0))&&dir>0) occlusion = 1.0;*/

        // plate
        //if(gl_FragCoord.x <(NX/16.0+1) && gl_FragCoord.x >(NX/16.0)
        //        && gl_FragCoord.y <(4.5*NY/8.0) && gl_FragCoord.y >(3.5*NY/8.0)) occlusion = 1.0;

        //Sphere
        //float dist = distance(gl_FragCoord.xy, vec2(NX/6.0, NY/2.0));
        //if(dist < 128.0*(NY/1024.0)) occlusion = 1.0;

        //NACA0015
        vec2 pivot = vec2(NX/6.0, NY/2.0);
       float dist = distance(gl_FragCoord.xy, pivot);
       float x1 = float(gl_FragCoord.x) - pivot.x; // 0 < dx < NX/12.0
       float y1 = float(gl_FragCoord.y) - pivot.y; // 0 < dy < f(x)
       float s = 0.01; // scale angular velocity
       float aoa = (25.0 * sin(s*nstep/180.0*3.1415926))/180.0*3.1415926;
       float omega = (25.0 * s * cos(s*nstep/180.0*3.1415926))/180.0*3.1415926/180.0*3.1415926;
       float dx = x1*cos(aoa) - y1*sin(aoa);
       float dy = x1*sin(aoa) + y1*cos(aoa);
       float c = NX/6.0; // chord
       if((dx >= 0.0f) && (dx <= c))
       {
           float x = dx/c; // chord
           float t = 0.15f; //  t gives the last two digits in the NACA 4-digit denomination divided by 100
           float surf = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*pow(x,3.0f) - 0.1015*pow(x,4.0f));
           if(abs(dy/c) < surf)
           {
               occlusion = 1.0;
               uw = -omega*y1;vw = omega*x1; // need to add relative velocity
               rhow=1.0f;
           }
       }


    }
    float robin_constant = 1.0; // 1.0 - dirchlet, -1.0 - neumann
    {
        // Right
        if(gl_FragCoord.x >(NX-1)) robin_constant = 0.0;
    }
    // Boundary velocity
    //if(occlusion > 0.0) {uw = 0.0f;vw = 0.0f;rhow=1.0f;}
    {
        // Right
        //if(gl_FragCoord.x >(NX-1)) {uw = 0.1; vw = 0.0;}

        // Left
        //if(gl_FragCoord.x <1) {uw = 0.1; vw = 0.0;}

        // Top
        //if(gl_FragCoord.y >(NY-1)) {init.x =u_max; init.y = 0.0;}

        // Bottom
        //if(gl_FragCoord.y <1) {uw =0.0; vw = 0.0;}

    }

    return vec4(robin_constant, 0.0f, 1.0f, occlusion);
}
