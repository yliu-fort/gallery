#version 430
// Kernel
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

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
layout(std140, binding = 4) buffer fFieldRead0
{
    mat3 _fr0[];
};
layout(std140, binding = 5) buffer fFieldWrite0
{
    mat3 _fw0[];
};


layout(rgba32f, binding = 0) uniform image3D vis;

// Visualization
layout(binding = 15) uniform sampler1D colormap;

// Uniform
layout (std140, binding = 0) uniform LatticeConstants
{
    mat3 wi;
    mat3 cx;
    mat3 cy;
    //mat3 cz;
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
    int nstep;
};

// Global variables
ivec3 gridPos = ivec3(0);
int index=0;

// Functions
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

// lattice boltzmann method
//void streaming(inout mat3 f, inout float r, inout float u, inout float v, inout float w, ivec3 gridPos);
//mat3 fetchLattice(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos);
//void fetchLatticev(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos);
mat3 compute_equilibrium(const float rho, const float u, const float v, const float w);
mat3 collision(const mat3 f, const mat3 feq, const float alpha, const float beta);

// entropy constraint
float compute_g(const mat3 f, const mat3 feq, const float a, const float b);
float compute_gradg(const mat3 f, const mat3 feq, const float a, const float b);
float constrain_entropy(const mat3 f, const mat3 feq, const float b, const float alpha_old);

// boundary
//vec4 setupBoundary(inout float uw,inout float vw,inout float rhow);

ivec3 relocUV(ivec3 v)
{
    return ivec3(v.x - ((v.x + NX)/NX - 1)*NX,
                 v.y - ((v.y + NY)/NY - 1)*NY,
                 v.z - ((v.z + NZ)/NZ - 1)*NZ);
}
int ind(ivec3 p)
{
    ivec3 _gridPos = relocUV(gridPos + p);
    return (_gridPos.x + _gridPos.y*NX + _gridPos.z*NX*NY);
}

// set a bi later...
//const int bi[28] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25,27}; // D3Q27
//const int bi[20] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,19}; // D3Q19
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


void streaming(inout mat3 f, inout float r, inout float u, inout float v, inout float w)
{
    // Streaming
    f = _fr0[index];

    //vec4 a = texelFetch(uvwr_old, gridPos, 0);
    vec4 a = _uvwr[index];
    u = a.x;
    v = a.y;
    w = 0.0f;
    r = a.w;

}

mat3 fetchLattice(inout float r, inout float u, inout float v, inout float w, ivec3 dir)
{
    mat3 f;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        { f[i][j] = _fr0[ind(ivec3(vec3(cx[i][j],cy[i][j],0)) + dir)][i][j];}

    //vec4 a = texelFetch(uvwr_old, relocUV(gridPos - dir), 0);
    vec4 a = _uvwr[ind(dir)];
    u = a.x;
    v = a.y;
    w = 0.0f;
    r = a.w;

    return f;
}

void fetchLatticev(inout float r, inout float u, inout float v, inout float w, ivec3 dir)
{
    vec4 a = _uvwr[ind(dir)];
    u = a.x;
    v = a.y;
    w = 0.0f;
    r = a.w;
}

int corner_detector()
{
    int flag=0,width=1;
    // Artificially increase tau
    if(gridPos.y <2+width || gridPos.y >(NY-3-width))       flag += 1;// Next to Bottom Face
    if(gridPos.x <2+width || gridPos.x >(NX-3-width))       flag += 1;// Next to Left Face
    //if(gridPos.z <2+width || gridPos.z >(NZ-3-width))       flag += 1;// Next to Front Face

    // Corner
    return flag;

}

void corner_stabilizer(inout float t)
{
    // Corner
    if(corner_detector() > 1) t = max(0.55f, t);

    // Edge
    if(corner_detector() > 0) t = max(0.505f,t);

}

mat3 vorticity()
{
    // x+ node
    mat3 p,m;
    p[0] = _uvwr[ind(ivec3(1,0,0))].xyz;
    p[1] = _uvwr[ind(ivec3(0,1,0))].xyz;
    p[2] = _uvwr[ind(ivec3(0,0,1))].xyz;
    m[0] = _uvwr[ind(ivec3(-1,0,0) )].xyz;
    m[1] = _uvwr[ind(ivec3(0,-1,0) )].xyz;
    m[2] = _uvwr[ind(ivec3(0,0,-1) )].xyz;

    // gradients
    return (0.5f*(p-m));
}

//mat3 vorticity_aux(ivec3 gridPos)
//{
//    // x+ node
//    dmat3 p;
//    vec3 normalizedGridPos = (gridPos + 0.5f)/vec3(NX,NY,NZ);
//    vec3 scale = 0.5f/vec3(NX,NY,NZ);
//    p[0] = texture(uvwr_old, (normalizedGridPos - vec3(-scale.x,0,0))).xyz - texture(uvwr_old, (normalizedGridPos - vec3(scale.x,0,0))).xyz; // dudn
//    p[1] = texture(uvwr_old, (normalizedGridPos - vec3(0,-scale.y,0))).xyz - texture(uvwr_old, (normalizedGridPos - vec3(0,scale.y,0))).xyz; // dvdn
//    p[2] = texture(uvwr_old, (normalizedGridPos - vec3(0,0,-scale.z))).xyz - texture(uvwr_old, (normalizedGridPos - vec3(0,0,scale.z))).xyz; // dwdn
//
//    // gradients, 1/0.02 = 500
//    return 0.5f*mat3(p[0]/scale,p[1]/scale,p[2]/scale);
//}

void getEigenValuesVectors ( in mat3 mat_data, out mat3 vectors, out vec3 values );
float lambda2()
{
    mat3 J = vorticity();
    mat3 S = J + transpose(J), Q = J - transpose(J);
    mat3 S2Q2 = mat3(
                S[0]*S[0] + Q[0]*Q[0],
            S[1]*S[1] + Q[1]*Q[1],
            S[2]*S[2] + Q[2]*Q[2]);

    mat3 eigenVector;
    vec3 eigenValue;
    getEigenValuesVectors(S2Q2, eigenVector, eigenValue);

    return min(max(eigenValue.x, eigenValue.y), eigenValue.z);
}
uniform float tau_gui;
void main()
{

    // get index in global work group i.e x,y position
    gridPos = ivec3(gl_GlobalInvocationID.xyz) + ntasks.xyz*ntasks.w;
    index = ind(ivec3(0));
    if(gridPos.x >= NX || gridPos.y >= NY || gridPos.z >= NZ) return;

    // Corner stabilizer
    float tau = tau;
    tau = tau_gui+0.5f;
    corner_stabilizer(tau);

    // Streaming
    mat3 f;
    float rho,u,v,w;
    streaming(f,rho,u,v,w);

    // Fetch boundary info
    //vec4 boundary = texelFetch(occl_old, gridPos, 0);
    vec4 boundary = _occl[index];
    float alpha = 2.0f*boundary.z; // will move to a new attachment later
    //boundary = setupBoundary(u,v,rho); // if dynamic

    // no-slip boundary
    //float occlusion = boundary.w;
    //float robin_constant = boundary.x;
    if(boundary.w > 0.0f) // boundary node
    {
        // Initial equilibrium relaxation
        float scale= min(1.0f, float(nstep)/max(256,max(NZ,max(NX,NY))));

        // Dirchlet
        if(boundary.x > 0.0f)
        {
            mat3 cidotu = cx*u*scale + cy*v*scale;
                        f -= mulMat3(2.0f*1.0f*csSqrInv*cidotu, wi);
                        bounce_back(f);
        }
        // Neumann
        if(boundary.x == 0.0f)
        {

            // interior node
            mat3 fi;
            float rhoi,ui,vi,wi;

            // Here is the tricky part: using outlet boundary at top & bottom may lead to reverse flow
            // which will cause instabilities. using value@ivec2(0,2) improves stability.
            ivec3 offset = ivec3(0,0,0);
            if(gridPos.y >(NY-2))  offset = ivec3( 0,-1, 0); // Top
            if(gridPos.y <1)       offset = ivec3( 0, 1, 0); // Bottom
            if(gridPos.x <1)       offset = ivec3( 1, 0, 0);// Left
            if(gridPos.x >(NX-2))  offset = ivec3(-1, 0, 0); // Right
            fi = fetchLattice(rhoi,ui,vi,wi,offset); // Right*/

            // Improved stability when reverse flow may appear
            //if(gridPos.y <2)       {vi = min(0.0, vi); }// Bottom
            //if(gridPos.y >(NY-3))  {vi = max(0.0, vi); } // Top
            //if(gridPos.x <2)       {ui = min(0.0, ui); }// Left
            //if(gridPos.x >(NX-3))  {ui = max(0.0, ui); } // Right

            // RBC
            //f = compute_equilibrium(rhoi, ui ,vi, wi);

            // CBC(2D)
            mat3 X = mat3(u,rho,0,csSqr/rho,u,0,0,0,u);
            mat3 Y = mat3(v,0,rho,0,v,0,csSqr/rho,0,v);

            float cs = sqrt(csSqr);
            float P1 = 0.5f*csSqrInv;
            float P2 = 0.5f/rho/cs;
            //mat3 Pxinv= mat3(P1,0,P1,
            //                 -P2,0,P2,
            //                 0,1,0);

            // interior* <- right node
            float rhoi_1,ui_1,vi_1,wi_1;
            fetchLatticev(rhoi_1,ui_1,vi_1,wi_1,ivec3(-2,0,0));

            // y+ node
            float rho_yp,u_yp,v_yp,w_yp;
            fetchLatticev(rho_yp,u_yp,v_yp,w_yp,ivec3(0,1,0));
            // y- node
            float rho_ym,u_ym,v_ym,w_ym;
            fetchLatticev(rho_ym,u_ym,v_ym,w_ym,ivec3(0,-1,0));

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
            //w = wi;

            f = compute_equilibrium(rho,u*scale,v*scale, 0.0f);

        }
    }
    else // interior node -> collision
    {
        // Compute rho, u and v
        rho = sumMat3(f);
        u = weightedSumMat3(cx, f)/rho;
        v = weightedSumMat3(cy, f)/rho;
        w = 0.0f;

        // External force
        //float Fx = 8f*csSqr*(tau-0.5f)*u_max/pow(float(NY-2), 2f);
        /*float Fx = 0;
        float Fy = 0;

        Fx = weightedSummat3(mulmat3(wi, cx)*Fx*csSqrInv, cx);
        Fy = weightedSummat3(mulmat3(wi, cy)*Fy*csSqrInv, cy);

        u += 0.5f*Fx/rho;
        v += 0.5f*Fy/rho;

        float ueq = u + (tau-0.5f)*Fx/rho;
        float veq = v + (tau-0.5f)*Fy/rho;*/

        //float ueq = u, veq = v, weq = w;

        // Compute equilibrium
        mat3 feq = compute_equilibrium(rho, u, v, w);

        // Entropic Correction
        float beta = 0.5f/(tau);
        alpha = constrain_entropy(f, feq, beta, alpha);
        //alpha = 2.0f;

        // Collision
        f = collision(f, feq, alpha, beta);

        // Alpha relaxation if needed
        boundary.z = min(1.0f,0.5f*alpha+0.005f); // update alpha value with relaxation as initial guess for next time step
        //boundary.z = 2.0f;
    }

    // Output
    //imageStore(uvwr, gridPos, vec4(u,v,w,rho));
    //imageStore(occl, gridPos, boundary);
    uvwr[index] = vec4(u,v,w,rho);
    occl[index] = boundary;
    for(int i = 0; i < 3; i++)
         for(int j = 0; j < 3; j++)
    { _fw0[ind(ivec3(vec3(cx[i][j],cy[i][j],0)))][i][j] = f[i][j]; }

    // Post-processing
    vec3 color;
    float isoSurf;
    color = texture(colormap, sqrt((u)*(u)+v*v)/u_max/2).xyz;
    //color = texture(colormap, sqrt((u-u_max)*(u-u_max)+v*v+w*w)/u_max).xyz;
    //color = texture(colormap, v/u_max+0.5f).xyz;
    //color = texture(colormap, abs(rho-1.0f)*10).xyz;
    //color = vec3(0.1*u*u,v*v,w*w)/u_max/u_max/0.05+0.0;
    //color = texture(colormap, -0.2f*lambda2(gridPos) + 0.1f).xyz;

    isoSurf = boundary.w;
    //isoSurf = sqrt(u*u+v*v)/u_max;
    //isoSurf = lambda2(gridPos)*1.0f+0.5f;
    //if(boundary.w > 0 && corner_detector(gridPos) == 0) {isoSurf = -100000.0f;color = vec3(0.7f);}

    imageStore(vis, gridPos, vec4(color, isoSurf));

    // if dynamic boundaries, setup new boundary
    //vec4 dynbc = setupBoundary(u,v,w,rho);

}

mat3 compute_equilibrium(const float rho, const float u, const float v, const float w)
{
    mat3 cidotu = cx*u + cy*v;
    mat3 feq = rho*(1.0f+3.0f*cidotu+4.5f*mulMat3(cidotu, cidotu)-1.5f*(u*u+v*v));
    return mulMat3(feq, wi);
}

mat3 collision(const mat3 f, const mat3 feq, const float alpha, const float beta)
{
    mat3 df = f - feq;
    float ab = alpha*beta;
    return (f - ab*df);
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

                //if(c[i][j] < 0.0f) {c[i][j] = 0.00001f;}
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
                //if(c[i][j] < 0.0f) {c[i][j] = 0.00001f;}
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

void getEigenValuesVectors ( in mat3 mat_data, out mat3 vectors, out vec3 values )
{
    vec3 e = vec3(0);

    int dimension = 3;
    int dimensionMinusOne = 2;

    for( int j = 0; j < dimension; ++j )
        values[ j ] =  mat_data[dimensionMinusOne][ j ];

    // Householder reduction to tridiagonal form.
    for( int i = dimensionMinusOne; i > 0 && i <= dimensionMinusOne; --i )
    {
        // Scale to avoid under/overflow.
        float scale = 0.0;
        float h =  0.0;
        for( int k = 0; k < i; ++k )
        {
            scale += abs( values[ k ] );
        }

        if( scale ==  0.0 )
        {
            e[ i ] = values[ i - 1 ];

            for( int j = 0; j < i; ++j )
            {
                values[ j ] = mat_data[ ( i - 1 ) ] [  j  ];
                mat_data[ i ][ j ] = 0.0;
                mat_data[ j ][ i ] = 0.0;
            }
        }
        else
        {
            // Generate Householder vector.
            for ( int k = 0; k < i; ++k )
            {
                values[ k ] /= scale;
                h += values[ k ] * values[ k ];
            }

            float f = values[ i - 1 ];
            float g = sqrt( h );

            if ( f >  0.0 )
            {
                g = -g;
            }

            e[ i ] = scale * g;
            h -= f * g;
            values[ i - 1 ] = f - g;

            for ( int j = 0; j < i; ++j)
            {
                e[ j ] =  0.0;
            }

            // Apply similarity transformation to remaining columns.
            for ( int j = 0; j < i; ++j )
            {
                f = values[ j ];
                mat_data[ j ][ i ] = f;
                g = e[ j ] +  mat_data[ j ][ j ] * f;

                for ( int k = j + 1; k <= i - 1; ++k )
                {
                    g +=  mat_data[ k ][ j ] * values[ k ];
                    e[ k ] +=  mat_data[ k ][ j ] * f;
                }

                e[ j ] = g;
            }

            f = 0.0;
            for ( int j = 0; j < i; ++j )
            {
                e[ j ] /= h;
                f += e[ j ] * values[ j ];
            }

            float hh = f / ( h + h );

            for ( int j = 0; j < i; ++j )
            {
                e[ j ] -= hh * values[ j ];
            }

            for ( int j = 0; j < i; ++j )
            {
                f = values[ j ];
                g = e[ j ];

                for ( int k = j; k <= i - 1; ++k )
                {
                    mat_data[ k ][ j ] =   mat_data[ k ][ j ] - ( f * e[ k ] + g * values[ k ] );
                }

                values[ j ] =  mat_data[ i - 1][ j ];
                mat_data[ i ][ j ] = 0.0;
            }
        }
        values[ i ] = h;
    }

    // Accumulate transformations.
    for ( int i = 0; i < dimensionMinusOne; ++i )
    {
        mat_data[dimensionMinusOne][ i ] =  mat_data[ i ][ i ];
        mat_data[ i ][ i ] = 1.0;
        float h = values[ i + 1 ];

        if ( h != 0.0 )
        {
            for ( int k = 0; k <= i; ++k )
            {
                values[ k ] =  mat_data[ k ][ i + 1 ] / h;
            }

            for ( int j = 0; j <= i; ++j )
            {
                float g = 0.0;

                for ( int k = 0; k <= i; ++k )
                {
                    g +=  mat_data[ k ][ i + 1 ] *  mat_data[ k ][ j ];
                }

                for ( int k = 0; k <= i; ++k )
                {
                    mat_data[ k ][ j ] =   mat_data[k][ j ] - ( g * values[ k ] );
                }
            }
        }
        for ( int k = 0; k <= i; ++k )
        {
            mat_data[ k ][ i + 1 ] =  0.0;
        }
    }

    for ( int j = 0; j < dimension; ++j )
    {
        values[ j ] =  mat_data[ dimensionMinusOne ][ j ];
        mat_data[ dimensionMinusOne ][ j ] = 0.0;
    }

    mat_data[ dimensionMinusOne ][ dimensionMinusOne ] =  1.0;
    e[ 0 ] =  0.0;

    for ( int i = 1; i < dimension; ++i )
        e[ i - 1 ] = e[ i ];

    e[ dimensionMinusOne ] = 0.0;

    float f = float( 0.0 );
    float tst1 = float( 0.0 );
    float eps = float( pow( 2.0, -52.0 ));
    for( int l = 0; l < dimension; ++l )
    {
        // Find small subdiagonal element
        tst1 = float( max( tst1, abs ( values[ l ] ) + abs( e[ l ] )));
        int m = l;
        while ( m < dimension )
        {
            if ( abs ( e[ m ] ) <= eps * tst1 ) break;
            ++m;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.
        if( m > l && l<2 )
        {
            int iter = 0;
            do
            {
                ++iter;  // (Could check iteration count here.)
                // Compute implicit shift
                float g = values[ l ];
                float p = ( values[ l + 1 ] - g ) / ( float( 2.0 ) * e[ l ] );
                float r = float( sqrt ( p * p + float( 1.0 ) * float( 1.0 )));
                if( p < 0 ) r = -r;
                values[ l ] = e[ l ] / ( p + r );
                values[ l + 1 ] = e[ l ] * ( p + r );
                float dl1 = values[ l + 1 ];
                float h = g - values[ l ];
                for( int i = l + 2; i < dimension; ++i )
                    values[ i ] -= h;
                f = f + h;

                // Implicit QL transformation.
                p = values[ m ];
                float c = float( 1.0 );
                float c2 = c;
                float c3 = c;
                float el1 = e[ l + 1 ];
                float s = float( 0.0 );
                float s2 = float( 0.0 );
                for ( int i = m - 1; i >= l && i <= m - 1; --i )
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[ i ];
                    h = c * p;
                    r = float( sqrt ( p * p + e[ i ] * e[ i ] ));
                    e[ i + 1 ] = s * r;
                    s = e[ i ] / r;
                    c = p / r;
                    p = c * values[ i ] - s * g;
                    values[ i + 1 ] = h + s * ( c * g + s * values[ i ] );

                    // Accumulate transformation.
                    for( int k = 0; k < dimension; ++k )
                    {
                        h =  mat_data[ k ][ i + 1 ];
                        mat_data[ k ][ i + 1 ] =  ( s *  mat_data[ k ][ i ] + c * h );
                        mat_data[ k ][ i ] = ( c *  mat_data[ k ][ i ] - s * h );
                    }
                }

                p = - s * s2 * c3 * el1 * e[ l ] / dl1;
                e[ l ] = s * p;
                values[ l ] = c * p;
                // Check for convergence.
            }
            while ( abs ( e[ l ] ) > eps * tst1 && iter < 30);
        }
        values[ l ] = values[ l ] + f;
        e[ l ] = float( 0.0 );
    }

    // Sort eigenvalues and corresponding vectors.
    for ( int i = 0; i < dimensionMinusOne; ++i )
    {
        int k = i;
        float p = values[ i ];

        for ( int j = i + 1; j < dimension; ++j )
        {
            if ( values[ j ] < p )
            {
                k = j;
                p = values[ j ];
            }
        }
        if ( k != i )
        {
            values[ k ] = values[ i ];
            values[ i ] = p;
            for ( int j = 0; j < dimension; ++j )
            {
                p =  mat_data[ j ][ i ];
                mat_data[ j ][ i ] =  mat_data[ j ][ k ];
                mat_data[ j ][ k ] = p;
            }
        }
    }

    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            vectors[i][j] = mat_data[i][j];
}
