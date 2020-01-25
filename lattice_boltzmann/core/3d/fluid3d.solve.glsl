#version 430
// Kernel
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

// Visualization
layout(binding = 15) uniform sampler1D colormap;

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
    int nstep;
};

// Functions
float sumL3(DiscreteLattice f)
{
    float result=0;
    for(int i = 0; i < 7; i++)
    {
        result += f.v[i].x;
        result += f.v[i].y;
        result += f.v[i].z;
        result += f.v[i].w;
    }
    return result;
}

float sumL3w(DiscreteLattice w, DiscreteLattice f)
{
    float result=0;
    for(int i = 0; i < 7; i++)
    {
        vec4 _t = w.v[i]*f.v[i];
        result += _t.x;
        result += _t.y;
        result += _t.z;
        result += _t.w;
    }
    return result;
}

DiscreteLattice mulL3(DiscreteLattice A, DiscreteLattice B)
{
    DiscreteLattice result;
    for(int i = 0; i < 7; i++)
    {
        result.v[i] = A.v[i]*B.v[i];
    }
    return result;
}
DiscreteLattice mulL3f(DiscreteLattice A, float b)
{
    DiscreteLattice result;
    for(int i = 0; i < 7; i++)
    {
        result.v[i] = A.v[i]*b;
    }
    return result;
}

DiscreteLattice plsL3(DiscreteLattice A, DiscreteLattice B)
{
    DiscreteLattice result;
    for(int i = 0; i < 7; i++)
    {
        result.v[i] = A.v[i] + B.v[i];
    }
    return result;
}

DiscreteLattice minusL3(DiscreteLattice A, DiscreteLattice B)
{
    DiscreteLattice result;
    for(int i = 0; i < 7; i++)
    {
        result.v[i] = A.v[i] - B.v[i];
    }
    return result;
}

// lattice boltzmann method
//void streaming(inout DiscreteLattice f, inout float r, inout float u, inout float v, inout float w, ivec3 gridPos);
//DiscreteLattice fetchLattice(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos);
//void fetchLatticev(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos);
DiscreteLattice compute_equilibrium( float rho,  float u,  float v,  float w);
DiscreteLattice collision( DiscreteLattice f,  DiscreteLattice feq,  float alpha,  float beta);

// entropy constraint
float compute_g( DiscreteLattice f, DiscreteLattice feq, const float a, const float b);
float compute_gradg( DiscreteLattice f, DiscreteLattice feq, const float a, const float b);
float constrain_entropy( DiscreteLattice f, DiscreteLattice feq, const float b, const float aold);

// boundary
//vec4 setupBoundary(inout float uw,inout float vw,inout float rhow);

ivec3 relocUV(ivec3 v)
{
    return ivec3(v.x - ((v.x + NX)/NX - 1)*NX,
                 v.y - ((v.y + NY)/NY - 1)*NY,
                 v.z - ((v.z + NZ)/NZ - 1)*NZ);
}

// set a bi later...
//const int bi[28] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25,27}; // D3Q27
//const int bi[28] = {1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,26,27}; // D3Q27
//const int bi[20] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,19}; // D3Q19
void bounce_back(inout DiscreteLattice f)
{

    DiscreteLattice result;
    for(int i = 0; i < 6;i++)
    {
        result.v[i].xyzw = f.v[i].yxwz;
    }
    result.v[6].xyzw = f.v[6].yxzw;
    f = result;
}
float one_sided_diff(float x1, float x2, float x3, float sign)
{
    return sign*0.5f*(-3.0f*x1+4.0f*x2-x3);
}
float central_diff(float x1, float x2)
{
    return 0.5f*(x1-x2);
}
int ind(ivec3 gridPos)
{
    return (gridPos.x + gridPos.y*NX + gridPos.z*NX*NY);
}

void streaming(inout DiscreteLattice f, inout float r, inout float u, inout float v, inout float w, ivec3 gridPos)
{
    // Streaming
    int index = ind(gridPos);
    f.v[0] = _fr0[index][0];
    f.v[1] = _fr0[index][1];
    f.v[2] = _fr1[index][0];
    f.v[3] = _fr1[index][1];
    f.v[4] = _fr2[index][0];
    f.v[5] = _fr2[index][1];
    f.v[6] = _fr3[index];

    //vec4 a = texelFetch(uvwr_old, gridPos, 0);
    vec4 a = _uvwr[index];

    u = a.x;
    v = a.y;
    w = a.z;
    r = a.w;

}

DiscreteLattice fetchLattice(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos)
{
    DiscreteLattice f;
    for(int j = 0; j < 4; j++){int i = 0;f.v[i][j] = _fr0[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j];}
    for(int j = 0; j < 4; j++){int i = 1;f.v[i][j] = _fr0[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j];}
    for(int j = 0; j < 4; j++){int i = 2;f.v[i][j] = _fr1[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j];}
    for(int j = 0; j < 4; j++){int i = 3;f.v[i][j] = _fr1[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j];}
    for(int j = 0; j < 4; j++){int i = 4;f.v[i][j] = _fr2[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j];}
    for(int j = 0; j < 4; j++){int i = 5;f.v[i][j] = _fr2[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j];}
    for(int j = 0; j < 4; j++){int i = 6;f.v[i][j] = _fr3[ind(relocUV(gridPos - dir + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))]   [j];}

    //vec4 a = texelFetch(uvwr_old, relocUV(gridPos - dir), 0);
    vec4 a = _uvwr[ind(relocUV(gridPos - dir))];
    u = a.x;
    v = a.y;
    w = a.z;
    r = a.w;

    return f;
}

void fetchLatticev(inout float r, inout float u, inout float v, inout float w, ivec3 dir, ivec3 gridPos)
{

    //vec4 a = texelFetch(uvwr_old, relocUV(gridPos - dir), 0);
    vec4 a = _uvwr[ind(relocUV(gridPos - dir))];
    u = a.x;
    v = a.y;
    w = a.z;
    r = a.w;
}

int corner_detector(ivec3 gridPos)
{
    int flag=0,width=1;
    // Artificially increase tau
    if(gridPos.y <2+width || gridPos.y >(NY-3-width))       flag += 1;// Next to Bottom Face
    if(gridPos.x <2+width || gridPos.x >(NX-3-width))       flag += 1;// Next to Left Face
    if(gridPos.z <2+width || gridPos.z >(NZ-3-width))       flag += 1;// Next to Front Face

    // Corner
    return flag;

}

void corner_stabilizer(inout float t, ivec3 gridPos)
{
    // Corner
    if(corner_detector(gridPos) > 2) t = max(0.55f, t);

    // Edge
    if(corner_detector(gridPos) > 1) t = max(0.505f,t);

}

mat3 vorticity(ivec3 gridPos)
{
    // x+ node
    mat3 p,m;
    //p[0] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(-1,0,0)), 0).xyz;
    //p[1] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(0,-1,0)), 0).xyz;
    //p[2] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(0,0,-1)), 0).xyz;
    //m[0] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(1,0,0)), 0).xyz;
    //m[1] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(0,1,0)), 0).xyz;
    //m[2] = texelFetch(uvwr_old, relocUV(gridPos - ivec3(0,0,1)), 0).xyz;
    p[0] = _uvwr[ind( relocUV(gridPos - ivec3(-1,0,0)))].xyz;
    p[1] = _uvwr[ind( relocUV(gridPos - ivec3(0,-1,0)))].xyz;
    p[2] = _uvwr[ind( relocUV(gridPos - ivec3(0,0,-1)))].xyz;
    m[0] = _uvwr[ind( relocUV(gridPos - ivec3(1,0,0) ))].xyz;
    m[1] = _uvwr[ind( relocUV(gridPos - ivec3(0,1,0) ))].xyz;
    m[2] = _uvwr[ind( relocUV(gridPos - ivec3(0,0,1) ))].xyz;

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
float lambda2(ivec3 gridPos)
{
    mat3 J = vorticity(gridPos);
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

void main()
{

    // get index in global work group i.e x,y position
    ivec3 gridPos = ivec3(gl_GlobalInvocationID.xyz) + ntasks.xyz*ntasks.w;
    if(gridPos.x >= NX || gridPos.y >= NY || gridPos.z >= NZ) return;

    // Corner stabilizer
    float tau = tau;
    corner_stabilizer(tau, gridPos);

    // Streaming
    DiscreteLattice f;
    float rho,u,v,w;
    streaming(f,rho,u,v,w, gridPos);

    // Fetch boundary info
    //vec4 boundary = texelFetch(occl_old, gridPos, 0);
    vec4 boundary = _occl[ind(gridPos)];
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
            DiscreteLattice cidotu =
                    plsL3(
                        plsL3(
                            mulL3f(mulL3f(cx,u),scale),
                            mulL3f(mulL3f(cy,v),scale)),
                        mulL3f(mulL3f(cz,w),scale));
            f = minusL3(f, mulL3(mulL3f(cidotu,2.0f*1.0f*csSqrInv), wi));

            bounce_back(f);
        }
        // Neumann
        if(boundary.x == 0.0f)
        {

            // interior node
            DiscreteLattice fi;
            float rhoi,ui,vi,wi;

            // Here is the tricky part: using outlet boundary at top & bottom may lead to reverse flow
            // which will cause instabilities. using value@ivec2(0,2) improves stability.
            ivec3 offset = ivec3(0,0,0);
            if(gridPos.y >(NY-2))  offset = ivec3( 0, 1, 0); // Top
            if(gridPos.y <1)       offset = ivec3( 0,-1, 0); // Bottom
            if(gridPos.x <1)       offset = ivec3(-1, 0, 0);// Left
            if(gridPos.x >(NX-2))  offset = ivec3( 1, 0, 0); // Right
            fi = fetchLattice(rhoi,ui,vi,wi,offset,gridPos); // Right*/

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

            // interior* node
            float rhoi_1,ui_1,vi_1,wi_1;
            fetchLatticev(rhoi_1,ui_1,vi_1,wi_1,ivec3(2,0,0),gridPos);

            // y+ node
            float rho_yp,u_yp,v_yp,w_yp;
            fetchLatticev(rho_yp,u_yp,v_yp,w_yp,ivec3(0,-1,0),gridPos);
            // y- node
            float rho_ym,u_ym,v_ym,w_ym;
            fetchLatticev(rho_ym,u_ym,v_ym,w_ym,ivec3(0,1,0),gridPos);

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

            f = compute_equilibrium(rho,u*scale,v*scale, w*scale);

        }
    }
    else // interior node -> collision
    {
        // Compute rho, u and v
        rho = sumL3(f);
        u = sumL3w(cx, f)/rho;
        v = sumL3w(cy, f)/rho;
        w = sumL3w(cz, f)/rho;

        // External force
        //float Fx = 8f*csSqr*(tau-0.5f)*u_max/pow(float(NY-2), 2f);
        /*float Fx = 0;
        float Fy = 0;
        Fx = weightedSumDiscreteLattice(mulDiscreteLattice(wi, cx)*Fx*csSqrInv, cx);
        Fy = weightedSumDiscreteLattice(mulDiscreteLattice(wi, cy)*Fy*csSqrInv, cy);
        u += 0.5f*Fx/rho;
        v += 0.5f*Fy/rho;
        float ueq = u + (tau-0.5f)*Fx/rho;
        float veq = v + (tau-0.5f)*Fy/rho;*/

        //float ueq = u, veq = v, weq = w;

        // Compute equilibrium
        DiscreteLattice feq = compute_equilibrium(rho, u, v, w);

        // Entropic Correction
        float beta = 0.5f/(tau);
        alpha = constrain_entropy(f, feq, beta, alpha);

        // Collision
        f = collision(f, feq, alpha, beta);

        // Alpha relaxation if needed
        boundary.z = min(1.0f,0.5f*alpha+0.005f); // update alpha value with relaxation as initial guess for next time step
        //boundary.z = 2.0f;
    }

    // Output
    uvwr[ind(gridPos)] = vec4(u,v,w,rho);
    occl[ind(gridPos)] = boundary;
    for(int j = 0; j < 4; j++){int i = 0; _fw0[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 1; _fw0[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 2; _fw1[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 3; _fw1[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 4; _fw2[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][0][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 5; _fw2[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))][1][j] = f.v[i][j];}
    for(int j = 0; j < 4; j++){int i = 6; _fw3[ind(relocUV(gridPos + ivec3(cx.v[i][j],cy.v[i][j],cz.v[i][j])))]   [j] = f.v[i][j];}

    // Post-processing
    vec3 color;
    float isoSurf;
    color = texture(colormap, sqrt((u)*(u)+v*v)/u_max/2).xyz;
    //color = texture(colormap, sqrt((u-u_max)*(u-u_max)+v*v+w*w)/u_max).xyz;
    //color = texture(colormap, v/u_max+0.5f).xyz;
    //color = texture(colormap, abs(rho-1.0f)*10).xyz;
    //color = vec3(0.1*u*u,v*v,w*w)/u_max/u_max/0.05+0.0;
    //color = texture(colormap, -0.2f*lambda2(gridPos) + 0.1f).xyz;
    //color.x = sqrt((u)*(u)+v*v)/u_max;
    //color.y = alpha;
    isoSurf = boundary.w;
    //isoSurf = sqrt(u*u+v*v)/u_max;
    //isoSurf = lambda2(gridPos)*1.0f+0.5f;
    //if(boundary.w > 0 && corner_detector(gridPos) == 0) {isoSurf = -100000.0f;color = vec3(0.7f);}

    imageStore(vis, gridPos, vec4(color, isoSurf));

    // if dynamic boundaries, setup new boundary
    //vec4 dynbc = setupBoundary(u,v,w,rho);

}

DiscreteLattice compute_equilibrium(float rho, float u, float v, float w)
{
    DiscreteLattice workMatrix;
    for(int i = 0; i < 7; i++)
    {
        vec4 cidotu = cx.v[i]*u + cy.v[i]*v + cz.v[i]*w;
        workMatrix.v[i] = wi.v[i]*rho*(1.0f+3.0f*cidotu+4.5f*cidotu*cidotu-1.5f*(u*u+v*v+w*w));
    }

    return workMatrix;
}

DiscreteLattice collision(DiscreteLattice f, DiscreteLattice feq, float alpha, float beta)
{
    DiscreteLattice df = minusL3(f, feq);
    float ab = alpha*beta;
    return minusL3(f, mulL3f(df, ab));
}


float compute_g(
        DiscreteLattice f,
        DiscreteLattice feq,
        const float a,const float b) {

    float result = 0.0f;
    DiscreteLattice c = collision(f,feq,a,b);

    for(int i = 0;i < 7; i++)
    {
        for(int j = 0;j < 4; j++)
        {
            if((4*i+j) > 26) continue;
            if(c.v[i][j] < 0.0f) {c.v[i][j] = 0.00001f;}
            result += c.v[i][j]*log(c.v[i][j]/wi.v[i][j]) - f.v[i][j]*log(f.v[i][j]/wi.v[i][j]);
        }
    }
    return result;
}

float compute_gradg(
        DiscreteLattice f,
        DiscreteLattice feq,
        const float a, const float b)
{

    float result = 0.0f;
    DiscreteLattice c = collision(f,feq,a,b);

    for(int i = 0;i < 7; i++)
    {   for(int j = 0;j < 4; j++)
        {
            if((4*i+j) > 26) continue;
            if(c.v[i][j] < 0.0f) {c.v[i][j] = 0.00001f;}
            result += -b*(f.v[i][j] - feq.v[i][j])*(log(c.v[i][j]/wi.v[i][j]) + 1.0f);
        }
    }
    return result;
}

float constrain_entropy(
        DiscreteLattice f,
        DiscreteLattice feq,
        const float b,
        const float alpha_old)
{
    // calculate deviation
    float amin=0.2f;
    float amax=alpha_old;
    float maxDeviation = 0.0f;
    for(int i = 0;i < 7; i++)
    {   for(int j = 0;j < 4; j++)
        {
            if((4*i+j) > 26) continue;
            float deviation = abs(f.v[i][j]-feq.v[i][j])/feq.v[i][j];
            if(deviation > maxDeviation)
                maxDeviation = deviation;
        }
    }

    // if deviation is too large
    //float stableDeviation = 0.2;
    if(maxDeviation < 0.2f) return amax;

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

    int maxIter = 10;
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
