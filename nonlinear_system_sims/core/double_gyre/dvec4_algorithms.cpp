#include "cartesian_field.h"
#include<cstdlib>
#include "shader.h"
#include "cmake_source_dir.h"

namespace algo
{

    static Shader functor[4];

    void reload_subroutines()
    {
        functor[0].reload_shader_program_from_files(FP("algo/double_gyre.glsl"));
        functor[1].reload_shader_program_from_files(FP("algo/double_gyre3d.glsl"));
        functor[2].reload_shader_program_from_files(FP("algo/ks.glsl"));
        functor[3].reload_shader_program_from_files(FP("algo/inertial.glsl"));
    };

    void ode45(simType T, float t0, float t_end, const Vec4& y, const Vec4& dydx)
    {
        // Bind arrays
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,y.d_data); // rhs
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,dydx.d_data); // rhs

        // Deploy kernel
        functor[T].use();

        functor[T].setFloat("t0",t0);
        functor[T].setFloat("t_end",t_end);
        functor[T].setFloat("tol",0.001f);
        functor[T].setFloat("tau_inv",1.0f/0.1f);
        functor[T].setVec4("g",glm::vec4(0,-9.81f,0,0));

        glDispatchCompute(y.mesh->grid<0>(), y.mesh->grid<1>(), y.mesh->grid<2>());

        // make sure writing to image has finished before read
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(0);

        return;
    }

    void kinetic_simulation_configuration(uint& ubo)
    {
        // Magic numbers
        float pi = 3.141592654f;
        float L = 10.0;
        float eta = L/91.0f; // 991: unstable when t > 5
        const uint nk = 64; // explicitly defined in kernels
        float E0 = 3.2f;
        float lambda = 0.5;

        //float LI = L/3.0f;
        //float TL = 0.2f*LI/sqrtf(E0);

        srand(777777);

        std::vector<float> theta,psi,alpha,k,E,dk,a,omega;
        for(uint i = 0; i < nk; i++){ theta.push_back(2*pi*((float)rand()/(float)RAND_MAX)); }
        for(uint i = 0; i < nk; i++){ psi.push_back(2*pi*((float)rand()/(float)RAND_MAX)); }
        // if 2D, let nz = 1 and psi = 0
        //for(uint i = 0; i < nk; i++){ psi.push_back(0); }

        float k1 = 2*pi/L;
        //float k_Nk = 2*pi/eta;

        //float alpha = (L/eta)^(1/(nk-1));
        for(uint i = 0; i < nk; i++){ alpha.push_back(powf(L/eta,1.0f/((float)nk-1))); }
        //float k = k1*alpha.^((1:nk).'-1);
        for(uint i = 0; i < nk; i++){ k.push_back(k1*powf(alpha[i],(float)i)); }
        //float E = E0*L*(k*L).^(-5/3);
        for(uint i = 0; i < nk; i++){ E.push_back(E0*L*powf((k[i]*L),-5.0f/3.0f)); }
        //float dk = gradient(k);
        for(uint i = 0; i < nk; i++){
            if(i == 0){ dk.push_back(k[1]-k[0]);continue; }
            if(i == nk-1){ dk.push_back(k[nk-1]-k[nk-2]);break; }
            dk.push_back(0.5f*(k[i+1]-k[i-1]));
        }

        //float a = sqrt(E.*dk);
        for(uint i = 0; i < nk; i++){ a.push_back(sqrtf(E[i]*dk[i])); }
        //float b = a;

        //float omega = lambda*sqrt(k.^3.*E).';
        for(uint i = 0; i < nk; i++)
        {
            omega.push_back( lambda*sqrtf(E[i]*powf(k[i],3.0f)) );
        }

        // Vectors
        //An = (a.*[ cos(theta) -sin(theta)]).';
        //Bn = (b.*[-cos(theta)  sin(theta)]).';
        //Kn = (k.*[ sin(theta)  cos(theta)]).';
        std::vector<float> a_xt, a_yt, a_zt, k_xt, k_yt, k_zt;
        std::vector<float> b_xt, b_yt, b_zt;
        for(uint i = 0; i < nk; i++){
            a_xt.push_back( a[i]*cosf(theta[i])*cosf(psi[i]) );
            a_yt.push_back(-a[i]*sinf(theta[i]) );
            a_zt.push_back( a[i]*cosf(theta[i])*sinf(psi[i]) );

            b_xt.push_back(-a[i]*cosf(theta[i])*sinf(psi[i]) ); // -a_zt
            b_yt.push_back( a[i]*sinf(theta[i]) ); // -a_yt
            b_zt.push_back(-a[i]*cosf(theta[i])*cosf(psi[i]) ); // -a_xt

            k_xt.push_back( k[i]*sinf(theta[i])*cosf(psi[i]) );
            k_yt.push_back( k[i]*cosf(theta[i]) );
            k_zt.push_back( k[i]*sinf(theta[i])*sinf(psi[i]) );
        }

        glBindBuffer(GL_UNIFORM_BUFFER, ubo);
        glBufferData(GL_UNIFORM_BUFFER, 7*nk*sizeof(float) , NULL, GL_STATIC_DRAW);
        glBufferSubData(GL_UNIFORM_BUFFER, 0*nk*sizeof(float), nk*sizeof(float), &a_xt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 1*nk*sizeof(float), nk*sizeof(float), &a_yt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 2*nk*sizeof(float), nk*sizeof(float), &a_zt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 3*nk*sizeof(float), nk*sizeof(float), &k_xt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 4*nk*sizeof(float), nk*sizeof(float), &k_yt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 5*nk*sizeof(float), nk*sizeof(float), &k_zt[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 6*nk*sizeof(float), nk*sizeof(float), &omega[0]);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        // binding to UBO0
        glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo); // Bind to UBO0
        return;
    }

}


