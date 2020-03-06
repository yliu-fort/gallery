#include "cartesian_field.h"
#include<cstdlib>
#include "shader.h"
#include "cmake_source_dir.h"

namespace algo
{
    static Shader functor[4];

    void reload_subroutines()
    {
        functor[0].reload_shader_program_from_files(FP("algo/ode45.glsl"));
        functor[1].reload_shader_program_from_files(FP("algo/ks.glsl"));
    };

    void ode45(float t0, float t_end, const Vec4& y, uint option)
    {
        // Bind arrays
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,y.d_data); // rhs

        // Deploy kernel
        functor[option].use();
        functor[option].setFloat("t0",t0);
        functor[option].setFloat("t_end",t_end);
        functor[option].setFloat("tol",0.001f);
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

        std::vector<float> theta,alpha,k,E,dk,a,omega;
        for(uint i = 0; i < nk; i++){ theta.push_back(2*pi*((float)rand()/(float)RAND_MAX)); }

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
        std::vector<float> a_st, a_ct, k_st, k_ct;
        for(uint i = 0; i < nk; i++){
            a_st.push_back(a[i]*sinf(theta[i]));
            a_ct.push_back(a[i]*cosf(theta[i]));
            k_st.push_back(k[i]*sinf(theta[i]));
            k_ct.push_back(k[i]*cosf(theta[i]));
        }

        glBindBuffer(GL_UNIFORM_BUFFER, ubo);
        glBufferData(GL_UNIFORM_BUFFER, 5*nk*sizeof(float) , NULL, GL_STATIC_DRAW);
        glBufferSubData(GL_UNIFORM_BUFFER, 0*nk*sizeof(float), nk*sizeof(float), &a_st[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 1*nk*sizeof(float), nk*sizeof(float), &a_ct[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 2*nk*sizeof(float), nk*sizeof(float), &k_st[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 3*nk*sizeof(float), nk*sizeof(float), &k_ct[0]);
        glBufferSubData(GL_UNIFORM_BUFFER, 4*nk*sizeof(float), nk*sizeof(float), &omega[0]);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        // binding to UBO0
        glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo); // Bind to UBO0
        return;
    }

}


