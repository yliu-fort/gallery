#include "cartesian_field.h"
#include<cstdlib>
#include "shader.h"
#include "cmake_source_dir.h"

namespace algo
{
static Shader functor[5];

void kinetic_simulation23_inertial(const Vec4&, const Vec4&, const Param&);
void kinetic_simulation23(const Vec4&, const Vec4&, const Param&);
void kinetic_simulation_configuration2D(uint&, const Param&);
void kinetic_simulation_configuration3D(uint&, const Param&);

void reload_subroutines()
{
    //functor[DOUBLE_GYRE].reload_shader_program_from_files(FP("algo/double_gyre.glsl"));
    //functor[DOUBLE_GYRE3D].reload_shader_program_from_files(FP("algo/double_gyre3d.glsl"));
    //functor[KINETIC_TURB2D].reload_shader_program_from_files(FP("algo/ks2D.glsl"));
    //functor[KINETIC_TURB3D].reload_shader_program_from_files(FP("algo/ks.glsl"));
    functor[INERTIAL_PARTICLE].reload_shader_program_from_files(FP("algo/inertial.glsl"));
}

void ode45_init(uint& ubo, const Param& param)
{
    if(param.sim_type == KINETIC_TURB2D)
        kinetic_simulation_configuration2D(ubo, param);
    if(param.sim_type == KINETIC_TURB3D)
        kinetic_simulation_configuration3D(ubo, param);
    if(param.sim_type == INERTIAL_PARTICLE)
        kinetic_simulation_configuration3D(ubo, param);
}

void ode45(const Vec4& y, const Vec4& dydx, const Param& param)
{
    if(param.sim_type == INERTIAL_PARTICLE)
        kinetic_simulation23_inertial(y, dydx, param);
    else
        kinetic_simulation23(y, dydx, param);
}

void kinetic_simulation23_inertial(const Vec4& y, const Vec4& dydx, const Param& param)
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,y.get_dptr()); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,dydx.get_dptr()); // rhs

    // Deploy kernel
    functor[param.sim_type].use();

    functor[param.sim_type].setFloat("t0",param.getCurrentTimestep());
    functor[param.sim_type].setFloat("t_end",param.getNextTimestep());
    functor[param.sim_type].setFloat("tol",param.tol);

    functor[param.sim_type].setFloat("tau_inv",param.tau_inv);
    functor[param.sim_type].setVec4("g",glm::vec4(0,param.g,0,0));

    y.dispatch();

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return;
}

void kinetic_simulation23(const Vec4& y, const Vec4& dydx, const Param& param)
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,y.get_dptr()); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,dydx.get_dptr()); // rhs

    // Deploy kernel
    functor[param.sim_type].use();

    functor[param.sim_type].setFloat("t0",param.getCurrentTimestep());
    functor[param.sim_type].setFloat("t_end",param.getNextTimestep());
    functor[param.sim_type].setFloat("tol",param.tol);

    y.dispatch();

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return;
}

void kinetic_simulation_configuration3D(uint& ubo, const Param& param)
{

    // Magic numbers
    //float pi = 3.141592654f;
    //float L = 10.0;
    //float eta = L/91.0f; // 991: unstable when t > 5
    //const uint nk = 64; // explicitly defined in kernels
    //float E0 = 3.2f;
    //float lambda = 0.5;

    //float LI = L/3.0f;
    //float TL = 0.2f*LI/sqrtf(E0);

    srand(param.randomSeed);

    std::vector<double> theta,psi,alpha,k,E,dk,a;

    for(uint i = 0; i < param.nk; i++){ theta.push_back(2.0*param.pi*(rand()/(double)RAND_MAX)); }
    for(uint i = 0; i < param.nk; i++){ psi.push_back(2.0*param.pi*(rand()/(double)RAND_MAX)); }

    double k1 = 2.0*param.pi/param.L;
    //float k_Nk = 2*pi/eta;

    //float alpha = (L/eta)^(1/(nk-1));
    for(uint i = 0; i < param.nk; i++){ alpha.push_back(pow(param.L/param.eta,1.0/(param.nk-1.0))); }
    //float k = k1*alpha.^((1:nk).'-1);
    for(uint i = 0; i < param.nk; i++){ k.push_back(k1*pow(alpha[i],(double)i)); }
    //float E = E0*L*(k*L).^(-5/3);
    for(uint i = 0; i < param.nk; i++){ E.push_back(param.E0*param.L*pow((k[i]*param.L),param.energy_cascade)); }
    //float dk = gradient(k);
    for(uint i = 0; i < param.nk; i++){
        if(i == 0){ dk.push_back(k[1]-k[0]);continue; }
        if(i == param.nk-1){ dk.push_back(k[param.nk-1]-k[param.nk-2]);break; }
        dk.push_back(0.5*(k[i+1]-k[i-1]));
    }

    //float a = sqrt(E.*dk);
    for(uint i = 0; i < param.nk; i++){ a.push_back(sqrt(E[i]*dk[i])); }
    //float b = a;

    //float omega = lambda*sqrt(k.^3.*E).';
    // Convert from double to float in order to squeeze into the buffer
    std::vector<float> a_xt, a_yt, a_zt, k_xt, k_yt, k_zt, omega;
    for(uint i = 0; i < param.nk; i++)
    {
        omega.push_back( param.lambda*sqrt(E[i]*pow(k[i],3.0)) );
    }

    // Vectors
    //An = (a.*[ cos(theta) -sin(theta)]).';
    //Bn = (b.*[-cos(theta)  sin(theta)]).';
    //Kn = (k.*[ sin(theta)  cos(theta)]).';
    //std::vector<float> b_xt, b_yt, b_zt;
    for(uint i = 0; i < param.nk; i++){
        a_xt.push_back( a[i]*cos(theta[i])*cos(psi[i]) );
        a_yt.push_back(-a[i]*sin(theta[i]) );
        a_zt.push_back( a[i]*cos(theta[i])*sin(psi[i]) );

        k_xt.push_back( k[i]*sin(theta[i])*cos(psi[i]) );
        k_yt.push_back( k[i]*cos(theta[i]) );
        k_zt.push_back( k[i]*sin(theta[i])*sin(psi[i]) );

    }

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, 7*param.nk*sizeof(float) , NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_UNIFORM_BUFFER, 0*param.nk*sizeof(float), param.nk*sizeof(float), &a_xt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 1*param.nk*sizeof(float), param.nk*sizeof(float), &a_yt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 2*param.nk*sizeof(float), param.nk*sizeof(float), &a_zt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 3*param.nk*sizeof(float), param.nk*sizeof(float), &k_xt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 4*param.nk*sizeof(float), param.nk*sizeof(float), &k_yt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 5*param.nk*sizeof(float), param.nk*sizeof(float), &k_zt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 6*param.nk*sizeof(float), param.nk*sizeof(float), &omega[0]);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    // binding to UBO0
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo); // Bind to UBO0
    return;
}

void kinetic_simulation_configuration2D(uint& ubo, const Param& param)
{
    // Magic numbers
    //float pi = 3.141592654f;
    //float L = 10.0;
    //float eta = L/91.0f; // 991: unstable when t > 5
    //const uint nk = 64; // explicitly defined in kernels
    //float E0 = 3.2f;
    //float lambda = 0.5;

    //float LI = L/3.0f;
    //float TL = 0.2f*LI/sqrtf(E0);

    srand(param.randomSeed);

    std::vector<double> theta,psi,alpha,k,E,dk,a;

    for(uint i = 0; i < param.nk; i++){ theta.push_back(2.0*param.pi*(rand()/(double)RAND_MAX)); }
    for(uint i = 0; i < param.nk; i++){ psi.push_back(2.0*param.pi*(rand()/(double)RAND_MAX)); }

    double k1 = 2.0*param.pi/param.L;
    //float k_Nk = 2*pi/eta;

    //float alpha = (L/eta)^(1/(nk-1));
    for(uint i = 0; i < param.nk; i++){ alpha.push_back(pow(param.L/param.eta,1.0/(param.nk-1.0))); }
    //float k = k1*alpha.^((1:nk).'-1);
    for(uint i = 0; i < param.nk; i++){ k.push_back(k1*pow(alpha[i],(double)i)); }
    //float E = E0*L*(k*L).^(-5/3);
    for(uint i = 0; i < param.nk; i++){ E.push_back(param.E0*param.L*pow((k[i]*param.L),param.energy_cascade)); }
    //float dk = gradient(k);
    for(uint i = 0; i < param.nk; i++){
        if(i == 0){ dk.push_back(k[1]-k[0]);continue; }
        if(i == param.nk-1){ dk.push_back(k[param.nk-1]-k[param.nk-2]);break; }
        dk.push_back(0.5*(k[i+1]-k[i-1]));
    }

    //float a = sqrt(E.*dk);
    for(uint i = 0; i < param.nk; i++){ a.push_back(sqrt(E[i]*dk[i])); }
    //float b = a;

    //float omega = lambda*sqrt(k.^3.*E).';
    // Convert from double to float in order to squeeze into the buffer
    std::vector<float> a_xt, a_yt, k_xt, k_yt, omega;
    for(uint i = 0; i < param.nk; i++)
    {
        omega.push_back( param.lambda*sqrt(E[i]*pow(k[i],3.0)) );
    }

    // Vectors
    //An = (a.*[ cos(theta) -sin(theta)]).';
    //Bn = (b.*[-cos(theta)  sin(theta)]).';
    //Kn = (k.*[ sin(theta)  cos(theta)]).';
    //std::vector<float> b_xt, b_yt, b_zt;
    for(uint i = 0; i < param.nk; i++){
        a_xt.push_back( a[i]*cos(theta[i]) );
        a_yt.push_back(-a[i]*sin(theta[i]) );

        k_xt.push_back( k[i]*sin(theta[i]) );
        k_yt.push_back( k[i]*cos(theta[i]) );
    }

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, 5*param.nk*sizeof(float) , NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_UNIFORM_BUFFER, 0*param.nk*sizeof(float), param.nk*sizeof(float), &a_xt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 1*param.nk*sizeof(float), param.nk*sizeof(float), &a_yt[0]);

    glBufferSubData(GL_UNIFORM_BUFFER, 2*param.nk*sizeof(float), param.nk*sizeof(float), &k_xt[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 3*param.nk*sizeof(float), param.nk*sizeof(float), &k_yt[0]);

    glBufferSubData(GL_UNIFORM_BUFFER, 4*param.nk*sizeof(float), param.nk*sizeof(float), &omega[0]);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    // binding to UBO0
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo); // Bind to UBO1
    return;
}


}
