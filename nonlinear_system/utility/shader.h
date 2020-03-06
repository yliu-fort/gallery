#ifndef SHADER_H
#define SHADER_H

#include <GL/glew.h>
#include "glm/glm.hpp"

#include <string>
#include <iostream>
#include <vector>

class Shader
{
public:
    unsigned int ID;
    // constructor generates the shader on the fly
    // ------------------------------------------------------------------------
    Shader():ID(0),isComputeShader(false) {}
    Shader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr);
    Shader(const char* computePath); // Compute shader, require OpenGL version >= 4.0
    void reload_shader_program_from_files(const char*,const char*,const char* = nullptr );
    void reload_shader_program_from_files(const char*);

    // activate the shader
    // ------------------------------------------------------------------------
    void use() const
    {
        glUseProgram(ID);
    }
    // utility uniform functions
    // ------------------------------------------------------------------------
    void setBool(const std::string &name, bool value) const
    {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), static_cast<int>(value));
    }
    // ------------------------------------------------------------------------
    void setInt(const std::string &name, int value) const
    {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
    }
    // ------------------------------------------------------------------------
    void setFloat(const std::string &name, float value) const
    {
        glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
    }
    // ------------------------------------------------------------------------
    void setVec2(const std::string &name, const glm::vec2 &value) const
    {
        glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec2(const std::string &name, float x, float y) const
    {
        glUniform2f(glGetUniformLocation(ID, name.c_str()), x, y);
    }
    // ------------------------------------------------------------------------
    void setVec3(const std::string &name, const glm::vec3 &value) const
    {
        glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec3(const std::string &name, float x, float y, float z) const
    {
        glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);
    }
    void setVec3i(const std::string &name, const glm::ivec3 &value) const
    {
        glUniform3iv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec3i(const std::string &name, int x, int y, int z) const
    {
        glUniform3i(glGetUniformLocation(ID, name.c_str()), x, y, z);
    }
    // ------------------------------------------------------------------------
    void setVec4(const std::string &name, const glm::vec4 &value) const
    {
        glUniform4fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec4(const std::string &name, float x, float y, float z, float w)
    {
        glUniform4f(glGetUniformLocation(ID, name.c_str()), x, y, z, w);
    }
    // ------------------------------------------------------------------------
    void setMat2(const std::string &name, const glm::mat2 &mat) const
    {
        glUniformMatrix2fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
    // ------------------------------------------------------------------------
    void setMat3(const std::string &name, const glm::mat3 &mat) const
    {
        glUniformMatrix3fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
    // ------------------------------------------------------------------------
    void setMat4(const std::string &name, const glm::mat4 &mat) const
    {
        glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
    // ------------------------------------------------------------------------
    void setUniformBlockBinding(const char *name, const int &binding) const
    {
        unsigned int index = glGetUniformBlockIndex(ID, name);
        glUniformBlockBinding(ID, index, binding);
    }
    // ------------------------------------------------------------------------
    const std::string& getVertexPath  () const {return _vertexPath;}
    const std::string& getFragmentPath() const {return _fragmentPath;}
    const std::string& getGeometryPath() const {return _geometryPath;}
    const std::string& getKernelPath  () const {return _vertexPath;}
    //void setVertexPath  (const char *name);
    //void setFragmentPath(const char *name);
    //void setGeometryPath(const char *name);
    //void setKernelPath  (const char *name);

private:
    std::string _vertexPath;
    std::string _fragmentPath;
    std::string _geometryPath;
    const bool isComputeShader;
    // utility function for checking shader compilation/linking errors.
    // ------------------------------------------------------------------------
    void checkCompileErrors(GLuint shader, std::string type);
};

/*class ShaderManagement
{
public:
    static void Create(){}
    static void Update(){}
    static void Destroy(){}
private:
    ShaderManagement* m_singleton;
    std::vector<Shader*> m_container;
    ShaderManagement(){}
    ~ShaderManagement(){}
};*/

#endif
