/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the documentation of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "openglwindow.h"

#include <QtGui/QGuiApplication>
#include <QtGui/QMatrix4x4>
#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QScreen>

#include <QtCore/qmath.h>

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>

#include <iostream>
#include <string>

#define CL_GL_INTEROP
#include "syclode45.hpp"


//! Do not using any namespace or it will cause including issue

//! [1]
class TriangleWindow : public OpenGLWindow
{
public:
    TriangleWindow();

    void initialize() override;
    void render() override;

private:
    std::shared_ptr<ode45<float, 4>> obj;
    GLuint m_posAttr;
    GLuint m_colAttr;
    GLuint m_matrixUniform;

    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vbo;
    GLuint m_count = 1<<18;
    GLuint m_size = 4;
    GLuint m_stride = 4;

    QOpenGLShaderProgram *m_program;
    int m_frame;

};

TriangleWindow::TriangleWindow()
    : m_program(0)
    , m_frame(0)
{
}
//! [1]

//! [2]
int main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    QSurfaceFormat format;
    format.setSamples(16);

    TriangleWindow window;
    window.setFormat(format);
    window.resize(640, 480);
    window.show();

    window.setAnimating(true);

    return app.exec();
}
//! [2]


//! [3]
static const char *vertexShaderSource =
        "attribute highp vec4 posAttr;\n"
        //        "attribute lowp vec4 colAttr;\n"
        //        "varying lowp vec4 col;\n"
        "uniform highp mat4 matrix;\n"
        "void main() {\n"
        //        "   col = colAttr;\n"
        "   gl_Position = matrix * vec4(vec2(posAttr), 0.0f,1.0f);\n"
        //"   gl_PointSize = 1.0f/( 0.1f + 0.055f*gl_Position.z + 0.035f*gl_Position.z*gl_Position.z );\n"
        "}\n";

static const char *fragmentShaderSource =
        //        "varying lowp vec4 col;\n"
        "void main() {\n"
        "   gl_FragColor = vec4(vec3(0.0f), 1.0f);\n"
        "}\n";
//! [3]

//! [4]
void TriangleWindow::initialize()
{
    glEnable(GL_POINT_SPRITE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    m_program = new QOpenGLShaderProgram(this);
    m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
    m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
    m_program->link();
    m_posAttr = m_program->attributeLocation("posAttr"); // should be equal to
    m_colAttr = m_program->attributeLocation("colAttr");
    m_matrixUniform = m_program->uniformLocation("matrix");

    m_camera.setToIdentity();
    m_camera.translate(0, 0, -1);

    // gl buffer module
    m_vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

    // Setup our vertex buffer object.
    m_vbo.create();
    m_vbo.bind();
    m_vbo.allocate(nullptr, m_count * m_size * sizeof(GLfloat));

    // Store the vertex attribute bindings for the program.
    m_vbo.bind();
    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
    f->glEnableVertexAttribArray(0);
    //f->glEnableVertexAttribArray(1);
    f->glVertexAttribPointer(0, m_stride, GL_FLOAT, GL_FALSE, m_size * sizeof(GLfloat), 0);
    //f->glVertexAttribPointer(1, m_stride, GL_FLOAT, GL_FALSE, m_size * sizeof(GLfloat), reinterpret_cast<void *>(m_stride * sizeof(GLfloat)));

    auto hptr = reinterpret_cast<float*>(m_vbo.map(QOpenGLBuffer::ReadWrite)); // cast back to STL type


    //>>>>>>>>>>>>> SYCL module Start
    obj.reset(new ode45<float, 4>);
#ifdef CL_GL_INTEROP
    //glFinish();
    obj->initialize(hptr, m_vbo.bufferId(), std::to_string(m_count).c_str());
#else
    obj->initialize(std::to_string(m_count).c_str());
#endif
    //>>>>>>>>>>>>> SYCL module End


    m_vbo.unmap();
    m_vbo.release();
}
//! [4]

//! [5]
void TriangleWindow::render()
{
    //>>>>>>>>>>>>> SYCL module Start
#ifdef CL_GL_INTEROP
    glFinish();
    obj->calc();
    obj->finalize();
#else
    m_vbo.bind();
    auto hptr = reinterpret_cast<float*>(m_vbo.map(QOpenGLBuffer::ReadWrite)); // cast back to STL type
    obj->calc();
    obj->finalize(hptr, m_size);
    m_vbo.unmap();
#endif
    //>>>>>>>>>>>>> SYCL module End

    // Draw
    const qreal retinaScale = devicePixelRatio();
    glViewport(0, 0, width() * retinaScale, height() * retinaScale);
    resizeGL(width(), height());

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    m_program->bind();

    m_world.setToIdentity();
    m_world.rotate(180.0f - (m_xRot / 16.0f), 1, 0, 0);
    m_world.rotate(m_yRot / 16.0f, 0, 1, 0);
    m_world.rotate(m_zRot / 16.0f, 0, 0, 1);
    m_world.scale(m_zoom);

    m_program->setUniformValue(m_matrixUniform, m_proj * m_camera * m_world);

    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    glDrawArrays(GL_POINTS, 0, m_count);

    m_program->release();

    ++m_frame;
}
//! [5]
