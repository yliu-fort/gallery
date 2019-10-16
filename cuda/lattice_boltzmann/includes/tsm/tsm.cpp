#include "tsm.h"
#include "camera.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <GLFW/glfw3.h>

#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace glm;

// Helper functions

void printMat4(const mat4& M)
{
    cout << "**Mat4***" << endl;
    for(int i = 0; i < 4; i++){
        cout << M[0][i] <<", "<< M[1][i] <<", " << M[2][i] <<", " << M[3][i] << endl;
    }
    cout << endl;
}

void printMat3(const mat3& M)
{
    cout << "**Mat3***" << endl;
    for(int i = 0; i < 3; i++){
        cout << M[0][i] <<", "<< M[1][i] <<", " << M[2][i] << endl;
    }
    cout << endl;
}

/////////////////////////////////////////////////////////////////
/// \brief area
/// \param v
/// \return
///
class Convexhull2D {
public:
        friend int indexcompare(const void* a, const void* b);
        friend struct VectorindexCompare;

        Convexhull2D(glm::vec2* a, size_t n);
        Convexhull2D();
        Convexhull2D(const Convexhull2D&);
        ~Convexhull2D();

        void construct(glm::vec2* points, size_t n);

        size_t size() const;

        std::vector<glm::vec2> vertices;

private:
        class Vectorindex {
                public:
                        Vectorindex() : flag(true) {};
                        glm::vec2 pv;
                        size_t i;
                        bool flag;
                        Vectorindex* pi;
        };
        int iSize;
};


class Line2D {
public:
        friend bool intersect(glm::vec2& i, const Line2D& l1, const Line2D& l2);
        Line2D(const glm::vec2& a, const glm::vec2& b) {this->a = a; this->b = b; this->genBuffer();}
        Line2D ortho();
        Line2D parallel(glm::vec2& q);
        glm::vec2 maxPointToLine(std::vector<glm::vec2>& v);
        float distpoint(glm::vec2& p) const;
        glm::vec2& getA() {return a;}
        glm::vec2& getB() {return b;}

        // Debug
        void Draw();
        void genBuffer();
        ~Line2D();
private:
        glm::vec2 a;
        glm::vec2 b;
        unsigned int VAO, VBO;
};

float area(std::vector<vec2>& v); // area of convex polygon define by contour v


inline size_t Convexhull2D::size() const {
    return iSize;
};

float ccw(const vec2& p, const vec2& q, const vec2& r) {
    return (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
}

#define areaTriangle(p, q, r) (ccw(p, q, r) * 0.5f)

Convexhull2D::Convexhull2D(vec2* a, size_t n) : iSize(0) {
    construct(a, n);
}

Convexhull2D::Convexhull2D() : iSize(0) {

}

Convexhull2D::Convexhull2D(const Convexhull2D& ch) : iSize(0) {
    *this = ch;
}

Convexhull2D::~Convexhull2D() {

}

float area(vector<vec2>& v) {
    float fArea = 0.0f;
    for(int i = 1; i < v.size() - 1; i++) {
        fArea += static_cast<float>(fabs(areaTriangle(v[0], v[i], v[i + 1])));
    }

    return fArea;
}

int indexcompare(const void* a, const void* b) {
    Convexhull2D::Vectorindex va = *(Convexhull2D::Vectorindex*)a;
    Convexhull2D::Vectorindex vb = *(Convexhull2D::Vectorindex*)b;

    if(fabs(va.pv[0] - vb.pv[0]) < numeric_limits<float>::epsilon()) {
        if(va.pv[1] >= vb.pv[1])
            return 1;
        else
            return  0;
    } else if(va.pv[0] >= vb.pv[0]) {
        return  1;
    }

    return -1;
}

struct VectorindexCompare {
    int operator() (const Convexhull2D::Vectorindex& va, const Convexhull2D::Vectorindex& vb) const {
        if(fabs(va.pv[0] - vb.pv[0]) < numeric_limits<float>::epsilon()) {
            if(va.pv[1] >= vb.pv[1])
                return 1;
            else
                return  0;
        } else if(va.pv[0] >= vb.pv[0]) {
            return 1;
        }

        return 0;
    }
};

void Convexhull2D::construct(vec2* points, size_t n) {
    vector<Convexhull2D::Vectorindex> sortedIndex(n);
    vector<Vectorindex> chpoints;

    size_t nsize;

    int i, j, k1, k2;

    float f;

    for(i = 0; i < n; i++) {
        sortedIndex[i].i = i;
        sortedIndex[i].pv = points[i];
    }

    sort(sortedIndex.begin(), sortedIndex.end(), VectorindexCompare());

    chpoints.resize(n + 1);

    sortedIndex[0].flag = true;
    sortedIndex[1].flag = false;

    chpoints[0] = sortedIndex[0];
    chpoints[1] = sortedIndex[1];

    chpoints[0].pi = &sortedIndex[0];
    chpoints[1].pi = &sortedIndex[1];

    k1 = 2;

    for(i = 2; i < n; i++) {
        for(j = k1 - 1; j > 0; j--) {
            f = ccw(chpoints[j - 1].pv, points[sortedIndex[i].i], chpoints[j].pv);
            if(fabs(f) < numeric_limits<float>::epsilon() || f >= 0.0f) {
                chpoints[--k1].pi->flag = true;
            } else
                break;
        }
        sortedIndex[i].flag = false;
        chpoints[k1] = sortedIndex[i];
        chpoints[k1++].pi = &sortedIndex[i];
    }

    k2 = 0;

    for(i = n - 2; i >= 0; i--) {
        if(sortedIndex[i].flag) {
            for(j = k2 - 1; j >= 0; j--) {
                f = ccw(chpoints[k1 + j - 1].pv, points[sortedIndex[i].i], chpoints[k1 + j].pv);
                if(fabs(f) < numeric_limits<float>::epsilon() || f >= 0.0f)
                    k2--;
                else
                    break;
            }
            chpoints[k1 + k2++] = sortedIndex[i];
        }
    }

    nsize = k1 + k2 - 1;

    chpoints.resize(nsize);

    vertices.resize(iSize = nsize);

    for(i = 0; i < iSize; i++) {
        vertices[i] = chpoints[i].pv;
    }
}

/////////////////////////////////////////////////////////////////
bool intersect(vec2& i, const Line2D& l1, const Line2D& l2);

Line2D Line2D::ortho() {
    vec2 u = b - a;
    vec2 v(-u[1], u[0]);
    return Line2D(a, v + a);
}

Line2D Line2D::parallel(vec2& q) {
    vec2 u = b - a;
    vec2 p = q + u;
    return Line2D(q, p);
}

vec2 Line2D::maxPointToLine(std::vector<vec2>& v) {
    float fMax = 0.0f, f;
    int iMax = 0;

    for(int i = 0; i < v.size(); i++) {
        f = fabs(distpoint(v[i]));
        if(f > fMax)  {
            fMax = f;
            iMax = i;
        }
    }
    return v[iMax];
}

float Line2D::distpoint(vec2& p) const {
    return ((b[0] - a[0]) * (a[1] - p[1]) - (a[0] - p[0]) * (b[1] - a[1])) / sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));
}

bool intersect(vec2& i, const Line2D& l1, const Line2D& l2) {
    float denominator = (l1.b[1] * (l2.a[0] - l2.b[0]) + l1.a[1] * (l2.b[0] - l2.a[0]) + (l1.a[0] - l1.b[0]) * (l2.a[1] - l2.b[1]));

    //if(denominator < 0.001)
    //	return false;

    i = vec2((-l1.b[0] * l2.a[1] * l2.b[0] + l1.a[1] * l1.b[0] * (l2.b[0] - l2.a[0]) + l1.b[0] * l2.a[0] * l2.b[1] + l1.a[0] * (l1.b[1] * l2.a[0] - l2.b[1] * l2.a[0] - l1.b[1] * l2.b[0] + l2.a[1] * l2.b[0])) / denominator,
            ( l1.b[1] * (-l2.a[1] * l2.b[0] + l1.a[0] * (l2.a[1] - l2.b[1]) + l2.a[0] * l2.b[1]) + l1.a[1] * (l2.a[1] * l2.b[0] - l2.a[0] * l2.b[1] + l1.b[0] * (l2.b[1] - l2.a[1]))) / denominator);

    return true;
}

// Debug
Line2D::~Line2D()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
}
void Line2D::genBuffer()
{
    vec4 c(a,b);
    glGenBuffers(1, &VBO);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec4), value_ptr(c), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), 0);
}
void Line2D::Draw()
{
    vec4 c(a,b);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec4), value_ptr(c), GL_STATIC_DRAW);
    // Draw mesh
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINE_STRIP, 0, 2);
    glBindVertexArray(0);

}

/////////////////////////////////////////////////////////////////

vec4* setProjection(float l, float r, float b, float t, float n, float f, const mat4& view) {

    float lFar = (f / n) * l;
    float rFar = (f / n) * r;
    float bFar = (f / n) * b;
    float tFar = (f / n) * t;

    vec4 m_vtVertices[8];
    m_vtVertices[0] = vec4(l, b, -n, 1.0f);
    m_vtVertices[1] = vec4(r, b, -n, 1.0f);
    m_vtVertices[2] = vec4(r, t, -n, 1.0f);
    m_vtVertices[3] = vec4(l, t, -n, 1.0f);
    m_vtVertices[4] = vec4(lFar, bFar, -f, 1.0f);
    m_vtVertices[5] = vec4(rFar, bFar, -f, 1.0f);
    m_vtVertices[6] = vec4(rFar, tFar, -f, 1.0f);
    m_vtVertices[7] = vec4(lFar, tFar, -f, 1.0f);


    //mat4 Cinv = m_C; // view matrix
    mat4 Cinv = glm::inverse(view);

    vec4* m_vtVerticesTransformed = new vec4[8];

    for(int i = 0; i < 8; i++)
        m_vtVerticesTransformed[i] = Cinv * m_vtVertices[i];

    return m_vtVerticesTransformed;

}

vec4* setProjection(float fovy, float aspect, float n, float f, const mat4& view) {

    float l, r, b, t;
    t = n * tan(fovy); // radian
    b = -t;
    l = b * aspect;
    r = t * aspect;

    return setProjection(l, r, b, t, n, f, view);
}

/////////////////////////////////////////////////////////////////

#define ASSIGN_MAT(M, u0, u3, u6, u1, u4, u7, u2, u5, u8) { \
    M[0][0] = u0; M[0][1] = u3; M[0][2] = u6; \
    M[1][0] = u1; M[1][1] = u4; M[1][2] = u7; \
    M[2][0] = u2; M[2][1] = u5; M[2][2] = u8; \
    }


#define DET2(a, b, c, d) ((a) * (d) - (b) * (c))

#define DOT2(u, v) (u[0] * v[0] + u[1] * v[1])

// Forward declarations
void intersect(float i[2], const vec2& g0, const vec2& g1, const vec2& h0, const vec2& h1);
void map_Trapezoid_To_Square(mat3& TR, const vec2& t0, const vec2& t1, const vec2& t2, const vec2& t3);
void setTrapezoidResolution(mat3& N,
                            float percentage, float d1,
                            Line2D& bottom, Line2D& top,
                            Line2D& lineofsight, vec2 eyeTrans[8]);

// Functions
void computeTsmMatrix(mat4& N_T, const Camera& camera, const mat4& Ml,
                      /*float tsmDistance, (not sure how to use this)*/
                      float percentage, int shadowmapHeight) {
    //Camera visEye = eye;
    //visEye.setFarPlaneDistance(tsmDistance); // what is tsmDistance?

    vec2 eyeTrans[8];
    vec2 eyeVisTrans[8];

    //mat4 Ml = light.m_P * light.m_C; // m_P frustum of light, m_C view of light

    //vec4* p1 = eye.m_vtVerticesTransformed;	// get transformed vertices for near and far region
    //vec4* p2 = visEye.m_vtVerticesTransformed;

    // Transfer frustum vetices to world coordinates
    vec4* p1 = setProjection(radians(camera.Zoom), camera.Aspect, camera.Near, 50, camera.GetViewMatrix()); // (fovy, aspect, near, far), also need view matrix
    vec4* p2 = setProjection(radians(camera.Zoom), camera.Aspect, camera.Near, 100, camera.GetViewMatrix());

    for(int i = 0; i < 8; i++) {
        vec4 v = Ml * p1[i]; // transform to lightspace
        eyeTrans[i][0] = v[0] / v[3]; // normalize
        eyeTrans[i][1] = v[1] / v[3];

        v = Ml * p2[i];
        eyeVisTrans[i][0] = v[0] / v[3];
        eyeVisTrans[i][1] = v[1] / v[3];
    }

    vec4 cu = 0.5f * (p1[2] - p1[0]) + p1[0]; // center point of near plane
    vec4 cv = 0.5f * (p1[6] - p1[4]) + p1[4]; // center point of  far plane

    // Convert to lightspace
    cu = Ml * cu;
    cv = Ml * cv;

    delete [] p1;
    delete [] p2;

    // Line of sight
    Line2D losEye(vec2(cu[0] / cu[3], cu[1] / cu[3]), vec2(cv[0] / cv[3], cv[1] / cv[3]));

    // Create convex hull for near&far planes
    Convexhull2D ch1(eyeTrans, 8);
    Convexhull2D ch2(eyeVisTrans, 8);

    // Get base lines
    Line2D orthoLosEye = losEye.ortho();
    vec2 maxPointA = orthoLosEye.maxPointToLine(ch1.vertices);

    Line2D bottom = orthoLosEye.parallel(maxPointA);
    vec2 maxPointB = bottom.maxPointToLine(ch1.vertices);

    Line2D top = bottom.parallel(maxPointB);
    vec2 maxPointC = top.maxPointToLine(ch2.vertices);

    float d1 = fabs(top.distpoint(maxPointC));

    // determine 80% rule Matrix
    float maxArea = 0.0f;
    float f;
    float fInc = 100.0f / shadowmapHeight;
    float g = percentage;

    mat3 N;
    for(g = percentage; g >= 0.0f; g -= fInc) {
        setTrapezoidResolution(N, g, d1, bottom, top, losEye, eyeTrans);

        vector<vec2> ch(ch2.size());

        for(int i = 0; i < ch2.size(); i++) {
            vec3 v = N * vec3(ch2.vertices[i][0], ch2.vertices[i][1], 1.0f);
            ch[i] = vec2(v[0] / v[2], v[1] / v[2]);
        }

        f = area(ch);

        if(f < maxArea) {
            g += fInc;
            break;
        } else {
            maxArea = f;
        }
    }
    //cout<< g << endl;
    //mat3 N;
    //setTrapezoidResolution(N, g, d1, bottom, top, losEye, eyeTrans);

    // Complete
    N_T[0] = vec4(N[0].x,N[0].y,0.0,N[0].z);
    N_T[1] = vec4(N[1].x,N[1].y,0.0,N[1].z);
    N_T[2] = vec4(0.0,   0.0,   1.0,   0.0);
    N_T[3] = vec4(N[2].x,N[2].y,0.0,N[2].z);

}
void setTrapezoidResolution(mat3& N,
                            float percentage, float d1,
                            Line2D& bottom, Line2D& top,
                            Line2D& lineofsight, vec2 eyeTrans[8]) {

    vec2 s = bottom.getA();
    float l = fabs(top.distpoint(s)); // bottom-top distance
    float d2 = -((2.0f * percentage) / 100.0f - 1.0f);
    float n = (l * d1 + l * d1 * d2) / (l - 2.0f * d1 - l * d2);

    vec2 u, v, w, i;
    intersect(u, top, lineofsight);
    intersect(v, bottom, lineofsight);
    w = v - u;

    float t = -n / glm::length(w);
    i = u - w * t; // orig: +

    float min1, min2, f;
    int i0 = 0, i1 = 0;
    min1 = min2 = numeric_limits<float>::max();
    w = i - u;

    w = glm::normalize(w);

    for(int j = 4; j < 8; j++) {
        vec2 x = eyeTrans[j] - i;

        x = glm::normalize(x);

        f = static_cast<float>(acos(glm::dot(x, w)));
        if(ccw(u, v, eyeTrans[j]) < 0.0f) {
            if(f < min1) {
                min1 = f;
                i0 = j;
            }
        } else {
            if(f < min2) {
                min2 = f;
                i1 = j;
            }
        }
    }

    vec2 left = eyeTrans[i0], right = eyeTrans[i1];

    Line2D h(i, right);
    Line2D g(i, left);

    vec2 t0, t1, t2, t3;
    if(!intersect(t0, bottom, h)) cerr << "couldn't calc t0" << endl;
    if(!intersect(t1, bottom, g)) cerr << "couldn't calc t1" << endl;
    if(!intersect(t2, top, g)) cerr << "couldn't calc t2" << endl;
    if(!intersect(t3, top, h)) cerr << "couldn't calc t3" << endl;

    //map_Trapezoid_To_Square(N, t0, t1, t2, t3);
    //N = transpose(N);

    {
        vec3 t_0, t_1, t_2, t_3;
        t_0 = vec3(t0, 1.0);t_1 = vec3(t1, 1.0);
        t_2 = vec3(t2, 1.0);t_3 = vec3(t3, 1.0);
        // 1.Translate
        vec3 u = 0.5f * (t_2 + t_3);
        mat3 T_1 = mat3(1,0,0,
                        0,1,0,
                        -u.x,-u.y,1);

        // 2.Rotate
        u = normalize(t_2 - t_3);
        mat3 R = mat3(-u.x,-u.y,0,
                      -u.y,u.x,0,
                      0,0,1);

        mat3 TR = R * T_1;

        // 3.Translate
        float i[2];
        intersect(i, t0, t3, t1, t2);

        u = TR * vec3(i[0],i[1],1.0);
        mat3 T_2 = mat3(1,0,0,
                        0,1,0,
                        -u.x,-u.y,1);
        // 4.Shear
        TR = T_2 * TR;
        u = (TR *(t_2 + t_3)) / 2.0f;
        mat3 H = mat3(1, 0, 0,
                      -u.x/u.y,1,0,
                      0,0,1);
        // 5.Scale
        TR = H * TR;
        u = TR * t_2;

        mat3 S_1 = mat3(1/u.x,0,0,
                        0,1/u.y,0,
                        0,0,1);
        // 6.Trapezoid
        mat3 Nt = mat3(1,0,0,
                       0,1,1,
                       0,1,0);

        // 7.Translating
        TR = Nt * S_1 * TR;
        u = TR * t_0;
        vec3 v = TR * t_2;

        mat3 T_3 = mat3(1, 0,  0 ,
                        0, 1,  0,
                        0, -(u.y/u.z+v.y/v.z)/2,  1);
        // 8.Scale to unit cube
        TR = T_3 * TR;
        u = TR * t_0;

        mat3 S_2 =mat3(1,     0,      0,
                       0, -u.z/u.y,   0,
                       0,     0 ,    1);


        N = S_2 * TR * mat3(1.0);
    }
}

void intersect(float i[2], const vec2& g0, const vec2& g1, const vec2& h0, const vec2& h1) {
    float a, b;

    i[0] = i[1] = 1.0f / DET2(g0[0] - g1[0], g0[1] - g1[1], h0[0] - h1[0], h0[1] - h1[1]);

    a = DET2(g0[0], g0[1], g1[0], g1[1]);
    b = DET2(h0[0], h0[1], h1[0], h1[1]);

    i[0] *=	DET2(a, g0[0] - g1[0], b, h0[0] - h1[0]);
    i[1] *=	DET2(a, g0[1] - g1[1], b, h0[1] - h1[1]);
}

void map_Trapezoid_To_Square(mat3& TR, 
                             const vec2& t0,
                             const vec2& t1,
                             const vec2& t2,
                             const vec2& t3) {
    float i[2], a, b, c, d;

    //M1 = R * T1
    a = 0.5f * (t2[0] - t3[0]);
    b = 0.5f * (t2[1] - t3[1]);

    ASSIGN_MAT(TR, a  ,  b  , a * a + b * b,
               b  , -a  , a * b - b * a,
               0.0f, 0.0f, 1.0f);

    //M2 = T2 * M1 = T2 * R * T1
    intersect(i, t0, t3, t1, t2);

    TR[0][2] = -DOT2(TR[0], i);
    TR[1][2] = -DOT2(TR[1], i);

    //M1 = H * M2 = H * T2 * R * T1
    a = DOT2(TR[0], t2) + TR[0][2];
    b = DOT2(TR[1], t2) + TR[1][2];
    c = DOT2(TR[0], t3) + TR[0][2];
    d = DOT2(TR[1], t3) + TR[1][2];

    a = -(a + c) / (b + d);

    TR[0][0] += TR[1][0] * a;
    TR[0][1] += TR[1][1] * a;
    TR[0][2] += TR[1][2] * a;

    //M2 = S1 * M1 = S1 * H * T2 * R * T1
    a = 1.0f / (DOT2(TR[0], t2) + TR[0][2]);
    b = 1.0f / (DOT2(TR[1], t2) + TR[1][2]);

    TR[0][0] *= a; TR[0][1] *= a; TR[0][2] *= a;
    TR[1][0] *= b; TR[1][1] *= b; TR[1][2] *= b;

    //M1 = N * M2 = N * S1 * H * T2 * R * T1
    TR[2][0] = TR[1][0]; TR[2][1] = TR[1][1]; TR[2][2] = TR[1][2];
    TR[1][2] += 1.0f;

    //M2 = T3 * M1 = T3 * N * S1 * H * T2 * R * T1
    a = DOT2(TR[1], t0) + TR[1][2];
    b = DOT2(TR[2], t0) + TR[2][2];
    c = DOT2(TR[1], t2) + TR[1][2];
    d = DOT2(TR[2], t2) + TR[2][2];

    a = -0.5f * (a / b + c / d);

    TR[1][0] += TR[2][0] * a;
    TR[1][1] += TR[2][1] * a;
    TR[1][2] += TR[2][2] * a;

    //M1 = S2 * M2 = S2 * T3 * N * S1 * H * T2 * R * T1
    a = DOT2(TR[1], t0) + TR[1][2];
    b = DOT2(TR[2], t0) + TR[2][2];

    c = -b / a;

    TR[1][0] *= c; TR[1][1] *= c; TR[1][2] *= c;
}

// Helper
void debugTSM(const Camera& camera, const glm::mat4& Ml, float percentage){

    // Helper class
    /*class _functor
    {
    public:
        float ccw(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) {
            return (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
        }

        void intersect(float i[2], const glm::vec2& g0, const glm::vec2& g1, const glm::vec2& h0, const glm::vec2& h1) {
            float a, b;

            i[0] = i[1] = 1.0f / DET2(g0[0] - g1[0], g0[1] - g1[1], h0[0] - h1[0], h0[1] - h1[1]);

            a = DET2(g0[0], g0[1], g1[0], g1[1]);
            b = DET2(h0[0], h0[1], h1[0], h1[1]);

            i[0] *=	DET2(a, g0[0] - g1[0], b, h0[0] - h1[0]);
            i[1] *=	DET2(a, g0[1] - g1[1], b, h0[1] - h1[1]);
        }
        glm::vec4* setProjection(float l, float r, float b, float t, float n, float f, const glm::mat4& view) {

            float lFar = (f / n) * l;
            float rFar = (f / n) * r;
            float bFar = (f / n) * b;
            float tFar = (f / n) * t;

            glm::vec4 m_vtVertices[8];
            m_vtVertices[0] = glm::vec4(l, b, -n, 1.0f);
            m_vtVertices[1] = glm::vec4(r, b, -n, 1.0f);
            m_vtVertices[2] = glm::vec4(r, t, -n, 1.0f);
            m_vtVertices[3] = glm::vec4(l, t, -n, 1.0f);
            m_vtVertices[4] = glm::vec4(lFar, bFar, -f, 1.0f);
            m_vtVertices[5] = glm::vec4(rFar, bFar, -f, 1.0f);
            m_vtVertices[6] = glm::vec4(rFar, tFar, -f, 1.0f);
            m_vtVertices[7] = glm::vec4(lFar, tFar, -f, 1.0f);


            //mat4 Cinv = m_C; // view matrix
            glm::mat4 Cinv = glm::inverse(view);

            glm::vec4* m_vtVerticesTransformed = new glm::vec4[8];

            for(int i = 0; i < 8; i++)
                m_vtVerticesTransformed[i] = Cinv * m_vtVertices[i];

            return m_vtVerticesTransformed;

        }

        glm::vec4* setProjection(float fovy, float aspect, float n, float f, const glm::mat4& view) {

            float l, r, b, t;
            t = n * tan(fovy); // radian
            b = -t;
            l = b * aspect;
            r = t * aspect;

            return setProjection(l, r, b, t, n, f, view);
        }

    private:
        float DET2(float a, float b, float c, float d)
        {
            return (a) * (d) - (b) * (c);
        }
    } _f;*/


    glm::vec2 eyeTrans[8];
    glm::vec2 eyeVisTrans[8];

    // Transfer frustum vetices to world coordinates
    glm::vec4* p1 = setProjection(glm::radians(camera.Zoom), camera.Aspect, camera.Near, 50, camera.GetViewMatrix()); // (fovy, aspect, near, far), also need view matrix
    glm::vec4* p2 = setProjection(glm::radians(camera.Zoom), camera.Aspect, camera.Near, 100, camera.GetViewMatrix());

    for(int i = 0; i < 8; i++) {
        glm::vec4 v = Ml * p1[i]; // transform to lightspace
        //std::cout << "p" << i << " " << v.x << ", "<< v.y << ", "<< v.z << ", "<< v.w << std::endl;
        eyeTrans[i][0] = v[0] / v[3]; // normalize
        eyeTrans[i][1] = v[1] / v[3];

        v = Ml * p2[i];
        eyeVisTrans[i][0] = v[0] / v[3];
        eyeVisTrans[i][1] = v[1] / v[3];
    }

    glm::vec4 u = 0.5f * (p1[2] - p1[0]) + p1[0]; // center point of near plane
    glm::vec4 v = 0.5f * (p1[6] - p1[4]) + p1[4]; // center point of  far plane

    // Convert to lightspace
    u = Ml * u;
    v = Ml * v;

    // Line of sight
    Line2D losEye(glm::vec2(u[0] / u[3], u[1] / u[3]), glm::vec2(v[0] / v[3], v[1] / v[3]));

    // Create convex hull for near&far planes
    Convexhull2D ch1(eyeTrans, 8);
    Convexhull2D ch2(eyeVisTrans, 8);

    // Get base lines
    Line2D orthoLosEye = losEye.ortho();
    glm::vec2 maxPointA = orthoLosEye.maxPointToLine(ch1.vertices);

    Line2D bottom = orthoLosEye.parallel(maxPointA);
    glm::vec2 maxPointB = bottom.maxPointToLine(ch1.vertices);

    Line2D top = bottom.parallel(maxPointB);
    glm::vec2 maxPointC = top.maxPointToLine(ch2.vertices);

    float d1 = fabs(top.distpoint(maxPointC));

    glm::vec2 s = bottom.getA();
    float l = fabs(top.distpoint(s)); // bottom-top distance
    float d2 = -((2.0f * percentage) / 100.0f - 1.0f);
    float n = (l * d1 + l * d1 * d2) / (l - 2.0f * d1 - l * d2);

    glm::vec2 ua, va, w, i;
    intersect(ua, top, losEye);
    intersect(va, bottom, losEye);
    w = va - ua;

    float t = -n / glm::length(w);
    i = ua - w * t;

    float min1, min2, f;
    int i0 = 0, i1 = 0;
    min1 = min2 = numeric_limits<float>::max();
    w = i - ua;

    w = glm::normalize(w);

    for(int j = 4; j < 8; j++) {
        glm::vec2 x = eyeTrans[j] - i;

        x = glm::normalize(x);

        f = static_cast<float>(acos(glm::dot(x, w)));

        if(ccw(ua, va, eyeTrans[j]) < 0.0f) {
            if(f < min1) {
                min1 = f;
                i0 = j;
            }
        } else {
            if(f < min2) {
                min2 = f;
                i1 = j;
            }
        }
    }

    glm::vec2 left = eyeTrans[i0], right = eyeTrans[i1];

    Line2D h(i, right);
    Line2D g(i, left);

    //h.Draw();
    //g.Draw();

    glm::vec2 t0, t1, t2, t3;
    if(!intersect(t0, bottom, h)) cerr << "couldn't calc t0" << endl;
    if(!intersect(t1, bottom, g)) cerr << "couldn't calc t1" << endl;
    if(!intersect(t2, top, g)) cerr << "couldn't calc t2" << endl;
    if(!intersect(t3, top, h)) cerr << "couldn't calc t3" << endl;

    glm::mat4 N_T;
    {
        glm::mat3 N;
        glm::vec3 t_0, t_1, t_2, t_3;
        t_0 = glm::vec3(t0, 1.0);t_1 = glm::vec3(t1, 1.0);
        t_2 = glm::vec3(t2, 1.0);t_3 = glm::vec3(t3, 1.0);
        // 1.Translate
        glm::vec3 u = 0.5f * (t_2 + t_3);
        glm::mat3 T_1 = glm::transpose(glm::mat3(1,0,-u.x,
                                                 0,1,-u.y,
                                                 0,0,1));

        // 2.Rotate
        u = glm::normalize(t_2 - t_3);
        glm::mat3 R = glm::mat3(-u.x,-u.y,0,
                                -u.y,u.x,0,
                                0,0,1);

        // 3.Translate
        float i[2];
        intersect(i, t0, t3, t1, t2);

        u = R * T_1 * glm::vec3(i[0],i[1],1.0);
        glm::mat3 T_2 = glm::transpose(glm::mat3(1,0,-u.x,
                                                 0,1,-u.y,
                                                 0,0,1));
        // 4.Shear
        u = (T_2 * R * T_1 *(t_2 + t_3)) / 2.0f;
        glm::mat3 H = glm::mat3(1, 0, 0,
                                -u.x/u.y,1,0,
                                0,0,1);
        // 5.Scale
        u = H * T_2 * R * T_1 * t_2;

        glm::mat3 S_1 = glm::mat3(1/u.x,0,0,
                                  0,1/u.y,0,
                                  0,0,1);
        // 6.Trapezoid
        glm::mat3 Nt = glm::mat3(1,0,0,
                                 0,1,1,
                                 0,1,0);

        // 7.Translating
        u = Nt * S_1 * H * T_2 * R * T_1 * t_0;
        glm::vec3 v = Nt * S_1 * H * T_2 * R * T_1 * t_2;

        glm::mat3 T_3 = glm::mat3(1, 0,  0 ,
                                  0, 1,  0,
                                  0, -(u.y/u.z+v.y/v.z)/2,  1);
        // 8.Scale to unit cube
        u = T_3 * Nt * S_1 * H * T_2 * R * T_1 * t_0;

        glm::mat3 S_2 =glm::mat3(1,     0,      0,
                                 0, -u.z/u.y,   0,
                                 0,     0 ,    1);


        N = S_2 * T_3 * Nt * S_1 * H * T_2 * R * T_1 * glm::mat3(1.0);

        N_T[0] = glm::vec4(N[0].x,N[0].y,0.0,N[0].z);
        N_T[1] = glm::vec4(N[1].x,N[1].y,0.0,N[1].z);
        N_T[2] = glm::vec4(0.0,   0.0,   1.0,   0.0);
        N_T[3] = glm::vec4(N[2].x,N[2].y,0.0,N[2].z);
    }

    {
        glm::vec4 t[4];
        t[0] = N_T * glm::vec4(t0,0,1);
        t[1] = N_T * glm::vec4(t1,0,1);
        t[2] = N_T * glm::vec4(t2,0,1);
        t[3] = N_T * glm::vec4(t3,0,1);
        for(int i = 0; i < 4; i++){
            t[i] /= t[i].w;
            std::cout << "t" << i << " " << t[i].x << ", "<< t[i].y << ", "<< t[i].z << ", "<< t[i].w << std::endl;
        }

        Line2D n1(t[0], t[1]);
        Line2D n2(t[1], t[2]);
        Line2D n3(t[2], t[3]);
        Line2D n4(t[3], t[0]);

        n1.Draw();
        n2.Draw();
        n3.Draw();
        n4.Draw();
    }

    for(int i = 0; i < 8; i++) {
        glm::vec4 v = N_T * Ml * p1[i]; // transform to lightspace
        //std::cout << "p" << i << " " << v.x << ", "<< v.y << ", "<< v.z << ", "<< v.w << std::endl;
        eyeTrans[i][0] = v[0] / v[3]; // normalize
        eyeTrans[i][1] = v[1] / v[3];
    }

    {
        Line2D n1(eyeTrans[0],eyeTrans[1]);
        n1.Draw();
        Line2D n2(eyeTrans[1],eyeTrans[2]);
        n2.Draw();
        Line2D n3(eyeTrans[2],eyeTrans[3]);
        n3.Draw();
        Line2D n4(eyeTrans[3],eyeTrans[0]);
        n4.Draw();

        Line2D f1(eyeTrans[4],eyeTrans[5]);
        f1.Draw();
        Line2D f2(eyeTrans[5],eyeTrans[6]);
        f2.Draw();
        Line2D f3(eyeTrans[6],eyeTrans[7]);
        f3.Draw();
        Line2D f4(eyeTrans[7],eyeTrans[4]);
        f4.Draw();

        Line2D l1(eyeTrans[0],eyeTrans[4]);
        l1.Draw();
        Line2D l2(eyeTrans[1],eyeTrans[5]);
        l2.Draw();
        Line2D l3(eyeTrans[2],eyeTrans[6]);
        l3.Draw();
        Line2D l4(eyeTrans[3],eyeTrans[7]);
        l4.Draw();
    }
    // Corss Lines
    {
        Line2D hori(glm::vec2(-1,0), glm::vec2(1,0)) ;
        Line2D vert(glm::vec2(0,-1), glm::vec2(0,1)) ;

        hori.Draw();
        vert.Draw();
    }

    delete [] p1;
    delete [] p2;
}
