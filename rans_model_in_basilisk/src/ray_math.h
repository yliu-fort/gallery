#ifndef RAY_MATH_H
#define RAY_MATH_H

typedef struct _Vec3{
    double x,y,z;
} Vec3;

Vec3 v3_sub(const Vec3 self, const Vec3 v){
    return (Vec3){self.x - v.x,
                self.y - v.y,
                self.z - v.z};
}

double v3_dot(const Vec3 self, const Vec3 v){
    return self.x * v.x + self.y * v.y + self.z * v.z;
}

Vec3 v3_cross(const Vec3 self, const Vec3 v){
    return (Vec3){self.y * v.z - self.z * v.y,
                self.z * v.x - self.x * v.z,
                self.x * v.y - self.y * v.x};
}

double v3_length(const Vec3 self){
    return sqrt(self.x * self.x +
                        self.y * self.y +
                        self.z * self.z);
}

void v3_normalize(Vec3 self){
    double l = v3_length(self);
    l = l > 0 ? l : 1;
    self.x /= l;
    self.y /= l;
    self.z /= l;
}

typedef struct _Ray{
    Vec3 orig;
    Vec3 direction;
} Ray;

typedef struct _Hit{
    double distance;
    double u; // hit point in uv coordinate
    double v; // hit point in uv coordinate
    bool hit;
} Hit;

Hit ray_triangle_intersect(const Ray r, Vec3 v0, Vec3 v1, const Vec3 v2){
    const Vec3 nvec = v3_cross(v3_sub(v1, v0), v3_sub(v2, v0));
    const bool flipped = v3_dot(r.direction, nvec) >= 0;
    if (flipped){
        Vec3 _v0 = v0;
        v0 = v1;
        v1 = _v0;
    }

    const Vec3 v0v1 = v3_sub(v1, v0);
    const Vec3 v0v2 = v3_sub(v2, v0);
    const Vec3 pvec = v3_cross(r.direction, v0v2);

    const double det = v3_dot(v0v1, pvec);

    if (det < 0.000001)
        return (Hit){-nodata,nodata,nodata,0};

    const double invDet = 1.0 / det;
    const Vec3 tvec = v3_sub(r.orig, v0);
    double u = v3_dot(tvec, pvec) * invDet;

    if (u < 0 || u > 1)
        return (Hit){-nodata,nodata,nodata,0};

    const Vec3 qvec = v3_cross(tvec, v0v1);
    double v = v3_dot(r.direction, qvec) * invDet;

    if (v < 0 || (u + v) > 1)
        return (Hit){-nodata,nodata,nodata,0};

    if (flipped)
    {
        double _u = u;
        u = v;
        v = _u;
    }

    return (Hit){v3_dot(v0v2, qvec) * invDet, u, v, 1};
}

#endif