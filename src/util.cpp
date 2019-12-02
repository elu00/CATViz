#include <iostream>
#include <cmath>

#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
const double PI_2 = 1.57079632679489661923;
const double PI_4 = 0.785398163397448309616;

using namespace glm;
using namespace std;

double angle(vec2 v1, vec2 v2) 
{
    return acos(dot(normalize(v1), normalize(v2)));
}
vec2 prod(vec2 a, vec2 b)
{
    return vec2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
}
vec2 div(vec2 a, vec2 b)
{
    return vec2(((a.x*b.x+a.y*b.y)/(b.x*b.x+b.y*b.y)),((a.y*b.x-a.x*b.y)/(b.x*b.x+b.y*b.y)));
}
double sc(double c) 
{
    return sqrt(c * c - 1.);
}

struct S2Triangle
{
    vec3 i, j, k;
    double area;
};
struct H2Triangle
{
    vec3 i, j, k;
    double area;
};

struct E2Triangle
{
    vec2 i, j, k;
};
struct CAT
{
    vec2 i, j, k;
    double a_ij, a_jk, a_ki;
};
struct MobiusInfo
{
    vec2 a, b, c, d;
};

vec2 ForwardMobius(vec2 z, MobiusInfo m)
{
    return div((prod(m.a, z) + m.b), (prod(m.c, z) + m.d));
}

vec2 InvMobius(vec2 z, MobiusInfo m)
{
    return div((prod(m.d, z) - m.b), (-prod(m.c, z) + m.a));
}

vec2 ComplexDet(vec2 a, vec2 b, vec2 c, vec2 d, vec2 e, vec2 f, vec2 g, vec2 h, vec2 i)
{
    return prod(a, (prod(e, i) - prod(f, h))) - prod(b, (prod(d, i) - prod(f, g))) + prod(c, (prod(d, h) - prod(e, g)));
}

vec2 ForwardSphericalProj(vec3 wh)
{
    return vec2(wh.x / (1. - wh.z), wh.y / (1. - wh.z));
}
vec3 InverseSphericalProj(vec2 z)
{
    //GLM PORT
    return vec3((2.f * z)/ (1.f + dot(z, z)), (-1.f + dot(z, z)) / (1.f + dot(z, z)));
}

double SphericalArea(vec3 A, vec3 B, vec3 C)
{
    return 2.*atan((dot(A,cross(B, C)))/(1. + dot(A, B) + dot(B, C) + dot(A,C)));
}
vec3 SphericalBary(vec3 p, S2Triangle T)
{
    double u = abs(SphericalArea(T.j, T.k, p) / T.area);
    double v = abs(SphericalArea(T.i, T.k, p) / T.area);
    double w = abs(SphericalArea(T.i, T.j, p) / T.area);
    if (u + v + w > 1.01)
    {
        return vec3(1., 1., 1.);
    }
    return vec3(u, v, w);
}
vec2 ForwardHyperbolicProj(vec3 wh)
{
    return vec2(wh.x / (1. + wh.z), wh.y / (1. + wh.z));
}
vec3 InverseHyperbolicProj(vec2 z)
{
    return vec3(2.f * z / (1.f - dot(z, z)),(1.f + dot(z, z)) / (1.f - dot(z, z)));
}
double HDist(vec3 i, vec3 j)
{
    return ((i.x * j.x + i.y * j.y) - i.z * j.z);
}

double HyperbolicArea(vec3 A, vec3 B, vec3 C)
{
    double cosha = -HDist(B, C);
    double coshb = -HDist(A, C);
    double coshc = -HDist(B, A);
    double alpha = acos((-cosha + coshb * coshc) / (sc(coshb) * sc(coshc)));
    double beta = acos((-coshb + cosha * coshc) / (sc(cosha) * sc(coshc)));
    double gamma = acos((-coshc + coshb * cosha) / (sc(coshb) * sc(cosha)));
    return PI - alpha - beta - gamma;
}
vec3 HyperbolicBary(vec3 p, H2Triangle T)
{
    double u = HyperbolicArea(T.j, T.k, p) / T.area;
    double v = HyperbolicArea(T.i, T.k, p) / T.area;
    double w = HyperbolicArea(T.i, T.j, p) / T.area;
    if (isnan(u) || isnan(v) || u + v + w > 1.01)
    {
        return vec3(1.,1.,1.);
    }
    return vec3(u, v, w);
}

vec3 BaryCalc(vec2 p, E2Triangle T)
{
    double u, v, w;
    vec2 v0 = T.j - T.i, v1 = T.k - T.i, v2 = p - T.i;
    double d00 = dot(v0, v0);
    double d01 = dot(v0, v1);
    double d11 = dot(v1, v1);
    double d20 = dot(v2, v0);
    double d21 = dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.f - v - w;
    if (u < 0. || v < 0. || w < 0.)
    {
        return vec3(1., 1., 1.);
    }
    return vec3(u, v, w);
}

MobiusInfo SolveTrans(vec2 z1, vec2 z2, vec2 z3, vec2 w1, vec2 w2, vec2 w3)
{
    vec2 z1w1 = prod(z1, w1);
    vec2 z2w2 = prod(z2, w2);
    vec2 z3w3 = prod(z3, w3);
    vec2 a = ComplexDet(z1w1, w1, vec2(1., 0), z2w2, w2, vec2(1, 0), z3w3, w3, vec2(1, 0));
    vec2 b = ComplexDet(z1w1, z1, w1, z2w2, z2, w2, z3w3, z3, w3);
    vec2 c = ComplexDet(z1, w1, vec2(1., 0), z2, w2, vec2(1., 0), z3, w3, vec2(1, 0));
    vec2 d = ComplexDet(z1w1, z1, vec2(1., 0), z2w2, z2, vec2(1., 0), z3w3, z3, vec2(1, 0));
    return MobiusInfo{a, b, c, d};
}

vec2 baryToPlane (vec2 i, vec2 j, vec2 k, double a_ij, double a_jk, double a_ki, float i_b, float j_b, float k_b)
{
    double alpha = angle(j - i, k - i) + a_ij + a_ki;
    double beta = angle(i - j, k - j) + a_ij + a_jk;
    double gamma = angle(i - k, j - k) + a_ki + a_jk;
    if (a_ij + a_jk + a_ki == 0.)
    {
        // Euclidean case
        vec2 ti, tj, tk;
        ti = vec2(0., 0.);
        double t_i = alpha;
        double t_j = beta;
        tj = vec2(cos(t_i) + sin(t_i) / tan(t_j), 0.);
        tk = vec2(cos(t_i), sin(t_i));
        MobiusInfo M = SolveTrans(ti, tj, tk, i, j, k);
        return ForwardMobius(i_b * ti + j_b * tj + k_b * tk, M);
    }
    else if (a_ij + a_jk + a_ki > 0.)
    {
        // spherical case
        vec3 ti, tj, tk;
        ti = vec3(0., 0., -1.);
   
        // law of cosines stuff to find side lenghts
        // a = length of jk, b = length of ik, c = length of ik
        double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
        double tj1 = sqrt(1. - c * c);
        tj = vec3(tj1, 0., -c);
        double tk1 = (a - b * c) / tj1;
        double tk2 = sqrt(1. - tk1 * tk1 - b * b);
        tk = vec3(tk1, tk2, -b);
        MobiusInfo M = SolveTrans(ForwardSphericalProj(ti), ForwardSphericalProj(tj), ForwardSphericalProj(tk), i, j, k);
        vec3 coord = i_b * ti + j_b * tj + k_b * tk;
        coord = coord/l2Norm(coord);
        return ForwardMobius(ForwardSphericalProj(coord), M);
    }
    else
    {
        // hyperbolic case
        vec3 ti, tj, tk;
        // law of cosines stuff to find side lenghts
        // a = length of jk, b = length of ik, c = length of ij
        double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
        double tj1 = sqrt(c * c - 1.);
        
        double tk1 = (-a + b * c) / tj1;
        double tk2 = sqrt(-1. - tk1 * tk1 + b * b);
        
        ti = vec3(0., 0., 1.);
        tj = vec3(tj1, 0., c);
        tk = vec3(tk1, tk2, b);
        MobiusInfo M = SolveTrans(ForwardHyperbolicProj(ti), ForwardHyperbolicProj(tj), ForwardHyperbolicProj(tk), i, j, k);
        vec3 coord = i_b * ti + j_b * tj + k_b * tk;
        coord = -coord/(coord.x * coord.x + coord.y * coord.y - coord.z * coord.z);
        return ForwardMobius(ForwardHyperbolicProj(coord), M);
    }

}
double l2DistSquared(double ij, double jk, double ki, double a_ij, double a_jk, double a_ki, double i1, double j1, double k1, double i2, double k2, double j2)
{
    vec2 i = vec2(0, 0);
    vec2 j = vec2(ij, 0);
    double angle = acos((-jk*jk + ki*ki + ij*ij)/(2*ki*ij));
    vec2 k = vec2(ki * cos(angle), ki * sin(angle));
    vec2 p1 = baryToPlane (i, j, k, a_ij, a_jk, a_ki, i1, j1, k1);
    vec2 p2 = baryToPlane (i, j, k, a_ij, a_jk, a_ki, i2, j2, k2);
    return ((p1-p2).x)*((p1-p2).x) + ((p1-p2).y)*((p1-p2).y);
}


