#include <iostream>
#include <cmath>

#include "glm/glm.hpp"

// Port of the written shader using GLM.
// Notable differences: doubles used instead of floats, no swizzling available so some hackier things are done
const double PI = 3.1415926535897932384626433832795;
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
float sc(float c) 
{
    return sqrt(c * c - 1.);
}

struct S2Triangle
{
    vec3 i, j, k;
    float area;
};
struct H2Triangle
{
    vec3 i, j, k;
    float area;
};

struct E2Triangle
{
    vec2 i, j, k;
};
struct CAT
{
    vec2 i, j, k;
    float a_ij, a_jk, a_ki;
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

float SphericalArea(vec3 A, vec3 B, vec3 C)
{
    return 2.*atan((dot(A,cross(B, C)))/(1. + dot(A, B) + dot(B, C) + dot(A,C)));
    /*
    float cosa = dot(B, C);
    float cosb = dot(A, C);
    float cosc = dot(B, A);
    float alpha = acos((cosa - cosb * cosc) / (sc(cosb) * sc(cosc)));
    float beta = acos((cosb - cosa * cosc) / (sc(cosa) * sc(cosc)));
    float gamma = acos((cosc - cosb * cosa) / (sc(cosb) * sc(cosa)));
    return alpha + beta + gamma - PI;
	*/
}
vec3 SphericalBary(vec3 p, S2Triangle T)
{
    float u = abs(SphericalArea(T.j, T.k, p) / T.area);
    float v = abs(SphericalArea(T.i, T.k, p) / T.area);
    float w = abs(SphericalArea(T.i, T.j, p) / T.area);
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
float HDist(vec3 i, vec3 j)
{
    return ((i.x * j.x + i.y * j.y) - i.z * j.z);
}

float HyperbolicArea(vec3 A, vec3 B, vec3 C)
{
    float cosha = -HDist(B, C);
    float coshb = -HDist(A, C);
    float coshc = -HDist(B, A);
    float alpha = acos((-cosha + coshb * coshc) / (sc(coshb) * sc(coshc)));
    float beta = acos((-coshb + cosha * coshc) / (sc(cosha) * sc(coshc)));
    float gamma = acos((-coshc + coshb * cosha) / (sc(coshb) * sc(cosha)));
    return PI - alpha - beta - gamma;
}
vec3 HyperbolicBary(vec3 p, H2Triangle T)
{
    float u = HyperbolicArea(T.j, T.k, p) / T.area;
    float v = HyperbolicArea(T.i, T.k, p) / T.area;
    float w = HyperbolicArea(T.i, T.j, p) / T.area;
    if (isnan(u) || isnan(v) || u + v + w > 1.01)
    {
        return vec3(1.,1.,1.);
    }
    //if (u + v + w > 1.)
    //{
     //   return vec3(1., 1., 1.);
    //}
    return vec3(u, v, w);
}

vec3 BaryCalc(vec2 p, E2Triangle T)
{
    float u, v, w;
    vec2 v0 = T.j - T.i, v1 = T.k - T.i, v2 = p - T.i;
    float d00 = dot(v0, v0);
    float d01 = dot(v0, v1);
    float d11 = dot(v1, v1);
    float d20 = dot(v2, v0);
    float d21 = dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
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

void mainImage(vec4 fragColor, vec2 fragCoord)
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord; // iResolution.xy;
    uv += vec2(-0.25, -0.25);
    uv *= 2.;
    // GLM PORT
    vec3 col;
    //float a1 = float(iFrame % 1500) / 1000.;// - 0.85;
    CAT CurCAT = CAT{vec2(0., 0.), vec2(1., 0.), vec2(0., 1.), -0.2, -0.2, 0.};
    float alpha = angle(CurCAT.j - CurCAT.i, CurCAT.k - CurCAT.i) + CurCAT.a_ij + CurCAT.a_ki;
    float beta = angle(CurCAT.i - CurCAT.j, CurCAT.k - CurCAT.j) + CurCAT.a_ij + CurCAT.a_jk;
    float gamma = angle(CurCAT.i - CurCAT.k, CurCAT.j - CurCAT.k) + CurCAT.a_ki + CurCAT.a_jk;
    if (CurCAT.a_ij + CurCAT.a_jk + CurCAT.a_ki == 0.)
    {
        // Euclidean case
        vec2 ti, tj, tk;
        ti = vec2(0., 0.);
        float t_i = alpha;
        float t_j = beta;
        tj = vec2(cos(t_i) + sin(t_i) / tan(t_j), 0.);
        tk = vec2(cos(t_i), sin(t_i));
        E2Triangle T = E2Triangle{ti, tj, tk};
        MobiusInfo M = SolveTrans(T.i, T.j, T.k, CurCAT.i, CurCAT.j, CurCAT.k);
        col = BaryCalc(InvMobius(uv, M), T);
    }
    else if (CurCAT.a_ij + CurCAT.a_jk + CurCAT.a_ki > 0.)
    {
        // spherical case
        vec3 ti, tj, tk;
        ti = vec3(0., 0., -1.);
   
        // law of cosines stuff to find side lenghts
        // a = length of jk, b = length of ik, c = length of ik
        float a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        float b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        float c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
        float tj1 = sqrt(1. - c * c);
        tj = vec3(tj1, 0., -c);
        float tk1 = (a - b * c) / tj1;
        float tk2 = sqrt(1. - tk1 * tk1 - b * b);
        tk = vec3(tk1, tk2, -b);
        S2Triangle T = S2Triangle{ti, tj, tk, 2.f*(CurCAT.a_ij + CurCAT.a_jk + CurCAT.a_ki)};
        MobiusInfo M = SolveTrans(ForwardSphericalProj(ti), ForwardSphericalProj(tj), ForwardSphericalProj(tk), CurCAT.i, CurCAT.j, CurCAT.k);
        col = SphericalBary(InverseSphericalProj(InvMobius(uv, M)), T);
    }
    else
    {
        // hyperbolic case
        vec3 ti, tj, tk;
        
        // law of cosines stuff to find side lenghts
        // a = length of jk, b = length of ik, c = length of ij
        float a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        float b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        float c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
        float tj1 = sqrt(c * c - 1.);
        
        float tk1 = (-a + b * c) / tj1;
        float tk2 = sqrt(-1. - tk1 * tk1 + b * b);
        
        ti = vec3(0., 0., 1.);
        tj = vec3(tj1, 0., c);
        tk = vec3(tk1, tk2, b);
        if (abs(HDist(tj, tk)+a)>1e-8) col = vec3(1,1,0);
        //col = tk;
        H2Triangle T = H2Triangle{ti, tj, tk, -2.f*(CurCAT.a_ij + CurCAT.a_jk + CurCAT.a_ki)};
        //cout << "THING:" << HyperbolicArea(T.i, T.j, T.k);
        MobiusInfo M = SolveTrans(ForwardHyperbolicProj(ti), ForwardHyperbolicProj(tj), ForwardHyperbolicProj(tk), CurCAT.i, CurCAT.j, CurCAT.k);
        col = vec3(abs(-1. - tk1 * tk1 + b * b), 0, 0);
        col = HyperbolicBary(InverseHyperbolicProj(InvMobius(uv, M)), T);
    }

    // Output to screen
    fragColor = vec4(col, 1.0);
}
int main()
{
    mainImage(vec4(0,0,0,0), vec2(0.5, 0.5));
}


