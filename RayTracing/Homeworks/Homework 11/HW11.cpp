#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

template <typename T>
T clamp(T x, T a, T b) { return std::max(a, std::min(x, b)); }

struct Vec3
{
    float x, y, z;
    Vec3 operator+(Vec3 b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator-(Vec3 b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    float dot(Vec3 b) const { return x * b.x + y * b.y + z * b.z; }
    Vec3 cross(Vec3 b) const { return {y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x}; }
    Vec3 normalize() const
    {
        float l = std::sqrt(x * x + y * y + z * z);
        return {x / l, y / l, z / l};
    }
};

struct Ray
{
    Vec3 o, d;
};

struct Material
{
    std::string type;
    Vec3 albedo;
    float ior;
    bool smooth;
};

struct Triangle
{
    Vec3 v0, v1, v2, n;
    int m;
};

struct Light
{
    Vec3 p;
    float i;
};

bool intersect(const Ray &r, const Triangle &t, float &u, float &v, float &dist)
{
    Vec3 e1 = t.v1 - t.v0, e2 = t.v2 - t.v0;
    Vec3 h = r.d.cross(e2);
    float a = e1.dot(h);
    if (std::fabs(a) < 1e-6f)
        return false;
    float f = 1.0f / a;
    Vec3 s = r.o - t.v0;
    u = f * s.dot(h);
    if (u < 0 || u > 1)
        return false;
    Vec3 q = s.cross(e1);
    v = f * r.d.dot(q);
    if (v < 0 || u + v > 1)
        return false;
    dist = f * e2.dot(q);
    return dist > 1e-6f;
}

Vec3 trace(const Ray &r, const std::vector<Triangle> &T, const std::vector<Material> &M, const Vec3 &bg, const std::vector<Light> &L, int d = 0)
{
    if (d > 5)
        return bg;
    float best = 1e9, u, v;
    int hit = -1;
    for (int i = 0; i < (int)T.size(); i++)
    {
        float t, u2, v2;
        if (intersect(r, T[i], u2, v2, t) && t < best)
        {
            best = t;
            u = u2;
            v = v2;
            hit = i;
        }
    }
    if (hit < 0)
        return bg;
    const Triangle &t = T[hit];
    const Material &m = M[t.m];
    Vec3 P = r.o + r.d * best;
    Vec3 N = t.n;
    if (r.d.dot(N) > 0)
        N = N * -1.f;
    if (m.type == "refractive")
    {
        float kr = 0.0f;
        Vec3 reflectDir = {r.d.x - 2 * N.x * r.d.dot(N), r.d.y - 2 * N.y * r.d.dot(N), r.d.z - 2 * N.z * r.d.dot(N)};
        Ray rr = {P + N * 1e-4f, reflectDir.normalize()};
        Vec3 rc = trace(rr, T, M, bg, L, d + 1);
        return rc;
    }
    Vec3 col = {0, 0, 0};
    for (auto &l : L)
    {
        Vec3 Ld = l.p - P;
        float dist2 = Ld.dot(Ld);
        Vec3 Ln = Ld.normalize();
        float diff = std::max(0.f, N.dot(Ln));
        float att = l.i / dist2 * 0.2f;
        col = col + Vec3{m.albedo.x * diff * att, m.albedo.y * diff * att, m.albedo.z * diff * att};
    }
    return {clamp(col.x, 0.f, 1.f), clamp(col.y, 0.f, 1.f), clamp(col.z, 0.f, 1.f)};
}

void save(const std::string &f, const std::vector<std::vector<Vec3>> &img)
{
    std::ofstream s(f);
    int h = img.size(), w = img[0].size();
    s << "P3\n"
      << w << " " << h << "\n255\n";
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            Vec3 c = img[y][x];
            s << int(c.x * 255) << " " << int(c.y * 255) << " " << int(c.z * 255) << "\n";
        }
}

int main(int c, char **v)
{
    using namespace rapidjson;
    if (c < 2)
        return 1;
    std::ifstream f(v[1]);
    IStreamWrapper is(f);
    Document d;
    d.ParseStream(is);

    auto &cam = d["camera"];
    Vec3 camP = {cam["position"][0].GetFloat(), cam["position"][1].GetFloat(), cam["position"][2].GetFloat()};
    int W = d["settings"]["image_settings"]["width"].GetInt();
    int H = d["settings"]["image_settings"]["height"].GetInt();
    Vec3 bg = {d["settings"]["background_color"][0].GetFloat(),
               d["settings"]["background_color"][1].GetFloat(),
               d["settings"]["background_color"][2].GetFloat()};

    std::vector<Light> lights;
    for (auto &L : d["lights"].GetArray())
        lights.push_back({{L["position"][0].GetFloat(), L["position"][1].GetFloat(), L["position"][2].GetFloat()}, L["intensity"].GetFloat()});

    std::vector<Material> mats;
    for (auto &M : d["materials"].GetArray())
    {
        mats.push_back({M["type"].GetString(),
                        {M.HasMember("albedo") ? M["albedo"][0].GetFloat() : 1,
                         M.HasMember("albedo") ? M["albedo"][1].GetFloat() : 1,
                         M.HasMember("albedo") ? M["albedo"][2].GetFloat() : 1},
                        M.HasMember("ior") ? M["ior"].GetFloat() : 1,
                        M["smooth_shading"].GetBool()});
    }

    std::vector<Triangle> tris;
    for (auto &O : d["objects"].GetArray())
    {
        int mi = O["material_index"].GetInt();
        auto &V = O["vertices"];
        std::vector<Vec3> v;
        for (int i = 0; i < (int)V.Size(); i += 3)
            v.push_back({V[i].GetFloat(), V[i + 1].GetFloat(), V[i + 2].GetFloat()});
        auto &T = O["triangles"];
        for (int i = 0; i < (int)T.Size(); i += 3)
        {
            Triangle t;
            t.v0 = v[T[i].GetInt()];
            t.v1 = v[T[i + 1].GetInt()];
            t.v2 = v[T[i + 2].GetInt()];
            t.n = (t.v1 - t.v0).cross(t.v2 - t.v0).normalize();
            t.m = mi;
            tris.push_back(t);
        }
    }

    std::vector<std::vector<Vec3>> img(H, std::vector<Vec3>(W));
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
        {
            float u = (x + 0.5f) / W * 2 - 1, v = 1 - (y + 0.5f) / H * 2;
            Vec3 dir = {u, v, -1};
            img[y][x] = trace({camP, dir.normalize()}, tris, mats, bg, lights);
        }

    save(std::string(v[1]) + ".ppm", img);
    return 0;
}
