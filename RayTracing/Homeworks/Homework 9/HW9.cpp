#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdio>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/en.h"

const float PI = 3.14159265f;
const int MAX_BOUNCES = 1;

struct Vector3
{
    float x, y, z;
    Vector3 operator-(const Vector3 &v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vector3 operator+(const Vector3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vector3 operator*(float t) const { return {x * t, y * t, z * t}; }
    Vector3 operator*(const Vector3 &v) const { return {x * v.x, y * v.y, z * v.z}; }
    float dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }
    Vector3 cross(const Vector3 &v) const { return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x}; }
    Vector3 normalize() const
    {
        float m = std::sqrt(x * x + y * y + z * z);
        return {x / m, y / m, z / m};
    }
};

enum MaterialType
{
    DIFFUSE,
    REFLECTIVE
};

struct Material
{
    MaterialType type;
    Vector3 albedo;
    bool smooth;
};

struct Triangle
{
    Vector3 v0, v1, v2;
    Vector3 n0, n1, n2;
    Material material;
};

struct Light
{
    Vector3 position;
    float intensity;
};

bool intersectRayTriangle(const Vector3 &origin, const Vector3 &dir, const Triangle &tri, float &tOut)
{
    Vector3 e0 = tri.v1 - tri.v0;
    Vector3 e1 = tri.v2 - tri.v0;
    Vector3 normal = e0.cross(e1).normalize();
    float dotNR = normal.dot(dir);
    if (std::abs(dotNR) < 1e-5f)
        return false;
    float dist = normal.dot(tri.v0 - origin);
    if (dist / dotNR < 0)
        return false;
    float t = dist / dotNR;
    Vector3 P = origin + dir * t;
    Vector3 V0P = P - tri.v0, V1P = P - tri.v1, V2P = P - tri.v2;
    Vector3 E0 = tri.v1 - tri.v0, E1 = tri.v2 - tri.v1, E2 = tri.v0 - tri.v2;
    if (normal.dot(E0.cross(V0P)) < 0)
        return false;
    if (normal.dot(E1.cross(V1P)) < 0)
        return false;
    if (normal.dot(E2.cross(V2P)) < 0)
        return false;
    tOut = t;
    return true;
}

Vector3 interpolateNormal(const Triangle &tri, const Vector3 &P)
{
    Vector3 v0 = tri.v1 - tri.v0;
    Vector3 v1 = tri.v2 - tri.v0;
    Vector3 v2 = P - tri.v0;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return (tri.n0 * u + tri.n1 * v + tri.n2 * w).normalize();
}

bool loadScene(const std::string &path,
               std::vector<Triangle> &triangles,
               std::vector<Material> &materials,
               Vector3 &camPos,
               Vector3 &bgColor,
               int &width,
               int &height,
               Vector3 &right,
               Vector3 &up,
               Vector3 &forward,
               std::vector<Light> &lights)
{
    using namespace rapidjson;
    FILE *fp = fopen(path.c_str(), "r");
    if (!fp)
    {
        std::cerr << "Cannot open " << path << "\n";
        return false;
    }
    char buf[1 << 20];
    FileReadStream is(fp, buf, sizeof(buf));
    Document d;
    d.ParseStream(is);
    fclose(fp);
    if (d.HasParseError())
    {
        std::cerr << "JSON parse error at offset " << d.GetErrorOffset()
                  << ": " << GetParseError_En(d.GetParseError()) << "\n";
        return false;
    }
    auto &b = d["settings"]["background_color"];
    bgColor = {b[0].GetFloat(), b[1].GetFloat(), b[2].GetFloat()};
    width = d["settings"]["image_settings"]["width"].GetInt();
    height = d["settings"]["image_settings"]["height"].GetInt();
    auto &cp = d["camera"]["position"];
    camPos = {cp[0].GetFloat(), cp[1].GetFloat(), cp[2].GetFloat()};
    auto &cm = d["camera"]["matrix"];
    right = {cm[0].GetFloat(), cm[1].GetFloat(), cm[2].GetFloat()};
    up = {cm[3].GetFloat(), cm[4].GetFloat(), cm[5].GetFloat()};
    forward = {cm[6].GetFloat(), cm[7].GetFloat(), cm[8].GetFloat()};
    forward = forward * -1;
    for (auto &lm : d["lights"].GetArray())
    {
        auto &p = lm["position"];
        lights.push_back({{p[0].GetFloat(), p[1].GetFloat(), p[2].GetFloat()},
                          lm["intensity"].GetFloat()});
    }
    for (auto &m : d["materials"].GetArray())
    {
        Material mat;
        std::string t = m["type"].GetString();
        mat.type = (t == "reflective" ? REFLECTIVE : DIFFUSE);
        auto &a = m["albedo"];
        mat.albedo = {a[0].GetFloat(), a[1].GetFloat(), a[2].GetFloat()};
        mat.smooth = m["smooth_shading"].GetBool();
        materials.push_back(mat);
    }
    for (auto &o : d["objects"].GetArray())
    {
        int mid = o["material_index"].GetInt();
        auto &v = o["vertices"];
        std::vector<Vector3> verts;
        for (int i = 0; i < v.Size(); i += 3)
            verts.push_back({v[i].GetFloat(), v[i + 1].GetFloat(), v[i + 2].GetFloat()});
        std::vector<Vector3> norms;
        if (o.HasMember("normals"))
        {
            auto &n = o["normals"];
            for (int i = 0; i < n.Size(); i += 3)
                norms.push_back({n[i].GetFloat(), n[i + 1].GetFloat(), n[i + 2].GetFloat()});
        }
        auto &tris = o["triangles"];
        for (int i = 0; i < tris.Size(); i += 3)
        {
            int i0 = tris[i].GetInt(), i1 = tris[i + 1].GetInt(), i2 = tris[i + 2].GetInt();
            Triangle tri;
            tri.v0 = verts[i0];
            tri.v1 = verts[i1];
            tri.v2 = verts[i2];
            tri.material = materials[mid];
            if (!norms.empty())
            {
                tri.n0 = norms[i0];
                tri.n1 = norms[i1];
                tri.n2 = norms[i2];
            }
            else
            {
                Vector3 N = (tri.v1 - tri.v0).cross(tri.v2 - tri.v0).normalize();
                tri.n0 = N;
                tri.n1 = N;
                tri.n2 = N;
            }
            triangles.push_back(tri);
        }
    }
    return true;
}

Vector3 traceRay(const Vector3 &orig, const Vector3 &dir,
                 const std::vector<Triangle> &trs,
                 const std::vector<Light> &lights,
                 const Vector3 &bg, int depth = 0)
{
    if (depth > MAX_BOUNCES)
        return bg;
    float best = std::numeric_limits<float>::max();
    const Triangle *hit = nullptr;
    for (auto &tri : trs)
    {
        float t;
        if (intersectRayTriangle(orig, dir, tri, t) && t < best)
        {
            best = t;
            hit = &tri;
        }
    }
    if (!hit)
        return bg;
    Vector3 P = orig + dir * best;
    Vector3 N;
    if (hit->material.smooth)
    {
        N = interpolateNormal(*hit, P);
    }
    else
    {
        N = (hit->v1 - hit->v0).cross(hit->v2 - hit->v0).normalize();
    }
    Vector3 col{0, 0, 0};
    for (auto &L : lights)
    {
        Vector3 ld = (L.position - P).normalize();
        float d2 = (L.position - P).dot(L.position - P);
        bool sh = false;
        for (auto &o : trs)
        {
            float st;
            if (&o != hit && intersectRayTriangle(P + N * 0.001f, ld, o, st))
            {
                sh = true;
                break;
            }
        }
        if (sh)
            continue;
        float nl = std::max(0.0f, N.dot(ld));
        float att = L.intensity / (4 * PI * d2);
        col = col + hit->material.albedo * nl * att;
    }
    if (hit->material.type == REFLECTIVE)
    {
        Vector3 rd = dir - N * 2 * dir.dot(N);
        rd = rd.normalize();
        Vector3 rc = traceRay(P + N * 0.001f, rd, trs, lights, bg, depth + 1);
        col = col + rc * hit->material.albedo;
    }
    return {std::min(1.0f, col.x), std::min(1.0f, col.y), std::min(1.0f, col.z)};
}

void renderScene(const std::string &in, const std::string &out)
{
    std::vector<Triangle> tris;
    std::vector<Material> mats;
    std::vector<Light> lights;
    Vector3 cam, bg, right, up, forward;
    int w, h;
    if (!loadScene(in, tris, mats, cam, bg, w, h, right, up, forward, lights))
    {
        std::cerr << "Error loading " << in << "\n";
        return;
    }
    std::ofstream img(out);
    img << "P3\n"
        << w << " " << h << "\n255\n";
    float ar = float(w) / h;
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            float px = (x + 0.5f) / w, py = (y + 0.5f) / h;
            float sx = (2 * px - 1) * ar, sy = 1 - 2 * py;
            Vector3 dir = (forward + right * sx + up * sy).normalize();
            Vector3 c = traceRay(cam, dir, tris, lights, bg, 0);
            int R = int(c.x * 255), G = int(c.y * 255), B = int(c.z * 255);
            img << R << " " << G << " " << B << "\n";
        }
    }
    std::cout << "Rendered: " << out << "\n";
}

int main()
{
    renderScene("scene2.crtscene", "scene2_output.ppm");
    renderScene("scene3.crtscene", "scene3_output.ppm");
    renderScene("scene4.crtscene", "scene4_output.ppm");
    renderScene("scene5.crtscene", "scene5_output.ppm");
    return 0;
}
