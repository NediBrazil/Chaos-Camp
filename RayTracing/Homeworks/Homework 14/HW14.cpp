#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdio>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <algorithm>
#include <condition_variable>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/en.h"

const float PI = 3.14159265f;
const int MAX_BOUNCES = 1;
int BUCKET_SIZE;
const int MAX_THREADS = std::thread::hardware_concurrency();

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

struct AABB
{
    Vector3 min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    Vector3 max = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
    void expand(const Vector3 &v)
    {
        min.x = std::min(min.x, v.x);
        min.y = std::min(min.y, v.y);
        min.z = std::min(min.z, v.z);
        max.x = std::max(max.x, v.x);
        max.y = std::max(max.y, v.y);
        max.z = std::max(max.z, v.z);
    }
    bool intersect(const Vector3 &orig, const Vector3 &dir, float &tmin, float &tmax) const
    {
        tmin = 0.0f;
        tmax = std::numeric_limits<float>::max();
        for (int i = 0; i < 3; ++i)
        {
            float invD = 1.0f / ((&dir.x)[i]);
            float t0 = ((&min.x)[i] - (&orig.x)[i]) * invD;
            float t1 = ((&max.x)[i] - (&orig.x)[i]) * invD;
            if (invD < 0.0f)
                std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax < tmin)
                return false;
        }
        return true;
    }
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
    float t = dist / dotNR;
    if (t < 0)
        return false;
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
    Vector3 v0 = tri.v1 - tri.v0, v1 = tri.v2 - tri.v0, v2 = P - tri.v0;
    float d00 = v0.dot(v0), d01 = v0.dot(v1), d11 = v1.dot(v1);
    float d20 = v2.dot(v0), d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return (tri.n0 * u + tri.n1 * v + tri.n2 * w).normalize();
}

struct KDNode
{
    AABB box;
    KDNode *left = nullptr, *right = nullptr;
    std::vector<const Triangle *> triangles;
    bool isLeaf() const { return left == nullptr && right == nullptr; }
};

KDNode *buildKDTree(std::vector<const Triangle *> &tris, int depth = 0)
{
    KDNode *node = new KDNode();
    for (auto *t : tris)
    {
        node->box.expand(t->v0);
        node->box.expand(t->v1);
        node->box.expand(t->v2);
    }
    if (tris.size() <= 4 || depth > 20)
    {
        node->triangles = tris;
        return node;
    }
    Vector3 size = node->box.max - node->box.min;
    int axis = (size.x > size.y && size.x > size.z) ? 0 : (size.y > size.z ? 1 : 2);
    float mid = 0.0f;
    for (auto *tri : tris)
        mid += ((&tri->v0.x)[axis] + (&tri->v1.x)[axis] + (&tri->v2.x)[axis]) / 3.0f;
    mid /= tris.size();
    std::vector<const Triangle *> left, right;
    for (auto *tri : tris)
    {
        float c = ((&tri->v0.x)[axis] + (&tri->v1.x)[axis] + (&tri->v2.x)[axis]) / 3.0f;
        ((c < mid) ? left : right).push_back(tri);
    }
    if (left.empty() || right.empty())
    {
        node->triangles = tris;
        return node;
    }
    node->left = buildKDTree(left, depth + 1);
    node->right = buildKDTree(right, depth + 1);
    return node;
}

bool traverseKD(const KDNode *node, const Vector3 &orig, const Vector3 &dir, const Triangle *&hit, float &tHit)
{
    float tmin, tmax;
    if (!node->box.intersect(orig, dir, tmin, tmax))
        return false;
    if (node->isLeaf())
    {
        bool any = false;
        for (auto *tri : node->triangles)
        {
            float t;
            if (intersectRayTriangle(orig, dir, *tri, t) && t < tHit)
            {
                tHit = t;
                hit = tri;
                any = true;
            }
        }
        return any;
    }
    bool h1 = traverseKD(node->left, orig, dir, hit, tHit);
    bool h2 = traverseKD(node->right, orig, dir, hit, tHit);
    return h1 || h2;
}

Vector3 traceRay(const Vector3 &orig, const Vector3 &dir, const KDNode *root, const std::vector<Light> &lights, const Vector3 &bg, int depth = 0)
{
    if (depth > MAX_BOUNCES)
        return bg;
    float tHit = std::numeric_limits<float>::max();
    const Triangle *hit = nullptr;
    if (!traverseKD(root, orig, dir, hit, tHit))
        return bg;
    Vector3 P = orig + dir * tHit;
    Vector3 N = hit->material.smooth ? interpolateNormal(*hit, P) : (hit->v1 - hit->v0).cross(hit->v2 - hit->v0).normalize();
    Vector3 col{0, 0, 0};
    for (auto &L : lights)
    {
        Vector3 ld = (L.position - P).normalize();
        float d2 = (L.position - P).dot(L.position - P);
        float st;
        const Triangle *dummy = nullptr;
        if (traverseKD(root, P + N * 0.001f, ld, dummy, st))
        {
            if (st < std::sqrt(d2))
                continue;
        }
        float nl = std::max(0.0f, N.dot(ld));
        float att = L.intensity / (4 * PI * d2);
        col = col + hit->material.albedo * nl * att;
    }
    if (hit->material.type == REFLECTIVE)
    {
        Vector3 rd = (dir - N * 2 * dir.dot(N)).normalize();
        Vector3 rc = traceRay(P + N * 0.001f, rd, root, lights, bg, depth + 1);
        col = col + rc * hit->material.albedo;
    }
    return {std::min(1.0f, col.x), std::min(1.0f, col.y), std::min(1.0f, col.z)};
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
               std::vector<Light> &lights,
               AABB &sceneBox)
{
    using namespace rapidjson;
    FILE *fp = fopen(path.c_str(), "r");
    if (!fp)
        return false;
    char buf[1 << 20];
    FileReadStream is(fp, buf, sizeof(buf));
    Document d;
    d.ParseStream(is);
    fclose(fp);
    if (d.HasParseError())
        return false;

    auto &b = d["settings"]["background_color"];
    bgColor = {b[0].GetFloat(), b[1].GetFloat(), b[2].GetFloat()};
    width = d["settings"]["image_settings"]["width"].GetInt();
    height = d["settings"]["image_settings"]["height"].GetInt();
    BUCKET_SIZE = d["settings"]["image_settings"]["bucket_size"].GetInt();
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
        {
            Vector3 vert = {v[i].GetFloat(), v[i + 1].GetFloat(), v[i + 2].GetFloat()};
            verts.push_back(vert);
            sceneBox.expand(vert);
        }
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
                tri.n0 = tri.n1 = tri.n2 = N;
            }
            triangles.push_back(tri);
        }
    }
    return true;
}

void renderScene(const std::string &in, const std::string &out)
{
    std::vector<Triangle> tris;
    std::vector<Material> mats;
    std::vector<Light> lights;
    Vector3 cam, bg, right, up, forward;
    AABB sceneBox;
    int w, h;
    if (!loadScene(in, tris, mats, cam, bg, w, h, right, up, forward, lights, sceneBox))
    {
        std::cerr << "Error loading " << in << "\n";
        return;
    }
    std::vector<std::vector<Vector3>> framebuffer(h, std::vector<Vector3>(w));
    std::queue<std::pair<int, int>> buckets;
    for (int y = 0; y < h; y += BUCKET_SIZE)
        for (int x = 0; x < w; x += BUCKET_SIZE)
            buckets.emplace(x, y);

    std::mutex mtx;
    std::vector<const Triangle *> triPtrs;
    for (auto &t : tris)
        triPtrs.push_back(&t);
    KDNode *kdRoot = buildKDTree(triPtrs);
    auto worker = [&]()
    {
        while (true)
        {
            std::pair<int, int> bucket;
            {
                std::lock_guard<std::mutex> lock(mtx);
                if (buckets.empty())
                    return;
                bucket = buckets.front();
                buckets.pop();
            }
            for (int y = bucket.second; y < std::min(h, bucket.second + BUCKET_SIZE); ++y)
            {
                for (int x = bucket.first; x < std::min(w, bucket.first + BUCKET_SIZE); ++x)
                {
                    float px = (x + 0.5f) / w, py = (y + 0.5f) / h;
                    float sx = (2 * px - 1) * (float(w) / h);
                    float sy = 1 - 2 * py;
                    Vector3 dir = (forward + right * sx + up * sy).normalize();
                    framebuffer[y][x] = traceRay(cam, dir, kdRoot, lights, bg, 0);
                }
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < MAX_THREADS; ++i)
        threads.emplace_back(worker);
    for (auto &t : threads)
        t.join();

    std::ofstream img(out);
    img << "P3\n"
        << w << " " << h << "\n255\n";
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
        {
            Vector3 c = framebuffer[y][x];
            img << int(c.x * 255) << " " << int(c.y * 255) << " " << int(c.z * 255) << "\n";
        }
    std::cout << "Rendered: " << out << "\n";
}

int main()
{
    renderScene("scene0.crtscene", "scene0_output.ppm");
    renderScene("scene1.crtscene", "scene1_output.ppm");
    return 0;
}