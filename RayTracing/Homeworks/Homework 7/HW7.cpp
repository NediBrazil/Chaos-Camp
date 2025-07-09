#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

struct Vec3
{
    float x, y, z;
    Vec3 operator-(const Vec3 &b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec3 operator+(const Vec3 &b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    float dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec3 cross(const Vec3 &b) const
    {
        return {
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x};
    }
    Vec3 normalize() const
    {
        float len = std::sqrt(x * x + y * y + z * z);
        return {x / len, y / len, z / len};
    }
};

struct Ray
{
    Vec3 origin, direction;
};

struct Triangle
{
    Vec3 v0, v1, v2;
};

bool intersect(const Ray &ray, const Triangle &tri, float &tOut)
{
    const float EPS = 1e-6f;
    Vec3 e1 = tri.v1 - tri.v0;
    Vec3 e2 = tri.v2 - tri.v0;
    Vec3 p = ray.direction.cross(e2);
    float det = e1.dot(p);
    if (fabs(det) < EPS)
        return false;
    float invDet = 1.0f / det;
    Vec3 tvec = ray.origin - tri.v0;
    float u = tvec.dot(p) * invDet;
    if (u < 0.0f || u > 1.0f)
        return false;
    Vec3 q = tvec.cross(e1);
    float v = ray.direction.dot(q) * invDet;
    if (v < 0.0f || u + v > 1.0f)
        return false;
    float t = e2.dot(q) * invDet;
    if (t > EPS)
    {
        tOut = t;
        return true;
    }
    return false;
}

void savePPM(const std::string &fn,
             const std::vector<std::vector<Vec3>> &img,
             int w, int h)
{
    std::ofstream ofs(fn);
    ofs << "P3\n"
        << w << " " << h << "\n255\n";
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            Vec3 c = img[y][x];
            int r = std::min(255, int(c.x * 255));
            int g = std::min(255, int(c.y * 255));
            int b = std::min(255, int(c.z * 255));
            ofs << r << " " << g << " " << b << "\n";
        }
    }
}

int main()
{
    using namespace rapidjson;
    std::ifstream ifs("scene0.crtscene");
    if (!ifs.is_open())
    {
        std::cerr << "Error opening scene file.\n";
        return 1;
    }
    IStreamWrapper isw(ifs);
    Document doc;
    doc.ParseStream(isw);
    if (doc.HasParseError())
    {
        std::cerr << "Error parsing JSON.\n";
        return 1;
    }
    const auto &settings = doc["settings"];
    const auto &imgSet = settings["image_settings"];
    int width = imgSet["width"].GetInt();
    int height = imgSet["height"].GetInt();
    const auto &camera = doc["camera"];
    Vec3 camPos = {
        camera["position"][0].GetFloat(),
        camera["position"][1].GetFloat(),
        camera["position"][2].GetFloat()};

    if (!doc.HasMember("objects") || !doc["objects"].IsArray())
    {
        std::cerr << "'objects' missing or invalid.\n";
        return 1;
    }
    std::vector<Triangle> trisList;
    for (const auto &obj : doc["objects"].GetArray())
    {
        const auto &verts = obj["vertices"].GetArray();
        if (verts.Size() % 3 != 0)
        {
            std::cerr << "Vertices array size not multiple of 3.\n";
            return 1;
        }
        std::vector<Vec3> vlist;
        vlist.reserve(verts.Size() / 3);
        for (SizeType i = 0; i < verts.Size(); i += 3)
        {
            vlist.push_back({verts[i].GetFloat(),
                             verts[i + 1].GetFloat(),
                             verts[i + 2].GetFloat()});
        }

        const auto &tidx = obj["triangles"].GetArray();
        for (SizeType i = 0; i < tidx.Size(); i += 3)
        {
            int i0 = tidx[i].GetInt();
            int i1 = tidx[i + 1].GetInt();
            int i2 = tidx[i + 2].GetInt();
            trisList.push_back({vlist[i0], vlist[i1], vlist[i2]});
        }
    }

    std::vector<std::vector<Vec3>> image(height, std::vector<Vec3>(width, {0, 0, 0}));
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float u = (x + 0.5f) / float(width);
            float v = (y + 0.5f) / float(height);
            Vec3 dir = Vec3{(u - 0.5f) * 2, (v - 0.5f) * 2, -1}.normalize();
            Ray ray{camPos, dir};

            float bestT = std::numeric_limits<float>::max();
            bool hit = false;
            for (auto &T : trisList)
            {
                float t;
                if (intersect(ray, T, t) && t < bestT)
                {
                    bestT = t;
                    hit = true;
                }
            }
            if (hit)
            {
                image[y][x] = {1.0f, 0.5f, 0.2f};
            }
        }
    }

    savePPM("output.ppm", image, width, height);
    return 0;
}
