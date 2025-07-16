#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

struct Vec3
{
    float x, y, z;
    Vec3 operator+(const Vec3 &b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator-(const Vec3 &b) const { return {x - b.x, y - b.y, z - b.z}; }
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
    int material_index;
};

struct Material
{
    std::string type;
    Vec3 albedo;
    bool smooth_shading;
};

bool intersect(const Ray &ray, const Triangle &tri, float &tOut)
{
    const float EPSILON = 1e-6f;
    Vec3 edge1 = tri.v1 - tri.v0;
    Vec3 edge2 = tri.v2 - tri.v0;
    Vec3 h = ray.direction.cross(edge2);
    float a = edge1.dot(h);
    if (fabs(a) < EPSILON)
        return false;
    float f = 1.0f / a;
    Vec3 s = ray.origin - tri.v0;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f)
        return false;
    Vec3 q = s.cross(edge1);
    float v = f * ray.direction.dot(q);
    if (v < 0.0f || u + v > 1.0f)
        return false;
    float t = f * edge2.dot(q);
    if (t > EPSILON)
    {
        tOut = t;
        return true;
    }
    return false;
}

void savePPM(const std::string &filename, const std::vector<std::vector<Vec3>> &image, int width, int height)
{
    std::ofstream ofs(filename);
    ofs << "P3\n"
        << width << " " << height << "\n255\n";
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
        {
            Vec3 c = image[y][x];
            int r = std::min(255, int(c.x * 255));
            int g = std::min(255, int(c.y * 255));
            int b = std::min(255, int(c.z * 255));
            ofs << r << " " << g << " " << b << "\n";
        }
}

int main(int argc, char **argv)
{
    using namespace rapidjson;

    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " scene_file.crtscene\n";
        return 1;
    }

    std::ifstream ifs(argv[1]);
    if (!ifs.is_open())
    {
        std::cerr << "Failed to open scene file.\n";
        return 1;
    }

    IStreamWrapper isw(ifs);
    Document doc;
    doc.ParseStream(isw);
    if (doc.HasParseError())
    {
        std::cerr << "JSON parse error.\n";
        return 1;
    }

    int width = doc["settings"]["image_settings"]["width"].GetInt();
    int height = doc["settings"]["image_settings"]["height"].GetInt();

    const auto &camPosArr = doc["camera"]["position"].GetArray();
    Vec3 camPos = {camPosArr[0].GetFloat(), camPosArr[1].GetFloat(), camPosArr[2].GetFloat()};

    std::vector<Material> materials;
    for (const auto &mat : doc["materials"].GetArray())
    {
        const auto &albedo = mat["albedo"].GetArray();
        materials.push_back({mat["type"].GetString(),
                             {albedo[0].GetFloat(), albedo[1].GetFloat(), albedo[2].GetFloat()},
                             mat["smooth_shading"].GetBool()});
    }

    std::vector<Triangle> triangles;
    for (const auto &obj : doc["objects"].GetArray())
    {
        int material_index = obj["material_index"].GetInt();

        const auto &verts = obj["vertices"].GetArray();
        std::vector<Vec3> vertexList;
        for (SizeType i = 0; i < verts.Size(); i += 3)
            vertexList.push_back({verts[i].GetFloat(), verts[i + 1].GetFloat(), verts[i + 2].GetFloat()});

        const auto &indices = obj["triangles"].GetArray();
        for (SizeType i = 0; i < indices.Size(); i += 3)
        {
            Triangle tri = {
                vertexList[indices[i].GetInt()],
                vertexList[indices[i + 1].GetInt()],
                vertexList[indices[i + 2].GetInt()],
                material_index};
            triangles.push_back(tri);
        }
    }

    std::vector<std::vector<Vec3>> image(height, std::vector<Vec3>(width, {0, 0, 0}));

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float u = (x + 0.5f) / width;
            float v = (y + 0.5f) / height;
            Vec3 dir = Vec3{(u - 0.5f) * 2, (v - 0.5f) * 2, -1}.normalize();
            Ray ray{camPos, dir};

            float closestT = std::numeric_limits<float>::max();
            int hitMaterial = -1;

            for (const auto &tri : triangles)
            {
                float t;
                if (intersect(ray, tri, t) && t < closestT)
                {
                    closestT = t;
                    hitMaterial = tri.material_index;
                }
            }

            if (hitMaterial != -1)
            {
                image[y][x] = materials[hitMaterial].albedo;
            }
        }
    }

    string outFile = string(argv[1]) + ".ppm";
    writePPM(outFile, framebuffer, width, height);
    cout << "Saved image: " << outFile << "\n";
    return 0;
}
