#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdio>
#include <string>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

using namespace rapidjson;
using namespace std;

struct Vec3
{
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }
    Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    float dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec3 cross(const Vec3 &b) const
    {
        return Vec3(
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x);
    }
    Vec3 normalize() const
    {
        float len = sqrt(dot(*this));
        return *this / len;
    }
};

struct Ray
{
    Vec3 origin, dir;
    Ray(const Vec3 &o, const Vec3 &d) : origin(o), dir(d) {}
};

struct Triangle
{
    Vec3 v0, v1, v2;
};

bool intersectTriangle(const Ray &ray, const Vec3 &v0, const Vec3 &v1, const Vec3 &v2,
                       float &t, float &u, float &v)
{
    const float EPSILON = 1e-6;
    Vec3 edge1 = v1 - v0;
    Vec3 edge2 = v2 - v0;
    Vec3 h = ray.dir.cross(edge2);
    float a = edge1.dot(h);
    if (abs(a) < EPSILON)
        return false;
    float f = 1.0f / a;
    Vec3 s = ray.origin - v0;
    u = f * s.dot(h);
    if (u < 0.0 || u > 1.0)
        return false;
    Vec3 q = s.cross(edge1);
    v = f * ray.dir.dot(q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    t = f * edge2.dot(q);
    return t > EPSILON;
}

void writePPM(const string &filename, const vector<Vec3> &framebuffer, int width, int height)
{
    ofstream ofs(filename);
    ofs << "P3\n"
        << width << " " << height << "\n255\n";
    for (const Vec3 &color : framebuffer)
    {
        int r = min(255, max(0, int(color.x * 255)));
        int g = min(255, max(0, int(color.y * 255)));
        int b = min(255, max(0, int(color.z * 255)));
        ofs << r << " " << g << " " << b << "\n";
    }
    ofs.close();
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Usage: ./HW9_task1 scene0.crtscene\n";
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (!fp)
    {
        cerr << "Cannot open scene file.\n";
        return 1;
    }

    char buffer[65536];
    FileReadStream is(fp, buffer, sizeof(buffer));
    Document doc;
    doc.ParseStream(is);
    fclose(fp);

    auto &settings = doc["settings"];
    Vec3 bg(
        settings["background_color"][0].GetFloat(),
        settings["background_color"][1].GetFloat(),
        settings["background_color"][2].GetFloat());
    int width = settings["image_settings"]["width"].GetInt();
    int height = settings["image_settings"]["height"].GetInt();
    vector<Vec3> framebuffer(width * height, bg);

    auto &cam = doc["camera"];
    Vec3 camPos(
        cam["position"][0].GetFloat(),
        cam["position"][1].GetFloat(),
        cam["position"][2].GetFloat());

    vector<Triangle> triangles;
    const auto &objects = doc["objects"];
    for (SizeType o = 0; o < objects.Size(); ++o)
    {
        const auto &obj = objects[o];
        vector<Vec3> verts;
        for (SizeType i = 0; i < obj["vertices"].Size(); i += 3)
        {
            verts.emplace_back(
                obj["vertices"][i].GetFloat(),
                obj["vertices"][i + 1].GetFloat(),
                obj["vertices"][i + 2].GetFloat());
        }
        for (SizeType i = 0; i < obj["triangles"].Size(); i += 3)
        {
            int i0 = obj["triangles"][i].GetInt();
            int i1 = obj["triangles"][i + 1].GetInt();
            int i2 = obj["triangles"][i + 2].GetInt();
            triangles.push_back({verts[i0], verts[i1], verts[i2]});
        }
    }

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float u = (2 * (x + 0.5f) / float(width) - 1) * (width / float(height));
            float v = 1 - 2 * (y + 0.5f) / float(height);
            Vec3 dir = Vec3(u, v, -1).normalize();
            Ray ray(camPos, dir);

            float closestT = numeric_limits<float>::max();
            Vec3 finalColor = bg;
            for (const Triangle &tri : triangles)
            {
                float t, bU, bV;
                if (intersectTriangle(ray, tri.v0, tri.v1, tri.v2, t, bU, bV))
                {
                    if (t < closestT)
                    {
                        closestT = t;
                        float bW = 1.0f - bU - bV;
                        finalColor = Vec3(bU, bV, bW);
                    }
                }
            }
            framebuffer[y * width + x] = finalColor;
        }
    }

    string outFile = string(argv[1]) + ".ppm";
    writePPM(outFile, framebuffer, width, height);
    cout << "Saved image: " << outFile << "\n";
    return 0;
}