#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <cstdlib>

struct Vector3
{
    float x, y, z;

    Vector3 operator-(const Vector3 &v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vector3 operator+(const Vector3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vector3 operator*(float t) const { return {x * t, y * t, z * t}; }
    Vector3 operator/(float t) const { return {x / t, y / t, z / t}; }

    float dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }

    Vector3 cross(const Vector3 &v) const
    {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x};
    }

    Vector3 normalize() const
    {
        float mag = std::sqrt(x * x + y * y + z * z);
        return {x / mag, y / mag, z / mag};
    }
};

struct Triangle
{
    Vector3 v0, v1, v2;
};

struct Camera
{
    Vector3 pos = {0, 0, 0};
    Vector3 dir = {0, 0, -1};
    Vector3 up = {0, 1, 0};

    void pan(float degrees)
    {
        float rad = degrees * M_PI / 180.0f;
        float cosA = std::cos(rad);
        float sinA = std::sin(rad);
        dir = {dir.x * cosA - dir.z * sinA, dir.y, dir.x * sinA + dir.z * cosA};
        dir = dir.normalize();
    }

    void tilt(float degrees)
    {
        float rad = degrees * M_PI / 180.0f;
        Vector3 right = dir.cross(up).normalize();
        float cosA = std::cos(rad);
        float sinA = std::sin(rad);
        dir = (dir * cosA + up * sinA).normalize();
        up = right.cross(dir).normalize();
    }

    void truck(float amount)
    {
        Vector3 right = dir.cross(up).normalize();
        pos = pos + right * amount;
    }

    void pedestal(float amount)
    {
        pos = pos + up * amount;
    }

    void dolly(float amount)
    {
        pos = pos + dir * amount;
    }

    void roll(float degrees)
    {
        float rad = degrees * M_PI / 180.0f;
        Vector3 right = dir.cross(up).normalize();
        up = (up * std::cos(rad) + right * std::sin(rad)).normalize();
    }
};

bool intersectRayTriangle(const Vector3 &origin, const Vector3 &dir, const Triangle &tri, float &tOut)
{
    Vector3 e0 = tri.v1 - tri.v0;
    Vector3 e1 = tri.v2 - tri.v0;
    Vector3 normal = e0.cross(e1).normalize();
    float dotNR = normal.dot(dir);
    if (std::abs(dotNR) < 1e-5)
        return false;
    float dist = normal.dot(tri.v0 - origin);
    if (dist / dotNR < 0)
        return false;
    float t = dist / dotNR;
    Vector3 P = origin + dir * t;
    Vector3 V0P = P - tri.v0;
    Vector3 V1P = P - tri.v1;
    Vector3 V2P = P - tri.v2;
    Vector3 E0 = tri.v1 - tri.v0;
    Vector3 E1 = tri.v2 - tri.v1;
    Vector3 E2 = tri.v0 - tri.v2;
    if (normal.dot(E0.cross(V0P)) < 0)
        return false;
    if (normal.dot(E1.cross(V1P)) < 0)
        return false;
    if (normal.dot(E2.cross(V2P)) < 0)
        return false;
    tOut = t;
    return true;
}

void renderFrame(const std::string &filename, const Camera &cam, const std::vector<Triangle> &tris, int width, int height)
{
    std::ofstream image(filename);
    image << "P3\n"
          << width << " " << height << "\n255\n";
    float aspect = float(width) / height;
    Vector3 right = cam.dir.cross(cam.up).normalize();
    Vector3 up = right.cross(cam.dir).normalize();

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            float u = (2.0f * x / width - 1.0f) * aspect;
            float v = 1.0f - 2.0f * y / height;
            Vector3 dir = (cam.dir + right * u + up * v).normalize();
            float minT = std::numeric_limits<float>::max();
            bool hit = false;
            for (const auto &tri : tris)
            {
                float t;
                if (intersectRayTriangle(cam.pos, dir, tri, t))
                {
                    if (t < minT)
                    {
                        minT = t;
                        hit = true;
                    }
                }
            }
            if (hit)
                image << "255 255 255\n";
            else
                image << "0 0 0\n";
        }
    }

    image.close();
}

int main()
{
    int width = 512;
    int height = 512;
    system("mkdir frames");

    std::vector<Triangle> tris = {
        {{-1.75, -1.75, -3}, {1.75, -1.75, -3}, {0, 1.75, -3}}};

    Camera cam;

    cam.pan(30);
    renderFrame("task1_pan.ppm", cam, tris, width, height);

    cam = {};
    cam.pos = {0, 0, 2};
    cam.dir = {0, 0, -1};
    renderFrame("task2_shifted.ppm", cam, tris, width, height);

    cam = {};
    renderFrame("task3_before.ppm", cam, tris, width, height);
    cam.dolly(-1.5f);
    renderFrame("task3_after.ppm", cam, tris, width, height);

    cam = {};
    renderFrame("task4_before.ppm", cam, tris, width, height);
    cam.pan(15);
    cam.tilt(-10);
    cam.truck(1);
    renderFrame("task4_after.ppm", cam, tris, width, height);

    cam = {};
    for (int i = 0; i < 72; i++)
    {
        cam.pan(5);
        std::string fname = "frames/frame_" + std::to_string(i) + ".ppm";
        renderFrame(fname, cam, tris, width, height);
    }

    return 0;
}
