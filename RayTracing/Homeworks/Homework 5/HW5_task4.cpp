#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

struct Vector3
{
    float x, y, z;

    Vector3 operator-(const Vector3 &v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vector3 operator+(const Vector3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vector3 operator*(float t) const { return {x * t, y * t, z * t}; }

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
    int r, g, b;
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

int main()
{
    const int width = 512;
    const int height = 512;
    std::ofstream image("task4_output.ppm");
    image << "P3\n"
          << width << " " << height << "\n255\n";

    Vector3 camPos{-2, 1, 0};
    float aspectRatio = float(width) / height;

    std::vector<Triangle> triangles = {
        {{-1.75, -1.75, -3}, {-1.75, -1.75, -5}, {0, 1.75, -4}, 0, 255, 0},
        {{1.75, -1.75, -3}, {1.75, -1.75, -5}, {0, 1.75, -4}, 0, 0, 255},
        {{-1.75, -1.75, -3}, {1.75, -1.75, -3}, {0, 1.75, -4}, 255, 0, 0}};

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float px = (x + 0.5f) / width;
            float py = (y + 0.5f) / height;

            float screenX = (2.0f * px - 1.0f) * aspectRatio;
            float screenY = 1.0f - 2.0f * py;

            Vector3 dir = Vector3{screenX, screenY, -1}.normalize();

            float closestT = std::numeric_limits<float>::max();
            int r = 0, g = 0, b = 0;
            for (const auto &tri : triangles)
            {
                float t;
                if (intersectRayTriangle(camPos, dir, tri, t))
                {
                    if (t < closestT)
                    {
                        closestT = t;
                        r = tri.r;
                        g = tri.g;
                        b = tri.b;
                    }
                }
            }

            image << r << " " << g << " " << b << "\n";
        }
    }

    image.close();
    return 0;
}
