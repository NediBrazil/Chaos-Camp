#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>

const float PI = 3.14159265f;

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

struct Camera
{
    Vector3 position = {0, 0, 0};
    Vector3 direction = {0, 0, -1};
    Vector3 up = {0, 1, 0};

    void dolly(float amount)
    {
        position = position + direction.normalize() * amount;
    }

    void truck(float amount)
    {
        Vector3 right = direction.cross(up).normalize();
        position = position + right * amount;
    }

    void pedestal(float amount)
    {
        position = position + up.normalize() * amount;
    }

    void pan(float degrees)
    {
        float radians = degrees * PI / 180.0f;
        float cosA = std::cos(radians);
        float sinA = std::sin(radians);
        Vector3 newDir = {
            direction.x * cosA + direction.z * sinA,
            direction.y,
            -direction.x * sinA + direction.z * cosA};
        direction = newDir.normalize();
    }

    void tilt(float degrees)
    {
        float radians = degrees * PI / 180.0f;
        Vector3 right = direction.cross(up).normalize();
        float cosA = std::cos(radians);
        float sinA = std::sin(radians);

        Vector3 newDir = {
            direction.x * cosA + up.x * sinA,
            direction.y * cosA + up.y * sinA,
            direction.z * cosA + up.z * sinA};
        direction = newDir.normalize();
    }

    void roll(float degrees)
    {
        float radians = degrees * PI / 180.0f;
        Vector3 right = direction.cross(up).normalize();
        float cosA = std::cos(radians);
        float sinA = std::sin(radians);

        Vector3 newUp = {
            up.x * cosA + right.x * sinA,
            up.y * cosA + right.y * sinA,
            up.z * cosA + right.z * sinA};
        up = newUp.normalize();
    }
};

void renderImage(const Camera &cam, const std::vector<Triangle> &triangles, const std::string &filename)
{
    const int width = 512, height = 512;
    float aspectRatio = float(width) / height;
    std::ofstream image(filename);
    image << "P3\n"
          << width << " " << height << "\n255\n";

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float px = (x + 0.5f) / width;
            float py = (y + 0.5f) / height;
            float screenX = (2.0f * px - 1.0f) * aspectRatio;
            float screenY = 1.0f - 2.0f * py;

            Vector3 right = cam.direction.cross(cam.up).normalize();
            Vector3 imagePoint = (cam.direction + right * screenX + cam.up * screenY).normalize();

            float closestT = std::numeric_limits<float>::max();
            int r = 0, g = 0, b = 0;
            for (const auto &tri : triangles)
            {
                float t;
                if (intersectRayTriangle(cam.position, imagePoint, tri, t))
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
}

int main()
{
    std::vector<Triangle> triangles = {
        {{-1.75f, -1.75f, -3}, {1.75f, -1.75f, -3}, {0, 1.75f, -3}, 255, 0, 0}};

    Camera cam;

    cam.pan(30);
    renderImage(cam, triangles, "task1_pan_30.ppm");

    cam.position = {0, 0, 2};
    cam.direction = {0, 0, -1};
    renderImage(cam, triangles, "task2_offset_cam.ppm");

    renderImage(cam, triangles, "task3_before.ppm");
    cam.dolly(-1.0f);
    cam.truck(0.5f);
    cam.tilt(15);
    renderImage(cam, triangles, "task3_after.ppm");

    renderImage(cam, triangles, "task4_before.ppm");
    cam.pan(15);
    cam.pedestal(0.3f);
    cam.roll(10);
    renderImage(cam, triangles, "task4_after.ppm");

    for (int i = 0; i < 72; ++i)
    {
        std::ostringstream fname;
        fname << "task5_frame_" << std::setw(2) << std::setfill('0') << i << ".ppm";
        renderImage(cam, triangles, fname.str());
        cam.pan(5.0f);
    }

    return 0;
}
