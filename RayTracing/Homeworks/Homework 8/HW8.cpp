#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

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

bool loadScene(const std::string &path,
               std::vector<Triangle> &trianglesOut,
               Vector3 &camPosOut,
               Vector3 &bgColorOut,
               int &widthOut,
               int &heightOut,
               Vector3 &rightOut,
               Vector3 &upOut,
               Vector3 &forwardOut,
               std::vector<Light> &lightsOut)
{
    using namespace rapidjson;

    FILE *fp = fopen(path.c_str(), "r");
    if (!fp)
        return false;

    char buffer[262144];
    FileReadStream is(fp, buffer, sizeof(buffer));
    fclose(fp);

    Document doc;
    doc.ParseStream(is);
    if (doc.HasParseError())
        return false;

    auto &bg = doc["settings"]["background_color"];
    bgColorOut = {bg[0].GetFloat(), bg[1].GetFloat(), bg[2].GetFloat()};

    widthOut = doc["settings"]["image_settings"]["width"].GetInt();
    heightOut = doc["settings"]["image_settings"]["height"].GetInt();

    auto &pos = doc["camera"]["position"];
    camPosOut = {pos[0].GetFloat(), pos[1].GetFloat(), pos[2].GetFloat()};

    auto &matrix = doc["camera"]["matrix"];
    rightOut = {matrix[0].GetFloat(), matrix[1].GetFloat(), matrix[2].GetFloat()};
    upOut = {matrix[3].GetFloat(), matrix[4].GetFloat(), matrix[5].GetFloat()};
    forwardOut = {matrix[6].GetFloat(), matrix[7].GetFloat(), matrix[8].GetFloat()};
    forwardOut = forwardOut * -1;

    auto &objects = doc["objects"];
    for (auto &obj : objects.GetArray())
    {
        auto &verts = obj["vertices"];
        auto &tris = obj["triangles"];

        std::vector<Vector3> vertices;
        for (SizeType i = 0; i < verts.Size(); i += 3)
            vertices.push_back({verts[i].GetFloat(), verts[i + 1].GetFloat(), verts[i + 2].GetFloat()});

        for (SizeType i = 0; i < tris.Size(); i += 3)
            trianglesOut.push_back({vertices[tris[i].GetInt()],
                                    vertices[tris[i + 1].GetInt()],
                                    vertices[tris[i + 2].GetInt()],
                                    255, 255, 255});
    }

    if (doc.HasMember("lights"))
    {
        for (auto &light : doc["lights"].GetArray())
        {
            auto &pos = light["position"];
            float intensity = light["intensity"].GetFloat();
            lightsOut.push_back({{pos[0].GetFloat(), pos[1].GetFloat(), pos[2].GetFloat()}, intensity});
        }
    }

    return true;
}

void renderScene(const std::string &jsonPath, const std::string &outputImage)
{
    std::vector<Triangle> triangles;
    std::vector<Light> lights;
    Vector3 camPos, bgColor, right, up, forward;
    int width = 0, height = 0;

    if (!loadScene(jsonPath, triangles, camPos, bgColor, width, height, right, up, forward, lights))
    {
        std::cerr << "Failed to load scene: " << jsonPath << "\n";
        return;
    }

    std::ofstream image(outputImage);
    image << "P3\n"
          << width << " " << height << "\n255\n";
    float aspectRatio = float(width) / height;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float px = (x + 0.5f) / width;
            float py = (y + 0.5f) / height;
            float screenX = (2.0f * px - 1.0f) * aspectRatio;
            float screenY = 1.0f - 2.0f * py;

            Vector3 rayDir = (forward + right * screenX + up * screenY).normalize();

            float closestT = std::numeric_limits<float>::max();
            const Triangle *hitTri = nullptr;

            for (const auto &tri : triangles)
            {
                float t;
                if (intersectRayTriangle(camPos, rayDir, tri, t) && t < closestT)
                {
                    closestT = t;
                    hitTri = &tri;
                }
            }

            int r = int(bgColor.x * 255), g = int(bgColor.y * 255), b = int(bgColor.z * 255);

            if (hitTri)
            {
                Vector3 hitPoint = camPos + rayDir * closestT;
                Vector3 normal = (hitTri->v1 - hitTri->v0).cross(hitTri->v2 - hitTri->v0).normalize();
                Vector3 color = {0, 0, 0};
                Vector3 albedo = {1, 1, 1};

                for (const auto &light : lights)
                {
                    Vector3 lightDir = (light.position - hitPoint).normalize();

                    bool inShadow = false;
                    float distToLight = (light.position - hitPoint).dot(light.position - hitPoint);
                    for (const auto &otherTri : triangles)
                    {
                        float shadowT;
                        if (&otherTri != hitTri &&
                            intersectRayTriangle(hitPoint + normal * 0.001f, lightDir, otherTri, shadowT) &&
                            shadowT * shadowT < distToLight)
                        {
                            inShadow = true;
                            break;
                        }
                    }

                    if (!inShadow)
                    {
                        float NdotL = std::max(0.0f, normal.dot(lightDir));
                        float distance2 = (light.position - hitPoint).dot(light.position - hitPoint);
                        float attenuation = light.intensity / (4.0f * PI * distance2);

                        Vector3 lightColor = albedo * NdotL * attenuation;
                        color = color + lightColor;
                    }
                }

                color.x = std::min(1.0f, color.x);
                color.y = std::min(1.0f, color.y);
                color.z = std::min(1.0f, color.z);

                r = int(color.x * 255);
                g = int(color.y * 255);
                b = int(color.z * 255);
            }

            image << r << " " << g << " " << b << "\n";
        }
    }

    image.close();
    std::cout << "Rendered: " << outputImage << "\n";
}

int main()
{
    renderScene("scene0.crtscene", "scene0_output.ppm");
    renderScene("scene1.crtscene", "scene1_output.ppm");
    renderScene("scene2.crtscene", "scene2_output.ppm");
    renderScene("scene3.crtscene", "scene3_output.ppm");
    return 0;
}
