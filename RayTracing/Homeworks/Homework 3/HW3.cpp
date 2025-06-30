#include <fstream>
#include <cmath>
#include <cstdint>

struct Vec3
{
    float x, y, z;

    Vec3 normalize() const
    {
        float len = std::sqrt(x * x + y * y + z * z);
        return {x / len, y / len, z / len};
    }
};

int main()
{
    const int width = 800;
    const int height = 600;
    const float aspectRatio = float(width) / height;

    std::ofstream fout("output.ppm", std::ios::binary);
    if (!fout)
        return -1;

    fout << "P6\n"
         << width << " " << height << "\n255\n";

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float ndc_x = (x + 0.5f) / width;
            float ndc_y = (y + 0.5f) / height;

            float screen_x = (2.0f * ndc_x - 1.0f) * aspectRatio;
            float screen_y = 1.0f - 2.0f * ndc_y;

            Vec3 rayDir = {screen_x, screen_y, -1.0f};
            rayDir = rayDir.normalize();

            uint8_t r = static_cast<uint8_t>((rayDir.x * 0.5f + 0.5f) * 255);
            uint8_t g = static_cast<uint8_t>((rayDir.y * 0.5f + 0.5f) * 255);
            uint8_t b = static_cast<uint8_t>((rayDir.z * 0.5f + 0.5f) * 255);

            fout.write(reinterpret_cast<char *>(&r), 1);
            fout.write(reinterpret_cast<char *>(&g), 1);
            fout.write(reinterpret_cast<char *>(&b), 1);
        }
    }

    fout.close();
    return 0;
}
