#include <fstream>
int main(int argc, const char *argv[])
{
    const int W = 1920;
    const int H = 1080;

    std::ofstream fout("output.ppm");
    if (fout.fail())
        return -1;

    fout << "P6\n";
    fout << W << " " << H << "\n";
    fout << "255\n";

    for (int y = 0; y < H; y++)
    {
        for (int x = 0; x < W; x++)
        {
            fout << uint8_t(rand() % 256);
            fout << uint8_t(rand() % 256);
            fout << uint8_t(rand() % 256);
        }
    }
    fout.close();
    return 0;
}