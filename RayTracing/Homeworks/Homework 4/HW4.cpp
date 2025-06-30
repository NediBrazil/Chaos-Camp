#include <iostream>
#include <cmath>

struct Vector3
{
    float x, y, z;

    Vector3 operator-(const Vector3 &v) const
    {
        return {x - v.x, y - v.y, z - v.z};
    }

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

    void print() const
    {
        std::cout << "(" << x << ", " << y << ", " << z << ")\n";
    }
};

float vectorMagnitude(const Vector3 &v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

struct Triangle
{
    Vector3 A, B, C;
};

Vector3 calculateNormal(const Triangle &tri)
{
    Vector3 U = tri.B - tri.A;
    Vector3 V = tri.C - tri.A;
    return U.cross(V).normalize();
}

float triangleArea(const Triangle &tri)
{
    Vector3 AB = tri.B - tri.A;
    Vector3 AC = tri.C - tri.A;
    return 0.5f * vectorMagnitude(AB.cross(AC));
}

int main()
{
    Vector3 A1{3.5f, 0.0f, 0.0f}, B1{1.75f, 3.5f, 0.0f};
    Vector3 cross1 = A1.cross(B1);
    std::cout << "Cross Product 1: ";
    cross1.print();

    Vector3 A2{3, -3, 1}, B2{4, 9, 3};
    Vector3 cross2 = A2.cross(B2);
    std::cout << "Cross Product 2: ";
    cross2.print();
    std::cout << "Parallelogram Area 1: " << vectorMagnitude(cross2) << "\n";

    Vector3 B3{-12, 12, -4};
    Vector3 cross3 = A2.cross(B3);
    std::cout << "Parallelogram Area 2: " << vectorMagnitude(cross3) << "\n";

    Triangle t1{{-1.75, -1.75, -3}, {1.75, -1.75, -3}, {0, 1.75, -3}};
    std::cout << "Normal 1: ";
    calculateNormal(t1).print();
    std::cout << "Area 1: " << triangleArea(t1) << "\n";

    Triangle t2{{0, 0, -1}, {1, 0, 1}, {-1, 0, 1}};
    std::cout << "Normal 2: ";
    calculateNormal(t2).print();
    std::cout << "Area 2: " << triangleArea(t2) << "\n";

    Triangle t3{{0.56, 1.11, 1.23}, {0.44, -2.368, -0.54}, {-1.56, 0.15, -1.92}};
    std::cout << "Normal 3: ";
    calculateNormal(t3).print();
    std::cout << "Area 3: " << triangleArea(t3) << "\n";

    return 0;
}
