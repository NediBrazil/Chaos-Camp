struct Payload { float4 color; };

cbuffer LightCB : register(b1, space0)
{
    float3 lightPos;
    float  lightIntensity;
    float3 lightColor;
    float  pad;
};

StructuredBuffer<float3> Vertices : register(t2, space0);
StructuredBuffer<uint>   Indices  : register(t3, space0);

[shader("closesthit")]
void ClosestHit(inout Payload payload, in BuiltInTriangleIntersectionAttributes attr)
{
    uint prim = PrimitiveIndex();

    uint i0 = Indices[prim * 3 + 0];
    uint i1 = Indices[prim * 3 + 1];
    uint i2 = Indices[prim * 3 + 2];

    float3 p0 = Vertices[i0];
    float3 p1 = Vertices[i1];
    float3 p2 = Vertices[i2];

    float3 N = normalize(cross(p1 - p0, p2 - p0));

    float3 hitPos = WorldRayOrigin() + RayTCurrent() * WorldRayDirection();
    float3 L = normalize(lightPos - hitPos);

    if (dot(N, L) < 0.0f)
        N = -N;

    float diff = saturate(dot(N, L));
    float ambient = 0.08f;

    uint id = InstanceID();

    float3 albedo = (id == 0)
        ? float3(1.0, 0.9, 0.2)
        : float3(0.8, 0.2, 0.2);

    float3 color = albedo * (ambient + diff * lightIntensity) * lightColor;

    payload.color = float4(color, 1.0);
}
