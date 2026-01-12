RWTexture2D<float4> Output : register(u0);
RaytracingAccelerationStructure Scene : register(t0);

struct Payload
{
    float4 color;
};

[shader("raygeneration")]
void RayGen()
{
    uint2 id = DispatchRaysIndex().xy;
    uint2 dim = DispatchRaysDimensions().xy;

    float2 uv = (id + 0.5) / dim;
    float2 d = uv * 2.0 - 1.0;

    RayDesc ray;
    ray.Origin = float3(0, 0, -1);
    ray.Direction = normalize(float3(d.x, -d.y, 1));
    ray.TMin = 0.001;
    ray.TMax = 10000.0;

    Payload payload;
    payload.color = float4(0, 0, 0, 1);

    TraceRay(
        Scene,
        RAY_FLAG_NONE,
        0xFF,
        0,
        1,
        0,
        ray,
        payload
    );

    Output[id] = payload.color;
}