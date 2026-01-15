RWTexture2D<float4> Output : register(u0);
RaytracingAccelerationStructure Scene : register(t0);

cbuffer CameraCB : register(b0)
{
    float3 camPos;
    float  pad0;
    float3 camForward;
    float  pad1;
    float3 camRight;
    float  pad2;
    float3 camUp;
    float  fov;
};

struct Payload
{
    float4 color;
};

[shader("raygeneration")]
void RayGen()
{
    uint2 id  = DispatchRaysIndex().xy;
    uint2 dim = DispatchRaysDimensions().xy;

    float2 uv = (id + 0.5) / dim;
    float2 ndc = uv * 2.0 - 1.0;

    float aspect = (float)dim.x / (float)dim.y;
    float fovScale = tan(fov * 0.5);

    float3 dir =
        camForward +
        ndc.x * aspect * fovScale * camRight -
        ndc.y * fovScale * camUp;

    RayDesc ray;
    ray.Origin = camPos;
    ray.Direction = normalize(dir);
    ray.TMin = 0.001;
    ray.TMax = 10000.0;

    Payload payload;
    payload.color = float4(0.1, 0.1, 0.1, 1.0);

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