struct Payload
{
    float4 color;
};

cbuffer LightCB : register(b1)
{
    float3 lightPos;
    float  lightIntensity;
    float3 lightColor;
    float  pad;
};

[shader("closesthit")]
void ClosestHit(inout Payload payload, in BuiltInTriangleIntersectionAttributes attr)
{
    float3 hitPos = WorldRayOrigin() + RayTCurrent() * WorldRayDirection();

    // Stable fake normal from barycentrics
    float3 N;
    N.xy = attr.barycentrics;
    N.z  = 1.0 - N.x - N.y;
    N = normalize(N);

    float3 L = normalize(lightPos - hitPos);

    float diff = saturate(dot(N, L));

    uint id = InstanceID();

    float3 albedo = (id == 0)
        ? float3(1.0, 0.9, 0.2) 
        : float3(0.8, 0.2, 0.2);  

    float3 color = albedo * diff * lightColor * lightIntensity;

    payload.color = float4(color, 1.0);
}
