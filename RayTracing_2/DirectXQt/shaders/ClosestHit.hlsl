struct Payload
{
    float4 color;
};

[shader("closesthit")]
void ClosestHit(
    inout Payload payload,
    in BuiltInTriangleIntersectionAttributes attr)
{
    float3 n = normalize(cross(
        WorldRayDirection(),
        float3(0.0, 1.0, 0.0)
    ));

    float3 lightDir = normalize(float3(-0.4, -1.0, -0.6));
    float NdotL = saturate(dot(n, -lightDir));

    float3 baseColor = float3(0.7, 0.7, 0.7);

    payload.color = float4(baseColor * (0.15 + NdotL), 1.0);
}