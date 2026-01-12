struct Payload
{
    float4 color;
};

[shader("closesthit")]
void ClosestHit(inout Payload payload, in BuiltInTriangleIntersectionAttributes attr)
{
    payload.color = float4(1.0, 0.0, 0.0, 1.0);
}
