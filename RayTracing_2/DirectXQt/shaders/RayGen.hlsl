RWTexture2D<float4> Output : register(u0);

[shader("raygeneration")]
void RayGen()
{
    uint2 id = DispatchRaysIndex().xy;
    Output[id] = float4(0.0, 0.0, 1.0, 1.0);
}