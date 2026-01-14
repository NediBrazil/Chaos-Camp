struct Payload
{
    float4 color;
};

[shader("miss")]
void Miss(inout Payload payload)
{
    payload.color = float4(0.02, 0.02, 0.025, 1.0);
}
