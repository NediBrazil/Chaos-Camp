struct Payload
{
    float4 color;
};

[shader("miss")]
void Miss(inout Payload payload)
{
    payload.color = float4(0.0, 1.0, 0.0, 1.0);
}
