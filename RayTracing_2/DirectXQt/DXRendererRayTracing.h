#pragma once
#include <windows.h>
#include <d3d12.h>
#include <dxgi1_4.h>
#include "d3dx12.h"

class DXRendererRayTracing
{
public:
    bool initialize(ID3D12Device* baseDevice, ID3D12CommandQueue* baseQueue, int w, int h);
    void renderFrame(ID3D12GraphicsCommandList4* cmd);
    void cleanup();
    ID3D12Resource* getOutputTexture();
    ~DXRendererRayTracing();

private:
    UINT sbtRecordSize = 0;
    bool accelBuilt = false;
    bool createOutputTexture();
    bool createDescriptorHeap();
    bool createGlobalRootSignature();
    bool createStateObject();
    bool createShaderBindingTable();
    bool createTriangleBuffers();
    bool createBLAS(ID3D12GraphicsCommandList4* cmd);
    bool createTLAS(ID3D12GraphicsCommandList4* cmd);
    void createTLASSRV();

    ID3D12Device* device = nullptr;
    ID3D12Device5* device5 = nullptr;
    ID3D12CommandQueue* queue = nullptr;

    ID3D12Resource* outputTexture = nullptr;
    ID3D12DescriptorHeap* descriptorHeap = nullptr;

    ID3D12RootSignature* globalRootSignature = nullptr;
    ID3D12StateObject* stateObject = nullptr;
    ID3D12StateObjectProperties* stateProps = nullptr;

    ID3D12Resource* sbtBuffer = nullptr;

    ID3D12Resource* vertexBuffer = nullptr;
    ID3D12Resource* indexBuffer = nullptr;

    ID3D12Resource* blasBuffer = nullptr;
    ID3D12Resource* blasScratch = nullptr;

    ID3D12Resource* tlasBuffer = nullptr;
    ID3D12Resource* tlasScratch = nullptr;
    ID3D12Resource* tlasInstanceDesc = nullptr;

    int width = 0;
    int height = 0;
};
