#pragma once
#include <windows.h>
#include <d3d12.h>
#include <dxgi1_4.h>
#include "d3dx12.h"

class DXRenderer
{
public:
    DXRenderer();
    ~DXRenderer();

    bool initialize(HWND hwnd, int w, int h);
    void renderFrame();
    unsigned char* getBackBufferCPU();
    void unmapBackBuffer();

    size_t getRowPitch() const { return rowPitch; }
    int getWidth() const { return width; }
    int getHeight() const { return height; }

private:
    bool createFactoryAndDevice();
    bool createCommandObjects();
    bool createSwapChainAndRTVs(HWND hwnd);
    bool createReadbackBuffer();
    bool createPipeline();
    bool createVertexBuffer();
    void waitForGpu();
    void safeRelease(IUnknown*& p);
    bool loadShader(const wchar_t* filename, ID3DBlob** blob);

    IDXGIFactory4* factory = nullptr;
    IDXGIAdapter1* adapter = nullptr;
    ID3D12Device* device = nullptr;

    ID3D12CommandQueue* queue = nullptr;
    ID3D12CommandAllocator* allocator = nullptr;
    ID3D12GraphicsCommandList* commandList = nullptr;

    IDXGISwapChain3* swapChain = nullptr;
    ID3D12DescriptorHeap* rtvHeap = nullptr;
    ID3D12Resource* renderTargets[2] = {};
    UINT rtvDescriptorSize = 0;

    ID3D12Fence* fence = nullptr;
    HANDLE fenceEvent = nullptr;
    UINT64 fenceValue = 0;

    ID3D12RootSignature* rootSignature = nullptr;
    ID3D12PipelineState* pipelineState = nullptr;

    ID3D12Resource* vertexBuffer = nullptr;
    D3D12_VERTEX_BUFFER_VIEW vertexBufferView = {};

    ID3D12Resource* readbackBuffer = nullptr;

    unsigned char* mappedPtr = nullptr;

    int width = 800;
    int height = 600;
    size_t rowPitch = 0;
    UINT backBufferIndex = 0;
};
