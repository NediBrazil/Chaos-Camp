#pragma once
#include <windows.h>
#include <d3d12.h>
#include <dxgi1_4.h>
#include <cstdint>

class DXRenderer
{
public:
    DXRenderer();
    ~DXRenderer();

    bool initialize(HWND hwnd, int w, int h);
    void renderFrame();
    unsigned char *getBackBufferCPU();
    void unmapBackBuffer();
    size_t getRowPitch() const { return rowPitch; }
    int getWidth() const { return width; }
    int getHeight() const { return height; }

private:
    bool createFactoryAndDevice();
    bool createCommandObjects();
    bool createSwapChainAndRTVs(HWND hwnd);
    bool createReadbackBuffer();
    void waitForGpu();
    void safeRelease(IUnknown *&p);

    IDXGIFactory4 *factory = nullptr;
    IDXGIAdapter1 *adapter = nullptr;
    ID3D12Device *device = nullptr;

    ID3D12CommandQueue *queue = nullptr;
    ID3D12CommandAllocator *allocator = nullptr;
    ID3D12GraphicsCommandList *commandList = nullptr;

    IDXGISwapChain3 *swapChain = nullptr;
    ID3D12DescriptorHeap *rtvHeap = nullptr;
    ID3D12Resource *renderTargets[2] = {nullptr, nullptr};
    UINT rtvDescriptorSize = 0;

    ID3D12Resource *readbackBuffer = nullptr;

    ID3D12Fence *fence = nullptr;
    HANDLE fenceEvent = nullptr;
    UINT64 fenceValue = 0;

    int width = 800;
    int height = 600;
    size_t rowPitch = 0;
    UINT backBufferIndex = 0;
    int frameIndex = 0;

    unsigned char *mappedPtr = nullptr;
};
