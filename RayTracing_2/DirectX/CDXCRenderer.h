#ifndef CDXCRENDERER_H
#define CDXCRENDERER_H

#include <windows.h>
#include <wrl.h>
#include <dxgi1_4.h>
#include <d3d12.h>

using Microsoft::WRL::ComPtr;

class CDXCRenderer
{
public:
    void prepareForRendering();
    void render();

private:
    void createDevice();
    void createCommandsManagers();
    void createRenderTarget();
    void performReadbackAndPrintFirstPixel();

    ComPtr<IDXGIFactory4> dxgiFactory;
    ComPtr<IDXGIAdapter1> adapter;
    ComPtr<ID3D12Device> d3d12Device;

    ComPtr<ID3D12CommandQueue> commandQueue;
    ComPtr<ID3D12CommandAllocator> commandAllocator;
    ComPtr<ID3D12GraphicsCommandList> commandList;

    ComPtr<ID3D12Fence> fence;
    UINT64 fenceValue = 0;
    HANDLE fenceEvent = nullptr;

    ComPtr<ID3D12Resource> renderTarget;
    ComPtr<ID3D12DescriptorHeap> rtvHeap;
    D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = {};

    UINT texWidth = 800;
    UINT texHeight = 600;
    DXGI_FORMAT texFormat = DXGI_FORMAT_R8G8B8A8_UNORM;
};

#endif
