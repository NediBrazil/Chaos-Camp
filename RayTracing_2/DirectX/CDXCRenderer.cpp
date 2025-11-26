#include "CDXCRenderer.h"
#include <iostream>
#pragma comment(lib, "d3d12.lib")
#pragma comment(lib, "dxgi.lib")

inline void PrintHResultFailure(const char *msg, HRESULT hr)
{
    std::cout << msg << " (0x" << std::hex << hr << ")\n";
}

#define CHECK_HR(hr, msg)             \
    if (FAILED(hr))                   \
    {                                 \
        PrintHResultFailure(msg, hr); \
        return;                       \
    }

void CDXCRenderer::prepareForRendering()
{
    HRESULT hr = CreateDXGIFactory1(IID_PPV_ARGS(&dxgiFactory));
    CHECK_HR(hr, "CreateDXGIFactory1 failed");

    for (UINT i = 0; dxgiFactory->EnumAdapters1(i, &adapter) != DXGI_ERROR_NOT_FOUND; ++i)
    {
        DXGI_ADAPTER_DESC1 desc;
        adapter->GetDesc1(&desc);

        if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)
        {
            adapter.Reset();
            continue;
        }

        if (SUCCEEDED(D3D12CreateDevice(adapter.Get(), D3D_FEATURE_LEVEL_11_0, __uuidof(ID3D12Device), nullptr)))
            break;

        adapter.Reset();
    }

    if (!adapter)
        return;

    createDevice();
    createCommandsManagers();
    createRenderTarget();

    d3d12Device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence));
    fenceValue = 1;
    fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
}

void CDXCRenderer::createDevice()
{
    HRESULT hr = D3D12CreateDevice(adapter.Get(), D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&d3d12Device));
    CHECK_HR(hr, "CreateDevice failed");
}

void CDXCRenderer::createCommandsManagers()
{
    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;

    HRESULT hr = d3d12Device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&commandQueue));
    CHECK_HR(hr, "CreateCommandQueue failed");

    hr = d3d12Device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&commandAllocator));
    CHECK_HR(hr, "CreateCommandAllocator failed");

    hr = d3d12Device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, commandAllocator.Get(), nullptr, IID_PPV_ARGS(&commandList));
    CHECK_HR(hr, "CreateCommandList failed");

    commandList->Close();
}

void CDXCRenderer::createRenderTarget()
{
    D3D12_HEAP_PROPERTIES heapProps = {};
    heapProps.Type = D3D12_HEAP_TYPE_DEFAULT;

    D3D12_RESOURCE_DESC texDesc = {};
    texDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
    texDesc.Width = texWidth;
    texDesc.Height = texHeight;
    texDesc.DepthOrArraySize = 1;
    texDesc.MipLevels = 1;
    texDesc.Format = texFormat;
    texDesc.SampleDesc.Count = 1;
    texDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET;

    D3D12_CLEAR_VALUE clearValue = {};
    clearValue.Format = texFormat;

    HRESULT hr = d3d12Device->CreateCommittedResource(
        &heapProps,
        D3D12_HEAP_FLAG_NONE,
        &texDesc,
        D3D12_RESOURCE_STATE_RENDER_TARGET,
        &clearValue,
        IID_PPV_ARGS(&renderTarget));
    CHECK_HR(hr, "CreateCommittedResource failed");

    D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
    heapDesc.NumDescriptors = 1;
    heapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;

    hr = d3d12Device->CreateDescriptorHeap(&heapDesc, IID_PPV_ARGS(&rtvHeap));
    CHECK_HR(hr, "CreateDescriptorHeap failed");

    rtvHandle = rtvHeap->GetCPUDescriptorHandleForHeapStart();
    d3d12Device->CreateRenderTargetView(renderTarget.Get(), nullptr, rtvHandle);
}

void CDXCRenderer::render()
{
    commandAllocator->Reset();
    commandList->Reset(commandAllocator.Get(), nullptr);

    FLOAT clearColor[] = {0.3f, 0.5f, 0.8f, 1.0f};
    commandList->OMSetRenderTargets(1, &rtvHandle, FALSE, nullptr);
    commandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);

    performReadbackAndPrintFirstPixel();

    std::cout << "Render + readback complete\n";
}

void CDXCRenderer::performReadbackAndPrintFirstPixel()
{
    D3D12_RESOURCE_BARRIER barrier = {};
    barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    barrier.Transition.pResource = renderTarget.Get();
    barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
    barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_COPY_SOURCE;

    commandList->ResourceBarrier(1, &barrier);

    D3D12_RESOURCE_DESC textureDesc = renderTarget->GetDesc();

    UINT64 totalBytes = 0;
    UINT numRows = 0;
    UINT64 rowSizeInBytes = 0;

    D3D12_PLACED_SUBRESOURCE_FOOTPRINT layout;
    UINT rows;
    UINT64 rowSize;
    d3d12Device->GetCopyableFootprints(
        &textureDesc, 0, 1, 0,
        &layout, &rows, &rowSize, &totalBytes);

    D3D12_HEAP_PROPERTIES readbackHeap = {};
    readbackHeap.Type = D3D12_HEAP_TYPE_READBACK;

    D3D12_RESOURCE_DESC bufferDesc = {};
    bufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    bufferDesc.Width = totalBytes;
    bufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    ComPtr<ID3D12Resource> readbackBuffer;
    d3d12Device->CreateCommittedResource(
        &readbackHeap,
        D3D12_HEAP_FLAG_NONE,
        &bufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&readbackBuffer));

    D3D12_TEXTURE_COPY_LOCATION src = {};
    src.pResource = renderTarget.Get();
    src.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;

    D3D12_TEXTURE_COPY_LOCATION dst = {};
    dst.pResource = readbackBuffer.Get();
    dst.Type = D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;
    dst.PlacedFootprint = layout;

    D3D12_BOX box = {};
    box.right = texWidth;
    box.bottom = texHeight;
    box.back = 1;

    commandList->CopyTextureRegion(&dst, 0, 0, 0, &src, &box);

    barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_SOURCE;
    barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
    commandList->ResourceBarrier(1, &barrier);

    commandList->Close();

    ID3D12CommandList *lists[] = {commandList.Get()};
    commandQueue->ExecuteCommandLists(1, lists);

    UINT64 signalValue = fenceValue++;
    commandQueue->Signal(fence.Get(), signalValue);

    if (fence->GetCompletedValue() < signalValue)
    {
        fence->SetEventOnCompletion(signalValue, fenceEvent);
        WaitForSingleObject(fenceEvent, INFINITE);
    }

    void *mapped = nullptr;
    D3D12_RANGE readRange = {0, (SIZE_T)totalBytes};
    readbackBuffer->Map(0, &readRange, &mapped);

    uint8_t *p = reinterpret_cast<uint8_t *>(mapped) + layout.Offset;
    std::cout << "RGBA: "
              << (int)p[0] << " "
              << (int)p[1] << " "
              << (int)p[2] << " "
              << (int)p[3] << "\n";

    D3D12_RANGE write = {0, 0};
    readbackBuffer->Unmap(0, &write);
}