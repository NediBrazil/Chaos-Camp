#include "DXRenderer.h"
#include <iostream>

#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d12.lib")

DXRenderer::DXRenderer() {}
DXRenderer::~DXRenderer()
{
    if (fenceEvent) CloseHandle(fenceEvent);
    if (mappedPtr && readbackBuffer) readbackBuffer->Unmap(0, nullptr);
    safeRelease((IUnknown*&)fence);
    safeRelease((IUnknown*&)readbackBuffer);
    safeRelease((IUnknown*&)renderTargets[0]);
    safeRelease((IUnknown*&)renderTargets[1]);
    safeRelease((IUnknown*&)rtvHeap);
    safeRelease((IUnknown*&)swapChain);
    safeRelease((IUnknown*&)commandList);
    safeRelease((IUnknown*&)allocator);
    safeRelease((IUnknown*&)queue);
    safeRelease((IUnknown*&)device);
    safeRelease((IUnknown*&)adapter);
    safeRelease((IUnknown*&)factory);
}

void DXRenderer::safeRelease(IUnknown*& p)
{
    if (p) { p->Release(); p = nullptr; }
}

bool DXRenderer::createFactoryAndDevice()
{
    if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&factory)))) return false;

    for (UINT i = 0;; ++i)
    {
        IDXGIAdapter1* a = nullptr;
        if (factory->EnumAdapters1(i, &a) == DXGI_ERROR_NOT_FOUND) break;
        DXGI_ADAPTER_DESC1 desc;
        a->GetDesc1(&desc);
        if (!(desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE))
        {
            if (SUCCEEDED(D3D12CreateDevice(a, D3D_FEATURE_LEVEL_11_0, __uuidof(ID3D12Device), nullptr)))
            {
                adapter = a;
                break;
            }
        }
        a->Release();
    }

    if (!adapter) return false;
    if (FAILED(D3D12CreateDevice(adapter, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&device)))) return false;
    return true;
}

bool DXRenderer::createCommandObjects()
{
    D3D12_COMMAND_QUEUE_DESC qd = {};
    qd.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    if (FAILED(device->CreateCommandQueue(&qd, IID_PPV_ARGS(&queue)))) return false;
    if (FAILED(device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&allocator)))) return false;
    if (FAILED(device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, allocator, nullptr, IID_PPV_ARGS(&commandList)))) return false;
    if (FAILED(device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence)))) return false;
    fenceValue = 1;
    fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    return fenceEvent != nullptr;
}

bool DXRenderer::createSwapChainAndRTVs(HWND hwnd)
{
    D3D12_DESCRIPTOR_HEAP_DESC hd = {};
    hd.NumDescriptors = 2;
    hd.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
    if (FAILED(device->CreateDescriptorHeap(&hd, IID_PPV_ARGS(&rtvHeap)))) return false;
    rtvDescriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);

    DXGI_SWAP_CHAIN_DESC1 sc = {};
    sc.Width = width;
    sc.Height = height;
    sc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sc.Stereo = FALSE;
    sc.SampleDesc.Count = 1;
    sc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sc.BufferCount = 2;
    sc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
    sc.Scaling = DXGI_SCALING_STRETCH;

    IDXGISwapChain1* temp = nullptr;
    if (FAILED(factory->CreateSwapChainForHwnd(queue, hwnd, &sc, nullptr, nullptr, &temp))) return false;
    if (FAILED(temp->QueryInterface(IID_PPV_ARGS(&swapChain)))) { temp->Release(); return false; }
    temp->Release();

    backBufferIndex = swapChain->GetCurrentBackBufferIndex();

    D3D12_CPU_DESCRIPTOR_HANDLE h = rtvHeap->GetCPUDescriptorHandleForHeapStart();
    for (UINT i = 0; i < 2; ++i)
    {
        if (FAILED(swapChain->GetBuffer(i, IID_PPV_ARGS(&renderTargets[i])))) return false;
        device->CreateRenderTargetView(renderTargets[i], nullptr, h);
        h.ptr += rtvDescriptorSize;
    }

    return true;
}

bool DXRenderer::createReadbackBuffer()
{
    D3D12_RESOURCE_DESC rd = renderTargets[0]->GetDesc();
    D3D12_PLACED_SUBRESOURCE_FOOTPRINT fp = {};
    UINT numRows = 0;
    UINT64 totalBytes = 0;
    device->GetCopyableFootprints(&rd, 0, 1, 0, &fp, &numRows, nullptr, &totalBytes);
    rowPitch = fp.Footprint.RowPitch;

    D3D12_HEAP_PROPERTIES hp = {};
    hp.Type = D3D12_HEAP_TYPE_READBACK;

    D3D12_RESOURCE_DESC bd = {};
    bd.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    bd.Width = totalBytes;
    bd.Height = 1;
    bd.DepthOrArraySize = 1;
    bd.MipLevels = 1;
    bd.SampleDesc.Count = 1;
    bd.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    if (FAILED(device->CreateCommittedResource(&hp, D3D12_HEAP_FLAG_NONE, &bd, D3D12_RESOURCE_STATE_COPY_DEST, nullptr, IID_PPV_ARGS(&readbackBuffer)))) return false;
    return true;
}

bool DXRenderer::initialize(HWND hwnd, int w, int h)
{
    width = w;
    height = h;
    mappedPtr = nullptr;
    if (!createFactoryAndDevice()) return false;
    if (!createCommandObjects()) return false;
    if (!createSwapChainAndRTVs(hwnd)) return false;
    if (!createReadbackBuffer()) return false;
    return true;
}

void DXRenderer::waitForGpu()
{
    UINT64 v = fenceValue;
    queue->Signal(fence, v);
    fenceValue++;
    if (fence->GetCompletedValue() < v)
    {
        fence->SetEventOnCompletion(v, fenceEvent);
        WaitForSingleObject(fenceEvent, INFINITE);
    }
}

void DXRenderer::renderFrame()
{
    allocator->Reset();
    commandList->Reset(allocator, nullptr);

    float r = (frameIndex % 255) / 255.0f;
    float g = ((frameIndex * 2) % 255) / 255.0f;
    float b = ((frameIndex * 3) % 255) / 255.0f;
    frameIndex++;

    D3D12_RESOURCE_BARRIER b1 = {};
    b1.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    b1.Transition.pResource = renderTargets[backBufferIndex];
    b1.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
    b1.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
    b1.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
    commandList->ResourceBarrier(1, &b1);

    D3D12_CPU_DESCRIPTOR_HANDLE h = rtvHeap->GetCPUDescriptorHandleForHeapStart();
    h.ptr += backBufferIndex * rtvDescriptorSize;
    float clearColor[4] = { r, g, b, 1.0f };
    commandList->ClearRenderTargetView(h, clearColor, 0, nullptr);

    D3D12_RESOURCE_BARRIER b2 = {};
    b2.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    b2.Transition.pResource = renderTargets[backBufferIndex];
    b2.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
    b2.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
    b2.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
    commandList->ResourceBarrier(1, &b2);

    commandList->Close();
    ID3D12CommandList* lists[] = { commandList };
    queue->ExecuteCommandLists(1, lists);
    swapChain->Present(0, 0);
    waitForGpu();

    backBufferIndex = swapChain->GetCurrentBackBufferIndex();

    allocator->Reset();
    commandList->Reset(allocator, nullptr);

    D3D12_RESOURCE_BARRIER cb = {};
    cb.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    cb.Transition.pResource = renderTargets[backBufferIndex];
    cb.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
    cb.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
    cb.Transition.StateAfter = D3D12_RESOURCE_STATE_COPY_SOURCE;
    commandList->ResourceBarrier(1, &cb);

    D3D12_TEXTURE_COPY_LOCATION src = {};
    src.pResource = renderTargets[backBufferIndex];
    src.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
    src.SubresourceIndex = 0;

    D3D12_TEXTURE_COPY_LOCATION dst = {};
    dst.pResource = readbackBuffer;
    dst.Type = D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;

    D3D12_PLACED_SUBRESOURCE_FOOTPRINT fp = {};
    UINT numRows = 0;
    UINT64 total = 0;
    D3D12_RESOURCE_DESC rd = renderTargets[0]->GetDesc();
    device->GetCopyableFootprints(&rd, 0, 1, 0, &fp, &numRows, nullptr, &total);
    dst.PlacedFootprint = fp;

    commandList->CopyTextureRegion(&dst, 0, 0, 0, &src, nullptr);

    D3D12_RESOURCE_BARRIER cb2 = {};
    cb2.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    cb2.Transition.pResource = renderTargets[backBufferIndex];
    cb2.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
    cb2.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_SOURCE;
    cb2.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
    commandList->ResourceBarrier(1, &cb2);

    commandList->Close();
    ID3D12CommandList* lists2[] = { commandList };
    queue->ExecuteCommandLists(1, lists2);
    waitForGpu();
}

unsigned char* DXRenderer::getBackBufferCPU()
{
    if (!readbackBuffer) return nullptr;
    if (mappedPtr) return mappedPtr;
    unsigned char* p = nullptr;
    D3D12_RANGE r = {0, 0};
    if (FAILED(readbackBuffer->Map(0, &r, reinterpret_cast<void**>(&p)))) return nullptr;
    mappedPtr = p;
    return mappedPtr;
}

void DXRenderer::unmapBackBuffer()
{
    if (!readbackBuffer || !mappedPtr) return;
    D3D12_RANGE r = {0, 0};
    readbackBuffer->Unmap(0, &r);
    mappedPtr = nullptr;
}
