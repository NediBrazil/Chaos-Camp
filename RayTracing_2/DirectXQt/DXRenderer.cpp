#include "DXRenderer.h"
#include "DXRendererRayTracing.h"
#include <d3dcompiler.h>
#include <cstring>
#include "d3dx12.h"
#include <dxgi1_6.h>
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d12.lib")
#pragma comment(lib, "d3dcompiler.lib")

struct Vertex
{
    float x;
    float y;
};


inline void PrintHR(const char* msg, HRESULT hr)
{
    if (FAILED(hr))
    {
        char buf[512];
        sprintf_s(buf, "%s\n", msg, static_cast<unsigned int>(hr));
        OutputDebugStringA(buf);
    }
}


DXRenderer::DXRenderer() {}

DXRenderer::~DXRenderer()
{
    if (fenceEvent) CloseHandle(fenceEvent);
    if (mappedPtr && readbackBuffer)
        readbackBuffer->Unmap(0, nullptr);

    safeRelease((IUnknown*&)vertexBuffer);
    safeRelease((IUnknown*&)pipelineState);
    safeRelease((IUnknown*&)rootSignature);
    safeRelease((IUnknown*&)readbackBuffer);
    safeRelease((IUnknown*&)renderTargets[0]);
    safeRelease((IUnknown*&)renderTargets[1]);
    safeRelease((IUnknown*&)rtvHeap);
    safeRelease((IUnknown*&)swapChain);
    safeRelease((IUnknown*&)commandList);
    safeRelease((IUnknown*&)allocator);
    safeRelease((IUnknown*&)queue);
    safeRelease((IUnknown*&)fence);
    safeRelease((IUnknown*&)device);
    safeRelease((IUnknown*&)adapter);
    safeRelease((IUnknown*&)factory);
}

void DXRenderer::safeRelease(IUnknown*& p)
{
    if (p)
    {
        p->Release();
        p = nullptr;
    }
}

void DXRenderer::setRenderMode(RenderMode mode){
    renderMode = mode;
}

bool DXRenderer::loadShader(const wchar_t* filename, ID3DBlob** blob)
{
    HANDLE h = CreateFileW(filename, GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE)
    {
        char buf[512];
        sprintf_s(buf, "loadShader: CreateFileW failed for %ls\n", filename);
        OutputDebugStringA(buf);
        return false;
    }

    LARGE_INTEGER size;
    if (!GetFileSizeEx(h, &size))
    {
        OutputDebugStringA("loadShader: GetFileSizeEx failed\n");
        CloseHandle(h);
        return false;
    }

    SIZE_T sz = static_cast<SIZE_T>(size.QuadPart);
    ID3DBlob* b = nullptr;
    HRESULT hr = D3DCreateBlob(sz, &b);
    if (FAILED(hr))
    {
        PrintHR("loadShader: D3DCreateBlob failed", hr);
        CloseHandle(h);
        return false;
    }

    DWORD read = 0;
    BOOL ok = ReadFile(h, b->GetBufferPointer(), (DWORD)sz, &read, nullptr);
    CloseHandle(h);

    if (!ok || read != (DWORD)sz)
    {
        OutputDebugStringA("loadShader: ReadFile failed or incomplete\n");
        b->Release();
        return false;
    }

    *blob = b;
    return true;
}
bool DXRenderer::createFactoryAndDevice()
{
    #if defined(_DEBUG)
    // Enable the debug layer
    ID3D12Debug* debugController = nullptr;
    if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debugController))))
    {
        debugController->EnableDebugLayer();
        debugController->Release();
    }
    #endif
    HRESULT hr = CreateDXGIFactory2(0, IID_PPV_ARGS(&factory));
    if (FAILED(hr))
        return false;

    IDXGIFactory6* factory6 = nullptr;
    hr = factory->QueryInterface(IID_PPV_ARGS(&factory6));
    if (FAILED(hr))
        return false;

    IDXGIAdapter1* adapter = nullptr;
    hr = factory6->EnumAdapterByGpuPreference(
        0,
        DXGI_GPU_PREFERENCE_HIGH_PERFORMANCE,
        IID_PPV_ARGS(&adapter));

    if (FAILED(hr))
    {
        factory6->Release();
        return false;
    }

    DXGI_ADAPTER_DESC1 desc;
    adapter->GetDesc1(&desc);

    char name[256];
    wcstombs_s(nullptr, name, desc.Description, 256);
    OutputDebugStringA("Using adapter: ");
    OutputDebugStringA(name);
    OutputDebugStringA("\n");

    hr = D3D12CreateDevice(
        adapter,
        D3D_FEATURE_LEVEL_12_1,
        IID_PPV_ARGS(&device));

    adapter->Release();
    factory6->Release();

    return SUCCEEDED(hr);
}

bool DXRenderer::createCommandObjects()
{
    D3D12_COMMAND_QUEUE_DESC q = {};
    HRESULT hr = device->CreateCommandQueue(&q, IID_PPV_ARGS(&queue));
    if (FAILED(hr)) { PrintHR("CreateCommandQueue FAILED", hr); return false; }

    hr = device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&allocator));
    if (FAILED(hr)) { PrintHR("CreateCommandAllocator FAILED", hr); return false; }

    hr = device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, allocator, nullptr, IID_PPV_ARGS(&commandList));
    if (FAILED(hr)) { PrintHR("CreateCommandList FAILED", hr); return false; }

    hr = commandList->Close();
    if (FAILED(hr)) { PrintHR("commandList->Close() FAILED", hr); return false; }

    hr = device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence));
    if (FAILED(hr)) { PrintHR("CreateFence FAILED", hr); return false; }

    fenceValue = 1;
    fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!fenceEvent) { OutputDebugStringA("CreateEvent FAILED\n"); return false; }

    return true;
}

bool DXRenderer::createSwapChainAndRTVs(HWND hwnd)
{
    if (!factory || !device || !queue)
    {
        OutputDebugStringA("createSwapChainAndRTVs: factory/device/queue missing\n");
        return false;
    }

    if (swapChain)
    {
        waitForGpu();

        for (int i = 0; i < 2; ++i)
            safeRelease((IUnknown*&)renderTargets[i]);

        safeRelease((IUnknown*&)rtvHeap);

        safeRelease((IUnknown*&)swapChain);
    }

    DXGI_SWAP_CHAIN_DESC1 sc = {};
    sc.BufferCount = 2;
    sc.Width = width;
    sc.Height = height;
    sc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
    sc.SampleDesc.Count = 1;
    sc.SampleDesc.Quality = 0;
    sc.AlphaMode = DXGI_ALPHA_MODE_IGNORE;
    sc.Scaling = DXGI_SCALING_STRETCH;

    IDXGISwapChain1* temp = nullptr;
    HRESULT hr = factory->CreateSwapChainForHwnd(queue, hwnd, &sc, nullptr, nullptr, &temp);
    if (FAILED(hr))
    {
        PrintHR("CreateSwapChainForHwnd FAILED", hr);
        return false;
    }

    hr = temp->QueryInterface(IID_PPV_ARGS(&swapChain));
    temp->Release();
    if (FAILED(hr))
    {
        PrintHR("QueryInterface for IDXGISwapChain3 FAILED", hr);
        return false;
    }
    factory->MakeWindowAssociation(hwnd, DXGI_MWA_NO_ALT_ENTER);

    D3D12_DESCRIPTOR_HEAP_DESC hd = {};
    hd.NumDescriptors = 2;
    hd.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
    hd.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    hr = device->CreateDescriptorHeap(&hd, IID_PPV_ARGS(&rtvHeap));
    if (FAILED(hr)) { PrintHR("CreateDescriptorHeap FAILED", hr); return false; }

    rtvDescriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
    auto h = rtvHeap->GetCPUDescriptorHandleForHeapStart();

    for (int i = 0; i < 2; i++)
    {
        hr = swapChain->GetBuffer(i, IID_PPV_ARGS(&renderTargets[i]));
        if (FAILED(hr)) { PrintHR("GetBuffer FAILED", hr); return false; }
        device->CreateRenderTargetView(renderTargets[i], nullptr, h);
        h.ptr = UINT64(h.ptr) + rtvDescriptorSize;
    }

    backBufferIndex = swapChain->GetCurrentBackBufferIndex();
    return true;
}



bool DXRenderer::createPipeline()
{
    D3D12_ROOT_SIGNATURE_DESC rs = {};
    rs.Flags = D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT;

    ID3DBlob* sig = nullptr;
    ID3DBlob* err = nullptr;
    HRESULT hr = D3D12SerializeRootSignature(&rs, D3D_ROOT_SIGNATURE_VERSION_1, &sig, &err);
    if (FAILED(hr))
    {
        PrintHR("D3D12SerializeRootSignature FAILED", hr);
        if (err) { OutputDebugStringA((char*)err->GetBufferPointer()); err->Release(); }
        return false;
    }

    hr = device->CreateRootSignature(0, sig->GetBufferPointer(), sig->GetBufferSize(), IID_PPV_ARGS(&rootSignature));
    if (FAILED(hr)) { PrintHR("CreateRootSignature FAILED", hr); if (sig) sig->Release(); return false; }
    if (sig) sig->Release();
    if (err) { err->Release(); err = nullptr; }

    ID3DBlob* vs = nullptr;
    ID3DBlob* ps = nullptr;
    if (!loadShader(L"TriangleVS.cso", &vs)) { OutputDebugStringA("Failed to load TriangleVS.cso\n"); return false; }
    if (!loadShader(L"TrianglePS.cso", &ps)) { OutputDebugStringA("Failed to load TrianglePS.cso\n"); if (vs) vs->Release(); return false; }

    if (!vs || !ps) { OutputDebugStringA("Shader blob missing\n"); if (vs) vs->Release(); if (ps) ps->Release(); return false; }

    D3D12_INPUT_ELEMENT_DESC layout[] =
        {
            { "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0,
             D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 }
        };

    D3D12_GRAPHICS_PIPELINE_STATE_DESC pso = {};
    pso.InputLayout = { layout, 1 };
    pso.pRootSignature = rootSignature;
    pso.VS = { vs->GetBufferPointer(), vs->GetBufferSize() };
    pso.PS = { ps->GetBufferPointer(), ps->GetBufferSize() };
    pso.RasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
    pso.BlendState = CD3DX12_BLEND_DESC(D3D12_DEFAULT);
    pso.SampleMask = UINT_MAX;
    pso.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
    pso.NumRenderTargets = 1;
    pso.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
    pso.SampleDesc.Count = 1;

    hr = device->CreateGraphicsPipelineState(&pso, IID_PPV_ARGS(&pipelineState));
    if (FAILED(hr)) { PrintHR("CreateGraphicsPipelineState FAILED", hr); vs->Release(); ps->Release(); return false; }

    vs->Release();
    ps->Release();

    return true;
}


bool DXRenderer::createVertexBuffer()
{
    Vertex verts[] =
        {
            {  0.0f,  0.5f },
            {  0.5f, -0.5f },
            { -0.5f, -0.5f }
        };

    auto heap = CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD);
    auto desc = CD3DX12_RESOURCE_DESC::Buffer(sizeof(verts));

    HRESULT hr = device->CreateCommittedResource(
        &heap,
        D3D12_HEAP_FLAG_NONE,
        &desc,
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&vertexBuffer));
    if (FAILED(hr)) { PrintHR("CreateCommittedResource(vertexBuffer) FAILED", hr); return false; }

    void* data;
    D3D12_RANGE r = {};
    hr = vertexBuffer->Map(0, &r, &data);
    if (FAILED(hr)) { PrintHR("vertexBuffer->Map FAILED", hr); return false; }

    memcpy(data, verts, sizeof(verts));

    vertexBuffer->Unmap(0, nullptr);

    vertexBufferView.BufferLocation = vertexBuffer->GetGPUVirtualAddress();
    vertexBufferView.StrideInBytes = sizeof(Vertex);
    vertexBufferView.SizeInBytes = sizeof(verts);

    return true;
}

bool DXRenderer::createReadbackBuffer()
{
    auto desc = renderTargets[0]->GetDesc();

    D3D12_PLACED_SUBRESOURCE_FOOTPRINT fp;
    UINT rows;
    UINT64 bytes;
    device->GetCopyableFootprints(&desc, 0, 1, 0, &fp, &rows, nullptr, &bytes);

    rowPitch = fp.Footprint.RowPitch;

    auto heap = CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_READBACK);
    auto buf = CD3DX12_RESOURCE_DESC::Buffer(bytes);

    HRESULT hr = device->CreateCommittedResource(
        &heap,
        D3D12_HEAP_FLAG_NONE,
        &buf,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&readbackBuffer));
    if (FAILED(hr)) { PrintHR("CreateCommittedResource(readback) FAILED", hr); return false; }

    return true;
}


bool DXRenderer::initialize(HWND hwnd, int w, int h)
{
    width = w;
    height = h;

    if (!createFactoryAndDevice()) return false;
    if (!createCommandObjects()) return false;
    if (!createSwapChainAndRTVs(hwnd)) return false;
    if (!createPipeline()) return false;
    if (!createVertexBuffer()) return false;
    if (!createReadbackBuffer()) return false;

    rayTracing = new DXRendererRayTracing();
    rayTracing->initialize(device, queue, hwnd, width, height);
    setRenderMode(RenderMode::RayTracing);
    return true;
}


void DXRenderer::waitForGpu()
{
    HRESULT hr = queue->Signal(fence, fenceValue);
    if (FAILED(hr)) { PrintHR("queue->Signal FAILED", hr); }

    if (fence->GetCompletedValue() < fenceValue)
    {
        fence->SetEventOnCompletion(fenceValue, fenceEvent);
        WaitForSingleObject(fenceEvent, INFINITE);
    }
    fenceValue++;
}


void DXRenderer::renderRaster()
{
    HRESULT hr = allocator->Reset();
    if (FAILED(hr)) { PrintHR("allocator->Reset FAILED", hr); return; }

    hr = commandList->Reset(allocator, pipelineState);
    if (FAILED(hr)) { PrintHR("commandList->Reset FAILED", hr); return; }

    auto barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex],
        D3D12_RESOURCE_STATE_PRESENT,
        D3D12_RESOURCE_STATE_RENDER_TARGET);
    commandList->ResourceBarrier(1, &barrier);

    auto rtv = rtvHeap->GetCPUDescriptorHandleForHeapStart();
    rtv.ptr += backBufferIndex * rtvDescriptorSize;

    float clear[4] = { 0.1f, 0.1f, 0.1f, 1.0f };
    commandList->ClearRenderTargetView(rtv, clear, 0, nullptr);

    commandList->OMSetRenderTargets(1, &rtv, FALSE, nullptr);

    D3D12_VIEWPORT viewport{};
    viewport.Width = (float)width;
    viewport.Height = (float)height;
    viewport.MaxDepth = 1.0f;

    D3D12_RECT scissor{ 0, 0, width, height };
    commandList->RSSetViewports(1, &viewport);
    commandList->RSSetScissorRects(1, &scissor);

    commandList->SetGraphicsRootSignature(rootSignature);
    commandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    commandList->IASetVertexBuffers(0, 1, &vertexBufferView);
    commandList->DrawInstanced(3, 1, 0, 0);

    barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex],
        D3D12_RESOURCE_STATE_RENDER_TARGET,
        D3D12_RESOURCE_STATE_COPY_SOURCE);
    commandList->ResourceBarrier(1, &barrier);

    D3D12_TEXTURE_COPY_LOCATION src{};
    src.pResource = renderTargets[backBufferIndex];
    src.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
    src.SubresourceIndex = 0;

    D3D12_TEXTURE_COPY_LOCATION dst{};
    dst.pResource = readbackBuffer;
    dst.Type = D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;

    D3D12_RESOURCE_DESC desc = renderTargets[backBufferIndex]->GetDesc();
    D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint{};
    UINT rows;
    UINT64 bytes;
    device->GetCopyableFootprints(
        &desc, 0, 1, 0,
        &footprint, &rows, nullptr, &bytes);

    dst.PlacedFootprint = footprint;

    commandList->CopyTextureRegion(&dst, 0, 0, 0, &src, nullptr);

    barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex],
        D3D12_RESOURCE_STATE_COPY_SOURCE,
        D3D12_RESOURCE_STATE_PRESENT);
    commandList->ResourceBarrier(1, &barrier);

    hr = commandList->Close();
    if (FAILED(hr)) { PrintHR("commandList->Close FAILED", hr); return; }

    ID3D12CommandList* lists[] = { commandList };
    queue->ExecuteCommandLists(1, lists);
    swapChain->Present(1, 0);

    waitForGpu();
    backBufferIndex = swapChain->GetCurrentBackBufferIndex();
}

void DXRenderer::renderRayTracing()
{
    allocator->Reset();
    commandList->Reset(allocator, nullptr);

    ID3D12GraphicsCommandList4* cmd4 = nullptr;
    commandList->QueryInterface(IID_PPV_ARGS(&cmd4));

    rayTracing->renderFrame(cmd4);

    auto barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        rayTracing->getOutputTexture(),
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
        D3D12_RESOURCE_STATE_COPY_SOURCE);
    commandList->ResourceBarrier(1, &barrier);

    barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex],
        D3D12_RESOURCE_STATE_PRESENT,
        D3D12_RESOURCE_STATE_COPY_DEST);
    commandList->ResourceBarrier(1, &barrier);


    D3D12_TEXTURE_COPY_LOCATION src{};
    src.pResource = rayTracing->getOutputTexture();
    src.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
    src.SubresourceIndex = 0;

    D3D12_TEXTURE_COPY_LOCATION dst{};
    dst.pResource = renderTargets[backBufferIndex];
    dst.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
    dst.SubresourceIndex = 0;

    commandList->CopyTextureRegion(&dst, 0, 0, 0, &src, nullptr);

    barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex],
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_PRESENT);
    commandList->ResourceBarrier(1, &barrier);

    barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        rayTracing->getOutputTexture(),
        D3D12_RESOURCE_STATE_COPY_SOURCE,
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
    commandList->ResourceBarrier(1, &barrier);

    commandList->Close();

    ID3D12CommandList* lists[] = { commandList };
    queue->ExecuteCommandLists(1, lists);
    swapChain->Present(1, 0);

    waitForGpu();
    backBufferIndex = swapChain->GetCurrentBackBufferIndex();

    cmd4->Release();
}

void DXRenderer::renderFrame()
{
    if (renderMode == RenderMode::RayTracing)
        renderRayTracing();
    else
        renderRaster();
}

unsigned char* DXRenderer::getBackBufferCPU()
{
    if (!mappedPtr)
    {
        D3D12_RANGE r = { 0, 0 };
        HRESULT hr = readbackBuffer->Map(0, &r, reinterpret_cast<void**>(&mappedPtr));
        if (FAILED(hr)) { PrintHR("readbackBuffer->Map FAILED", hr); return nullptr; }
    }
    return mappedPtr;
}

void DXRenderer::unmapBackBuffer()
{
    if (mappedPtr)
    {
        readbackBuffer->Unmap(0, nullptr);
        mappedPtr = nullptr;
    }
}
