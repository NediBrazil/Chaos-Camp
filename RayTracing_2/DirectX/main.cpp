#include <iostream>
#include <dxgi1_4.h>
#include <d3d12.h>
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d12.lib")

class CDXCRenderer
{
public:
    void prepareForRendering();
    void render();

private:
    void createDevice();
    void createCommandsManagers();
    void createRenderTarget();

    IDXGIFactory4 *dxgiFactory = nullptr;
    IDXGIAdapter1 *adapter = nullptr;
    ID3D12Device *d3d12Device = nullptr;

    ID3D12CommandQueue *commandQueue = nullptr;
    ID3D12CommandAllocator *commandAllocator = nullptr;
    ID3D12GraphicsCommandList *commandList = nullptr;

    ID3D12Resource *renderTarget = nullptr;
    ID3D12DescriptorHeap *rtvHeap = nullptr;
    D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = {};
};

void CDXCRenderer::prepareForRendering()
{
    if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&dxgiFactory))))
    {
        std::cout << "Failed to create DXGIFactory\n";
        return;
    }

    for (UINT i = 0; dxgiFactory->EnumAdapters1(i, &adapter) != DXGI_ERROR_NOT_FOUND; ++i)
    {
        DXGI_ADAPTER_DESC1 desc;
        adapter->GetDesc1(&desc);
        if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)
            continue;

        if (SUCCEEDED(D3D12CreateDevice(adapter, D3D_FEATURE_LEVEL_11_0, _uuidof(ID3D12Device), nullptr)))
        {
            std::wcout << L"Using adapter: " << desc.Description << std::endl;
            break;
        }
        adapter->Release();
    }

    if (!adapter)
    {
        std::cout << "No suitable GPU found\n";
        return;
    }

    createDevice();
    createCommandsManagers();
    createRenderTarget();
}

void CDXCRenderer::createDevice()
{
    if (FAILED(D3D12CreateDevice(adapter, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&d3d12Device))))
    {
        std::cout << "Failed to create D3D12 device\n";
        return;
    }
}

void CDXCRenderer::createCommandsManagers()
{
    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;

    if (FAILED(d3d12Device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&commandQueue))))
    {
        std::cout << "Failed to create Command Queue\n";
        return;
    }

    if (FAILED(d3d12Device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&commandAllocator))))
    {
        std::cout << "Failed to create Command Allocator\n";
        return;
    }

    if (FAILED(d3d12Device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, commandAllocator, nullptr, IID_PPV_ARGS(&commandList))))
    {
        std::cout << "Failed to create Command List\n";
        return;
    }
}

void CDXCRenderer::createRenderTarget()
{
    D3D12_HEAP_PROPERTIES heapProps = {};
    heapProps.Type = D3D12_HEAP_TYPE_DEFAULT;

    D3D12_RESOURCE_DESC texDesc = {};
    texDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
    texDesc.Width = 800;
    texDesc.Height = 600;
    texDesc.DepthOrArraySize = 1;
    texDesc.MipLevels = 1;
    texDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    texDesc.SampleDesc.Count = 1;
    texDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
    texDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET;

    D3D12_CLEAR_VALUE clearValue = {};
    clearValue.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    clearValue.Color[0] = 0.1f;
    clearValue.Color[1] = 0.2f;
    clearValue.Color[2] = 0.4f;
    clearValue.Color[3] = 1.0f;

    if (FAILED(d3d12Device->CreateCommittedResource(
            &heapProps,
            D3D12_HEAP_FLAG_NONE,
            &texDesc,
            D3D12_RESOURCE_STATE_RENDER_TARGET,
            &clearValue,
            IID_PPV_ARGS(&renderTarget))))
    {
        std::cout << "Failed to create render target resource\n";
        return;
    }

    D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
    heapDesc.NumDescriptors = 1;
    heapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;

    if (FAILED(d3d12Device->CreateDescriptorHeap(&heapDesc, IID_PPV_ARGS(&rtvHeap))))
    {
        std::cout << "Failed to create RTV heap\n";
        return;
    }

    rtvHandle = rtvHeap->GetCPUDescriptorHandleForHeapStart();
    d3d12Device->CreateRenderTargetView(renderTarget, nullptr, rtvHandle);
}

void CDXCRenderer::render()
{
    FLOAT clearColor[] = {0.3f, 0.5f, 0.8f, 1.0f};

    commandList->OMSetRenderTargets(1, &rtvHandle, FALSE, nullptr);
    commandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);

    std::cout << "Render target cleared successfully\n";
}

int main()
{
    CDXCRenderer renderer;
    renderer.prepareForRendering();
    renderer.render();
    return 0;
}
