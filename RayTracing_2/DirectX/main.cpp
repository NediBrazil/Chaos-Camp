#include <iostream>
#include <dxgi1_4.h>
#include <d3d12.h>
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d12.lib")

int main()
{
    IDXGIFactory4 *pFactory = nullptr;
    if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&pFactory))))
    {
        std::cout << "Failed to create DXGIFactory\n";
        return -1;
    }

    IDXGIAdapter1 *pAdapter = nullptr;
    for (UINT i = 0; DXGI_ERROR_NOT_FOUND != pFactory->EnumAdapters1(i, &pAdapter); ++i)
    {
        DXGI_ADAPTER_DESC1 desc;
        pAdapter->GetDesc1(&desc);

        if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)
            continue;

        if (SUCCEEDED(D3D12CreateDevice(pAdapter, D3D_FEATURE_LEVEL_11_0, _uuidof(ID3D12Device), nullptr)))
        {
            std::wcout << L"Using adapter: " << desc.Description << std::endl;
            break;
        }
        pAdapter->Release();
    }

    if (!pAdapter)
    {
        std::cout << "No suitable GPU found\n";
        pFactory->Release();
        return -1;
    }

    ID3D12Device *pDevice = nullptr;
    if (FAILED(D3D12CreateDevice(pAdapter, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&pDevice))))
    {
        std::cout << "Failed to create D3D12 device\n";
        pAdapter->Release();
        pFactory->Release();
        return -1;
    }

    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    ID3D12CommandQueue *pCommandQueue = nullptr;
    if (FAILED(pDevice->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&pCommandQueue))))
    {
        std::cout << "Failed to create Command Queue\n";
        pDevice->Release();
        pAdapter->Release();
        pFactory->Release();
        return -1;
    }

    ID3D12CommandAllocator *pCommandAllocator = nullptr;
    if (FAILED(pDevice->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&pCommandAllocator))))
    {
        std::cout << "Failed to create Command Allocator\n";
        pCommandQueue->Release();
        pDevice->Release();
        pAdapter->Release();
        pFactory->Release();
        return -1;
    }

    ID3D12GraphicsCommandList *pCommandList = nullptr;
    if (FAILED(pDevice->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, pCommandAllocator, nullptr, IID_PPV_ARGS(&pCommandList))))
    {
        std::cout << "Failed to create Command List\n";
        pCommandAllocator->Release();
        pCommandQueue->Release();
        pDevice->Release();
        pAdapter->Release();
        pFactory->Release();
        return -1;
    }

    std::cout << "Direct3D 12 initialization successful\n";

    pCommandList->Release();
    pCommandAllocator->Release();
    pCommandQueue->Release();
    pDevice->Release();
    pAdapter->Release();
    pFactory->Release();
    return 0;
}
