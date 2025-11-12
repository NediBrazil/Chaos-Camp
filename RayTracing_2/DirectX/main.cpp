#include <iostream>
#include <dxgi.h>
#pragma comment(lib, "dxgi.lib")

int main()
{
    IDXGIFactory *pFactory = nullptr;
    if (FAILED(CreateDXGIFactory(__uuidof(IDXGIFactory), (void **)&pFactory)))
    {
        std::cout << "Failed to create DXGIFactory\n";
        return -1;
    }

    IDXGIAdapter *pAdapter = nullptr;
    UINT index = 0;
    while (pFactory->EnumAdapters(index, &pAdapter) != DXGI_ERROR_NOT_FOUND)
    {
        DXGI_ADAPTER_DESC desc;
        pAdapter->GetDesc(&desc);

        std::wcout << L"Adapter " << index << L": " << desc.Description << std::endl;

        switch (desc.VendorId)
        {
        case 0x10DE:
            std::wcout << L"Vendor: NVIDIA" << std::endl;
            break;
        case 0x1002:
            std::wcout << L"Vendor: AMD" << std::endl;
            break;
        case 0x8086:
            std::wcout << L"Vendor: Intel" << std::endl;
            break;
        default:
            std::wcout << L"Vendor: Unknown (0x" << std::hex << desc.VendorId << L")" << std::endl;
            break;
        }

        std::wcout << L"Dedicated Video Memory: " << (desc.DedicatedVideoMemory / (1024 * 1024)) << L" MB" << std::endl;
        std::wcout << std::endl;

        pAdapter->Release();
        index++;
    }

    pFactory->Release();
    return 0;
}
