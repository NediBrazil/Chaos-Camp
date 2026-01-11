#include "DXRendererRayTracing.h"
#include <d3dcompiler.h>
#include <vector>

bool DXRendererRayTracing::initialize(ID3D12Device* baseDevice, ID3D12CommandQueue* baseQueue, int w, int h)
{

    device = baseDevice;
    queue = baseQueue;
    width = w;
    height = h;

    if (FAILED(device->QueryInterface(IID_PPV_ARGS(&device5))))
        return false;

    if (!createOutputTexture()) return false;
    if (!createDescriptorHeap()) return false;
    if (!createGlobalRootSignature()) return false;
    if (!createStateObject()) return false;
    if (!createShaderBindingTable()) return false;

    return true;
}

bool DXRendererRayTracing::createOutputTexture()
{
    D3D12_RESOURCE_DESC d = {};
    d.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
    d.Width = width;
    d.Height = height;
    d.DepthOrArraySize = 1;
    d.MipLevels = 1;
    d.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    d.SampleDesc.Count = 1;
    d.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
    d.Flags = D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;

    D3D12_HEAP_PROPERTIES heap = {};
    heap.Type = D3D12_HEAP_TYPE_DEFAULT;

    HRESULT hr = device->CreateCommittedResource(
        &heap,
        D3D12_HEAP_FLAG_NONE,
        &d,
        D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
        nullptr,
        IID_PPV_ARGS(&outputTexture)
        );

    return SUCCEEDED(hr);
}

bool DXRendererRayTracing::createDescriptorHeap()
{
    D3D12_DESCRIPTOR_HEAP_DESC h = {};
    h.NumDescriptors = 1;
    h.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    h.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;

    if (FAILED(device->CreateDescriptorHeap(&h, IID_PPV_ARGS(&descriptorHeap))))
        return false;

    D3D12_UNORDERED_ACCESS_VIEW_DESC u = {};
    u.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    u.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE2D;

    device->CreateUnorderedAccessView(
        outputTexture,
        nullptr,
        &u,
        descriptorHeap->GetCPUDescriptorHandleForHeapStart()
        );

    return true;
}

bool DXRendererRayTracing::createGlobalRootSignature()
{
    CD3DX12_DESCRIPTOR_RANGE1 range;
    range.Init(D3D12_DESCRIPTOR_RANGE_TYPE_UAV, 1, 0);

    CD3DX12_ROOT_PARAMETER1 param;
    param.InitAsDescriptorTable(1, &range);

    CD3DX12_VERSIONED_ROOT_SIGNATURE_DESC desc;
    desc.Init_1_1(1, &param, 0, nullptr);

    ID3DBlob* blob = nullptr;
    if (FAILED(D3D12SerializeVersionedRootSignature(&desc, &blob, nullptr)))
        return false;

    bool ok = SUCCEEDED(device->CreateRootSignature(
        0,
        blob->GetBufferPointer(),
        blob->GetBufferSize(),
        IID_PPV_ARGS(&globalRootSignature)
        ));

    blob->Release();
    return ok;
}

bool DXRendererRayTracing::createStateObject()
{
    D3D12_FEATURE_DATA_D3D12_OPTIONS5 opt = {};
    if (FAILED(device->CheckFeatureSupport(
            D3D12_FEATURE_D3D12_OPTIONS5,
            &opt,
            sizeof(opt))))
        return false;

    if (opt.RaytracingTier == D3D12_RAYTRACING_TIER_NOT_SUPPORTED)
    {
        MessageBoxA(0, "DXR NOT SUPPORTED ON THIS GPU", "DXR", MB_OK);
        return false;
    }

    ID3DBlob* raygenBlob = nullptr;
    ID3DBlob* missBlob = nullptr;

    if (FAILED(D3DReadFileToBlob(L"RayGen.cso", &raygenBlob))) return false;
    if (FAILED(D3DReadFileToBlob(L"Miss.cso", &missBlob))) return false;

    std::vector<D3D12_EXPORT_DESC> exports;
    exports.reserve(2);
    exports.push_back({ L"RayGen", nullptr, D3D12_EXPORT_FLAG_NONE });
    exports.push_back({ L"Miss", nullptr, D3D12_EXPORT_FLAG_NONE });

    std::vector<D3D12_DXIL_LIBRARY_DESC> libs;
    libs.resize(2);
    libs[0].DXILLibrary.pShaderBytecode = raygenBlob->GetBufferPointer();
    libs[0].DXILLibrary.BytecodeLength = raygenBlob->GetBufferSize();
    libs[0].NumExports = 1;
    libs[0].pExports = &exports[0];
    libs[1].DXILLibrary.pShaderBytecode = missBlob->GetBufferPointer();
    libs[1].DXILLibrary.BytecodeLength = missBlob->GetBufferSize();
    libs[1].NumExports = 1;
    libs[1].pExports = &exports[1];

    D3D12_RAYTRACING_SHADER_CONFIG shaderConfig = {};
    shaderConfig.MaxPayloadSizeInBytes = sizeof(float) * 4;
    shaderConfig.MaxAttributeSizeInBytes = 0;

    D3D12_RAYTRACING_PIPELINE_CONFIG pipelineConfig = {};
    pipelineConfig.MaxTraceRecursionDepth = 1;

    LPCWSTR exportNames[] = { L"RayGen", L"Miss" };

    std::vector<D3D12_STATE_SUBOBJECT> subs;
    subs.reserve(8);

    D3D12_STATE_SUBOBJECT raygenLibSub = {};
    raygenLibSub.Type = D3D12_STATE_SUBOBJECT_TYPE_DXIL_LIBRARY;
    raygenLibSub.pDesc = &libs[0];
    subs.push_back(raygenLibSub);

    D3D12_STATE_SUBOBJECT missLibSub = {};
    missLibSub.Type = D3D12_STATE_SUBOBJECT_TYPE_DXIL_LIBRARY;
    missLibSub.pDesc = &libs[1];
    subs.push_back(missLibSub);

    D3D12_STATE_SUBOBJECT shaderConfigSub = {};
    shaderConfigSub.Type = D3D12_STATE_SUBOBJECT_TYPE_RAYTRACING_SHADER_CONFIG;
    shaderConfigSub.pDesc = &shaderConfig;
    subs.push_back(shaderConfigSub);

    D3D12_SUBOBJECT_TO_EXPORTS_ASSOCIATION shaderConfigAssoc = {};
    shaderConfigAssoc.NumExports = 2;
    shaderConfigAssoc.pExports = exportNames;
    shaderConfigAssoc.pSubobjectToAssociate = &subs[2];
    D3D12_STATE_SUBOBJECT shaderConfigAssocSub = {};
    shaderConfigAssocSub.Type = D3D12_STATE_SUBOBJECT_TYPE_SUBOBJECT_TO_EXPORTS_ASSOCIATION;
    shaderConfigAssocSub.pDesc = &shaderConfigAssoc;
    subs.push_back(shaderConfigAssocSub);

    D3D12_STATE_SUBOBJECT pipelineConfigSub = {};
    pipelineConfigSub.Type = D3D12_STATE_SUBOBJECT_TYPE_RAYTRACING_PIPELINE_CONFIG;
    pipelineConfigSub.pDesc = &pipelineConfig;
    subs.push_back(pipelineConfigSub);

    D3D12_STATE_SUBOBJECT rootSigSub = {};
    rootSigSub.Type = D3D12_STATE_SUBOBJECT_TYPE_GLOBAL_ROOT_SIGNATURE;
    rootSigSub.pDesc = &globalRootSignature;
    subs.push_back(rootSigSub);

    D3D12_STATE_OBJECT_DESC desc = {};
    desc.Type = D3D12_STATE_OBJECT_TYPE_RAYTRACING_PIPELINE;
    desc.NumSubobjects = (UINT)subs.size();
    desc.pSubobjects = subs.data();

    HRESULT hr = device5->CreateStateObject(&desc, IID_PPV_ARGS(&stateObject));
    if (FAILED(hr))
    {
        raygenBlob->Release();
        missBlob->Release();
        return false;
    }

    stateObject->QueryInterface(IID_PPV_ARGS(&stateProps));
    raygenBlob->Release();
    missBlob->Release();
    return true;
}

bool DXRendererRayTracing::createShaderBindingTable()
{
    UINT idSize = D3D12_SHADER_IDENTIFIER_SIZE_IN_BYTES;
    UINT align = D3D12_RAYTRACING_SHADER_RECORD_BYTE_ALIGNMENT;

    UINT recordSize = (idSize + align - 1) & ~(align - 1);

    UINT tableAlign = D3D12_RAYTRACING_SHADER_TABLE_BYTE_ALIGNMENT;
    UINT raygenSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);
    UINT missSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);

    UINT sbtSize = raygenSize + missSize;

    D3D12_HEAP_PROPERTIES heap = {};
    heap.Type = D3D12_HEAP_TYPE_UPLOAD;

    D3D12_RESOURCE_DESC desc = {};
    desc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    desc.Width = sbtSize;
    desc.Height = 1;
    desc.DepthOrArraySize = 1;
    desc.MipLevels = 1;
    desc.SampleDesc.Count = 1;
    desc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    if (FAILED(device->CreateCommittedResource(
            &heap,
            D3D12_HEAP_FLAG_NONE,
            &desc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&sbtBuffer))))
        return false;

    uint8_t* p = nullptr;
    sbtBuffer->Map(0, nullptr, reinterpret_cast<void**>(&p));
    memcpy(p, stateProps->GetShaderIdentifier(L"RayGen"), idSize);
    memcpy(p + raygenSize, stateProps->GetShaderIdentifier(L"Miss"), idSize);
    sbtBuffer->Unmap(0, nullptr);

    return true;
}

void DXRendererRayTracing::renderFrame(ID3D12GraphicsCommandList4* cmd)
{
    ID3D12DescriptorHeap* heaps[] = { descriptorHeap };
    cmd->SetDescriptorHeaps(1, heaps);
    cmd->SetComputeRootSignature(globalRootSignature);
    cmd->SetComputeRootDescriptorTable(0, descriptorHeap->GetGPUDescriptorHandleForHeapStart());
    cmd->SetPipelineState1(stateObject);

    UINT idSize = D3D12_SHADER_IDENTIFIER_SIZE_IN_BYTES;
    UINT align = D3D12_RAYTRACING_SHADER_RECORD_BYTE_ALIGNMENT;
    UINT recordSize = (idSize + align - 1) & ~(align - 1);

    UINT tableAlign = D3D12_RAYTRACING_SHADER_TABLE_BYTE_ALIGNMENT;
    UINT raygenSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);
    UINT missSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);

    D3D12_GPU_VIRTUAL_ADDRESS addr = sbtBuffer->GetGPUVirtualAddress();

    D3D12_DISPATCH_RAYS_DESC d = {};
    d.RayGenerationShaderRecord.StartAddress = addr;
    d.RayGenerationShaderRecord.SizeInBytes = recordSize;

    d.MissShaderTable.StartAddress = addr + raygenSize;
    d.MissShaderTable.SizeInBytes = missSize;
    d.MissShaderTable.StrideInBytes = recordSize;

    d.Width = width;
    d.Height = height;
    d.Depth = 1;

    cmd->DispatchRays(&d);

    D3D12_RESOURCE_BARRIER b = {};
    b.Type = D3D12_RESOURCE_BARRIER_TYPE_UAV;
    b.UAV.pResource = outputTexture;
    cmd->ResourceBarrier(1, &b);
}

ID3D12Resource* DXRendererRayTracing::getOutputTexture()
{
    return outputTexture;
}

DXRendererRayTracing::~DXRendererRayTracing()
{
    if (stateProps) {
        stateProps->Release();
        stateProps = nullptr;
    }
    if (stateObject) {
        stateObject->Release();
        stateObject = nullptr;
    }
    if (globalRootSignature) {
        globalRootSignature->Release();
        globalRootSignature = nullptr;
    }
    if (descriptorHeap) {
        descriptorHeap->Release();
        descriptorHeap = nullptr;
    }
    if (outputTexture) {
        outputTexture->Release();
        outputTexture = nullptr;
    }
    if (device5) {
        device5->Release();
        device5 = nullptr;
    }
    if (sbtBuffer) {
        sbtBuffer->Release();
        sbtBuffer = nullptr;
    }
}
