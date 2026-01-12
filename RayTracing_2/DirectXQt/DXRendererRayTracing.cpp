#include "DXRendererRayTracing.h"
#include <d3dcompiler.h>
#include <vector>
#include <DirectXMath.h>

struct CameraCB
{
    DirectX::XMFLOAT3 pos;
    float pad0;
    DirectX::XMFLOAT3 forward;
    float pad1;
    DirectX::XMFLOAT3 right;
    float pad2;
    DirectX::XMFLOAT3 up;
    float fov;
};
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
    if (!createCameraCB()) return false;
    createCameraCBV();
    if (!createTriangleBuffers()) return false;
    if (!createGlobalRootSignature()) return false;
    if (!createStateObject()) return false;
    if (!createShaderBindingTable()) return false;
    return true;
}
bool DXRendererRayTracing::createCameraCB()
{
    D3D12_HEAP_PROPERTIES heap = {};
    heap.Type = D3D12_HEAP_TYPE_UPLOAD;

    D3D12_RESOURCE_DESC desc = {};
    desc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    desc.Width = (sizeof(CameraCB) + 255) & ~255;
    desc.Height = 1;
    desc.DepthOrArraySize = 1;
    desc.MipLevels = 1;
    desc.SampleDesc.Count = 1;
    desc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    return SUCCEEDED(device->CreateCommittedResource(
        &heap,
        D3D12_HEAP_FLAG_NONE,
        &desc,
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&cameraCB)));
}

void DXRendererRayTracing::createCameraCBV()
{
    D3D12_CONSTANT_BUFFER_VIEW_DESC cbv = {};
    cbv.BufferLocation = cameraCB->GetGPUVirtualAddress();
    cbv.SizeInBytes = (sizeof(CameraCB) + 255) & ~255;

    UINT inc = device->GetDescriptorHandleIncrementSize(
        D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);

    D3D12_CPU_DESCRIPTOR_HANDLE cpu =
        descriptorHeap->GetCPUDescriptorHandleForHeapStart();

    cpu.ptr += inc * 2;

    device->CreateConstantBufferView(&cbv, cpu);
}

void DXRendererRayTracing::updateCamera()
{
    CameraCB cam = {};
    cam.pos     = { 0.0f, 0.0f, -3.0f };
    cam.forward = { 0.0f, 0.0f,  1.0f };
    cam.right   = { 1.0f, 0.0f,  0.0f };
    cam.up      = { 0.0f, 1.0f,  0.0f };
    cam.fov     = 60.0f;

    void* p = nullptr;
    cameraCB->Map(0, nullptr, &p);
    memcpy(p, &cam, sizeof(CameraCB));
    cameraCB->Unmap(0, nullptr);
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
    h.NumDescriptors = 3;
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
    CD3DX12_DESCRIPTOR_RANGE1 ranges[3];
    ranges[0].Init(D3D12_DESCRIPTOR_RANGE_TYPE_UAV, 1, 0);
    ranges[1].Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 0);
    ranges[2].Init(D3D12_DESCRIPTOR_RANGE_TYPE_CBV, 1, 0);

    CD3DX12_ROOT_PARAMETER1 params[3];
    params[0].InitAsDescriptorTable(1, &ranges[0]);
    params[1].InitAsDescriptorTable(1, &ranges[1]);
    params[2].InitAsDescriptorTable(1, &ranges[2]);

    CD3DX12_VERSIONED_ROOT_SIGNATURE_DESC desc;
    desc.Init_1_1(3, params, 0, nullptr);

    ID3DBlob* blob = nullptr;
    if (FAILED(D3D12SerializeVersionedRootSignature(&desc, &blob, nullptr)))
        return false;

    bool ok = SUCCEEDED(device->CreateRootSignature(
        0,
        blob->GetBufferPointer(),
        blob->GetBufferSize(),
        IID_PPV_ARGS(&globalRootSignature)));

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
    ID3DBlob* chBlob = nullptr;

    if (FAILED(D3DReadFileToBlob(L"RayGen.cso", &raygenBlob))) return false;
    if (FAILED(D3DReadFileToBlob(L"Miss.cso", &missBlob))) return false;
    if (FAILED(D3DReadFileToBlob(L"ClosestHit.cso", &chBlob))) return false;

    std::vector<D3D12_EXPORT_DESC> exports;
    exports.reserve(3);
    exports.push_back({ L"RayGen", nullptr, D3D12_EXPORT_FLAG_NONE });
    exports.push_back({ L"Miss", nullptr, D3D12_EXPORT_FLAG_NONE });
    exports.push_back({ L"ClosestHit", nullptr, D3D12_EXPORT_FLAG_NONE });

    std::vector<D3D12_DXIL_LIBRARY_DESC> libs;
    libs.resize(3);
    libs[0].DXILLibrary.pShaderBytecode = raygenBlob->GetBufferPointer();
    libs[0].DXILLibrary.BytecodeLength = raygenBlob->GetBufferSize();
    libs[0].NumExports = 1;
    libs[0].pExports = &exports[0];
    libs[1].DXILLibrary.pShaderBytecode = missBlob->GetBufferPointer();
    libs[1].DXILLibrary.BytecodeLength = missBlob->GetBufferSize();
    libs[1].NumExports = 1;
    libs[1].pExports = &exports[1];
    libs[2].DXILLibrary.pShaderBytecode = chBlob->GetBufferPointer();
    libs[2].DXILLibrary.BytecodeLength = chBlob->GetBufferSize();
    libs[2].NumExports = 1;
    libs[2].pExports = &exports[2];

    D3D12_RAYTRACING_SHADER_CONFIG shaderConfig = {};
    shaderConfig.MaxPayloadSizeInBytes = sizeof(float) * 4;
    shaderConfig.MaxAttributeSizeInBytes = 8;

    D3D12_RAYTRACING_PIPELINE_CONFIG pipelineConfig = {};
    pipelineConfig.MaxTraceRecursionDepth = 1;

    LPCWSTR exportNames[] = { L"RayGen", L"Miss",L"HitGroup" };

    std::vector<D3D12_STATE_SUBOBJECT> subs;
    subs.reserve(10);

    D3D12_STATE_SUBOBJECT raygenLibSub = {};
    raygenLibSub.Type = D3D12_STATE_SUBOBJECT_TYPE_DXIL_LIBRARY;
    raygenLibSub.pDesc = &libs[0];
    subs.push_back(raygenLibSub);

    D3D12_STATE_SUBOBJECT missLibSub = {};
    missLibSub.Type = D3D12_STATE_SUBOBJECT_TYPE_DXIL_LIBRARY;
    missLibSub.pDesc = &libs[1];
    subs.push_back(missLibSub);

    D3D12_STATE_SUBOBJECT chSub = {};
    chSub.Type = D3D12_STATE_SUBOBJECT_TYPE_DXIL_LIBRARY;
    chSub.pDesc = &libs[2];
    subs.push_back(chSub);

    D3D12_HIT_GROUP_DESC hitGroup = {};
    hitGroup.HitGroupExport = L"HitGroup";
    hitGroup.ClosestHitShaderImport = L"ClosestHit";
    hitGroup.Type = D3D12_HIT_GROUP_TYPE_TRIANGLES;

    D3D12_STATE_SUBOBJECT hitGroupSub = {};
    hitGroupSub.Type = D3D12_STATE_SUBOBJECT_TYPE_HIT_GROUP;
    hitGroupSub.pDesc = &hitGroup;
    subs.push_back(hitGroupSub);

    D3D12_STATE_SUBOBJECT shaderConfigSub = {};
    shaderConfigSub.Type = D3D12_STATE_SUBOBJECT_TYPE_RAYTRACING_SHADER_CONFIG;
    shaderConfigSub.pDesc = &shaderConfig;
    subs.push_back(shaderConfigSub);

    D3D12_SUBOBJECT_TO_EXPORTS_ASSOCIATION shaderConfigAssoc = {};
    shaderConfigAssoc.NumExports = 3;
    shaderConfigAssoc.pExports = exportNames;
    shaderConfigAssoc.pSubobjectToAssociate = &subs.back();

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
        chBlob->Release();
        return false;
    }

    stateObject->QueryInterface(IID_PPV_ARGS(&stateProps));
    raygenBlob->Release();
    missBlob->Release();
    chBlob->Release();
    return true;
}

bool DXRendererRayTracing::createShaderBindingTable()
{
    const UINT shaderIdSize = D3D12_SHADER_IDENTIFIER_SIZE_IN_BYTES;
    const UINT shaderRecordAlignment = D3D12_RAYTRACING_SHADER_RECORD_BYTE_ALIGNMENT;
    const UINT shaderTableAlignment  = D3D12_RAYTRACING_SHADER_TABLE_BYTE_ALIGNMENT;

    UINT recordSize = (shaderIdSize + shaderRecordAlignment - 1)
      & ~(shaderRecordAlignment - 1);

    UINT raygenTableSize = (recordSize + shaderTableAlignment - 1)
      & ~(shaderTableAlignment - 1);
    UINT missTableSize   = (recordSize + shaderTableAlignment - 1)
     & ~(shaderTableAlignment - 1);
    UINT hitTableSize    = (recordSize + shaderTableAlignment - 1)
       & ~(shaderTableAlignment - 1);

    UINT sbtSize = raygenTableSize + missTableSize + hitTableSize;

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

    uint8_t* pData = nullptr;
    sbtBuffer->Map(0, nullptr, reinterpret_cast<void**>(&pData));

    memcpy(
        pData,
        stateProps->GetShaderIdentifier(L"RayGen"),
        shaderIdSize
        );

    memcpy(
        pData + raygenTableSize,
        stateProps->GetShaderIdentifier(L"Miss"),
        shaderIdSize
        );

    memcpy(
        pData + raygenTableSize + missTableSize,
        stateProps->GetShaderIdentifier(L"HitGroup"),
        shaderIdSize
        );

    sbtBuffer->Unmap(0, nullptr);

    return true;
}

bool DXRendererRayTracing::createTriangleBuffers()
{
    struct Vertex
    {
        float x, y, z;
    };

    Vertex vertices[] =
        {
            {  0.0f,  0.5f, 0.0f },
            {  0.5f, -0.5f, 0.0f },
            { -0.5f, -0.5f, 0.0f }
        };

    uint32_t indices[] = { 0, 1, 2 };

    UINT vbSize = sizeof(vertices);
    UINT ibSize = sizeof(indices);

    D3D12_HEAP_PROPERTIES heap = {};
    heap.Type = D3D12_HEAP_TYPE_UPLOAD;

    D3D12_RESOURCE_DESC vbDesc = {};
    vbDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    vbDesc.Width = vbSize;
    vbDesc.Height = 1;
    vbDesc.DepthOrArraySize = 1;
    vbDesc.MipLevels = 1;
    vbDesc.SampleDesc.Count = 1;
    vbDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    if (FAILED(device->CreateCommittedResource(
            &heap,
            D3D12_HEAP_FLAG_NONE,
            &vbDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&vertexBuffer))))
        return false;

    void* pData = nullptr;
    vertexBuffer->Map(0, nullptr, &pData);
    memcpy(pData, vertices, vbSize);
    vertexBuffer->Unmap(0, nullptr);

    D3D12_RESOURCE_DESC ibDesc = vbDesc;
    ibDesc.Width = ibSize;

    if (FAILED(device->CreateCommittedResource(
            &heap,
            D3D12_HEAP_FLAG_NONE,
            &ibDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&indexBuffer))))
        return false;

    indexBuffer->Map(0, nullptr, &pData);
    memcpy(pData, indices, ibSize);
    indexBuffer->Unmap(0, nullptr);

    return true;
}

bool DXRendererRayTracing::createBLAS(ID3D12GraphicsCommandList4* cmd)
{
    D3D12_RAYTRACING_GEOMETRY_DESC geom = {};
    geom.Type = D3D12_RAYTRACING_GEOMETRY_TYPE_TRIANGLES;
    geom.Flags = D3D12_RAYTRACING_GEOMETRY_FLAG_OPAQUE;

    geom.Triangles.VertexBuffer.StartAddress = vertexBuffer->GetGPUVirtualAddress();
    geom.Triangles.VertexBuffer.StrideInBytes = sizeof(float) * 3;
    geom.Triangles.VertexCount = 3;
    geom.Triangles.VertexFormat = DXGI_FORMAT_R32G32B32_FLOAT;

    geom.Triangles.IndexBuffer = indexBuffer->GetGPUVirtualAddress();
    geom.Triangles.IndexCount = 3;
    geom.Triangles.IndexFormat = DXGI_FORMAT_R32_UINT;

    D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_INPUTS inputs = {};
    inputs.Type = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL;
    inputs.DescsLayout = D3D12_ELEMENTS_LAYOUT_ARRAY;
    inputs.NumDescs = 1;
    inputs.pGeometryDescs = &geom;
    inputs.Flags = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_BUILD_FLAG_NONE;

    D3D12_RAYTRACING_ACCELERATION_STRUCTURE_PREBUILD_INFO info = {};
    device5->GetRaytracingAccelerationStructurePrebuildInfo(&inputs, &info);

    if (info.ResultDataMaxSizeInBytes == 0)
        return false;

    D3D12_HEAP_PROPERTIES heap = {};
    heap.Type = D3D12_HEAP_TYPE_DEFAULT;

    D3D12_RESOURCE_DESC bufferDesc = {};
    bufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    bufferDesc.Width = info.ResultDataMaxSizeInBytes;
    bufferDesc.Height = 1;
    bufferDesc.DepthOrArraySize = 1;
    bufferDesc.MipLevels = 1;
    bufferDesc.SampleDesc.Count = 1;
    bufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    bufferDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;

    if (FAILED(device->CreateCommittedResource(
            &heap,
            D3D12_HEAP_FLAG_NONE,
            &bufferDesc,
            D3D12_RESOURCE_STATE_RAYTRACING_ACCELERATION_STRUCTURE,
            nullptr,
            IID_PPV_ARGS(&blasBuffer))))
        return false;

    bufferDesc.Width = info.ScratchDataSizeInBytes;

    if (FAILED(device->CreateCommittedResource(
            &heap,
            D3D12_HEAP_FLAG_NONE,
            &bufferDesc,
            D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
            nullptr,
            IID_PPV_ARGS(&blasScratch))))
        return false;

    D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_DESC build = {};
    build.Inputs = inputs;
    build.ScratchAccelerationStructureData = blasScratch->GetGPUVirtualAddress();
    build.DestAccelerationStructureData = blasBuffer->GetGPUVirtualAddress();

    cmd->BuildRaytracingAccelerationStructure(&build, 0, nullptr);

    D3D12_RESOURCE_BARRIER barrier = {};
    barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_UAV;
    barrier.UAV.pResource = blasBuffer;
    cmd->ResourceBarrier(1, &barrier);

    return true;
}

bool DXRendererRayTracing::createTLAS(ID3D12GraphicsCommandList4* cmd)
{
    D3D12_RAYTRACING_INSTANCE_DESC instance = {};
    instance.Transform[0][0] = 1.0f;
    instance.Transform[1][1] = 1.0f;
    instance.Transform[2][2] = 1.0f;
    instance.InstanceMask = 1;
    instance.AccelerationStructure = blasBuffer->GetGPUVirtualAddress();

    D3D12_HEAP_PROPERTIES uploadHeap = {};
    uploadHeap.Type = D3D12_HEAP_TYPE_UPLOAD;

    D3D12_RESOURCE_DESC instDesc = {};
    instDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    instDesc.Width = sizeof(D3D12_RAYTRACING_INSTANCE_DESC);
    instDesc.Height = 1;
    instDesc.DepthOrArraySize = 1;
    instDesc.MipLevels = 1;
    instDesc.SampleDesc.Count = 1;
    instDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;

    if (FAILED(device->CreateCommittedResource(
            &uploadHeap,
            D3D12_HEAP_FLAG_NONE,
            &instDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&tlasInstanceDesc))))
        return false;

    void* pData;
    tlasInstanceDesc->Map(0, nullptr, &pData);
    memcpy(pData, &instance, sizeof(instance));
    tlasInstanceDesc->Unmap(0, nullptr);

    D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_INPUTS inputs = {};
    inputs.Type = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL;
    inputs.DescsLayout = D3D12_ELEMENTS_LAYOUT_ARRAY;
    inputs.NumDescs = 1;
    inputs.InstanceDescs = tlasInstanceDesc->GetGPUVirtualAddress();

    D3D12_RAYTRACING_ACCELERATION_STRUCTURE_PREBUILD_INFO info = {};
    device5->GetRaytracingAccelerationStructurePrebuildInfo(&inputs, &info);

    D3D12_HEAP_PROPERTIES defaultHeap = {};
    defaultHeap.Type = D3D12_HEAP_TYPE_DEFAULT;

    D3D12_RESOURCE_DESC buffer = instDesc;
    buffer.Width = info.ResultDataMaxSizeInBytes;
    buffer.Flags = D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;

    if (FAILED(device->CreateCommittedResource(
            &defaultHeap,
            D3D12_HEAP_FLAG_NONE,
            &buffer,
            D3D12_RESOURCE_STATE_RAYTRACING_ACCELERATION_STRUCTURE,
            nullptr,
            IID_PPV_ARGS(&tlasBuffer))))
        return false;

    buffer.Width = info.ScratchDataSizeInBytes;

    if (FAILED(device->CreateCommittedResource(
            &defaultHeap,
            D3D12_HEAP_FLAG_NONE,
            &buffer,
            D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
            nullptr,
            IID_PPV_ARGS(&tlasScratch))))
        return false;

    D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_DESC build = {};
    build.Inputs = inputs;
    build.DestAccelerationStructureData = tlasBuffer->GetGPUVirtualAddress();
    build.ScratchAccelerationStructureData = tlasScratch->GetGPUVirtualAddress();

    cmd->BuildRaytracingAccelerationStructure(&build, 0, nullptr);

    D3D12_RESOURCE_BARRIER barrier = {};
    barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_UAV;
    barrier.UAV.pResource = tlasBuffer;
    cmd->ResourceBarrier(1, &barrier);

    return true;
}

void DXRendererRayTracing::createTLASSRV()
{
    UINT inc = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);

    D3D12_CPU_DESCRIPTOR_HANDLE cpu = descriptorHeap->GetCPUDescriptorHandleForHeapStart();
    cpu.ptr += inc;

    D3D12_SHADER_RESOURCE_VIEW_DESC srv = {};
    srv.ViewDimension = D3D12_SRV_DIMENSION_RAYTRACING_ACCELERATION_STRUCTURE;
    srv.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    srv.RaytracingAccelerationStructure.Location = tlasBuffer->GetGPUVirtualAddress();

    device->CreateShaderResourceView(nullptr, &srv, cpu);
}



void DXRendererRayTracing::renderFrame(ID3D12GraphicsCommandList4* cmd)
{
    updateCamera();
    if (!accelBuilt)
    {
        createBLAS(cmd);
        createTLAS(cmd);
        createTLASSRV();
        accelBuilt = true;
    }

    ID3D12DescriptorHeap* heaps[] = { descriptorHeap };
    cmd->SetDescriptorHeaps(1, heaps);
    cmd->SetPipelineState1(stateObject);
    cmd->SetComputeRootSignature(globalRootSignature);
    UINT inc = device->GetDescriptorHandleIncrementSize(
        D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV
        );

    D3D12_GPU_DESCRIPTOR_HANDLE gpu = descriptorHeap->GetGPUDescriptorHandleForHeapStart();

    cmd->SetComputeRootDescriptorTable(0, gpu);

    gpu.ptr += inc;
    cmd->SetComputeRootDescriptorTable(1, gpu);

    gpu.ptr += inc;
    cmd->SetComputeRootDescriptorTable(2, gpu);

    UINT idSize = D3D12_SHADER_IDENTIFIER_SIZE_IN_BYTES;
    UINT align = D3D12_RAYTRACING_SHADER_RECORD_BYTE_ALIGNMENT;
    UINT recordSize = (idSize + align - 1) & ~(align - 1);

    UINT tableAlign = D3D12_RAYTRACING_SHADER_TABLE_BYTE_ALIGNMENT;
    UINT raygenSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);
    UINT missSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);
    UINT hitSize = (recordSize + tableAlign - 1) & ~(tableAlign - 1);

    D3D12_GPU_VIRTUAL_ADDRESS addr = sbtBuffer->GetGPUVirtualAddress();

    D3D12_DISPATCH_RAYS_DESC d = {};
    d.RayGenerationShaderRecord.StartAddress = addr;
    d.RayGenerationShaderRecord.SizeInBytes = recordSize;

    d.MissShaderTable.StartAddress = addr + raygenSize;
    d.MissShaderTable.SizeInBytes = missSize;
    d.MissShaderTable.StrideInBytes = recordSize;

    d.HitGroupTable.StartAddress = addr + raygenSize + missSize;
    d.HitGroupTable.SizeInBytes = hitSize;
    d.HitGroupTable.StrideInBytes = recordSize;

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
