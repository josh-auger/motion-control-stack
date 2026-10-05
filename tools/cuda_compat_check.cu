#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cuda.h>
#include <cuda_runtime.h>

#include <dlfcn.h>
#include <limits.h>
#include <link.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

constexpr size_t kAllocationBytes = 1024 * 1024;
constexpr const char* kLogPath = "/tmp/share/cuda_compat_check.log";

class DiagnosticOutput {
 public:
  explicit DiagnosticOutput(const char* log_path) {
    log_.open(log_path, std::ios::out | std::ios::trunc);
    if (!log_.is_open()) {
      std::cerr << "WARNING: unable to open CUDA diagnostic logfile " << log_path
                << "; continuing with console output only\n";
    }
  }

  ~DiagnosticOutput() {
    Flush();
    if (log_.is_open()) {
      log_.close();
    }
  }

  DiagnosticOutput(const DiagnosticOutput&) = delete;
  DiagnosticOutput& operator=(const DiagnosticOutput&) = delete;

  template <typename T>
  DiagnosticOutput& operator<<(const T& value) {
    std::cout << value;
    if (log_.is_open()) {
      log_ << value;
      log_.flush();
    }
    return *this;
  }

  void Flush() {
    std::cout.flush();
    if (log_.is_open()) {
      log_.flush();
    }
  }

 private:
  std::ofstream log_;
};

struct FailureTracker {
  bool failed = false;
  std::string operation;
  long long code = 0;
  std::string detail;

  void Record(const std::string& failed_operation, long long failed_code,
              const std::string& failed_detail) {
    if (failed) {
      return;
    }
    failed = true;
    operation = failed_operation;
    code = failed_code;
    detail = failed_detail;
  }
};

std::string VersionString(int version) {
  if (version <= 0) {
    return "unknown";
  }

  const int major = version / 1000;
  const int minor = (version % 1000) / 10;
  const int patch = version % 10;
  std::ostringstream result;
  result << major << '.' << minor;
  if (patch != 0) {
    result << '.' << patch;
  }
  return result.str();
}

std::string RuntimeErrorDetail(cudaError_t result) {
  std::ostringstream detail;
  detail << cudaGetErrorName(result) << ": " << cudaGetErrorString(result);
  return detail.str();
}

bool ReportRuntimeResult(const char* operation, cudaError_t result,
                         FailureTracker* failures,
                         DiagnosticOutput* output) {
  const bool passed = result == cudaSuccess;
  const std::string detail = RuntimeErrorDetail(result);
  *output << (passed ? "[PASS] " : "[FAIL] ") << operation
          << ": code=" << static_cast<int>(result) << " (" << detail << ")\n";
  if (!passed) {
    failures->Record(operation, static_cast<int>(result), detail);
  }
  return passed;
}

using CuInitFunction = CUresult(CUDAAPI*)(unsigned int);
using CuDriverGetVersionFunction = CUresult(CUDAAPI*)(int*);
using CuGetErrorNameFunction = CUresult(CUDAAPI*)(CUresult, const char**);
using CuGetErrorStringFunction = CUresult(CUDAAPI*)(CUresult, const char**);

struct DriverFunctions {
  CuInitFunction init = nullptr;
  CuDriverGetVersionFunction get_version = nullptr;
  CuGetErrorNameFunction get_error_name = nullptr;
  CuGetErrorStringFunction get_error_string = nullptr;
};

std::string DriverErrorDetail(CUresult result, const DriverFunctions& driver) {
  const char* name = nullptr;
  const char* description = nullptr;
  if (driver.get_error_name != nullptr) {
    driver.get_error_name(result, &name);
  }
  if (driver.get_error_string != nullptr) {
    driver.get_error_string(result, &description);
  }

  std::ostringstream detail;
  detail << (name != nullptr ? name : "name unavailable") << ": "
         << (description != nullptr ? description : "description unavailable");
  return detail.str();
}

bool ReportDriverResult(const char* operation, CUresult result,
                        const DriverFunctions& driver,
                        FailureTracker* failures,
                        DiagnosticOutput* output) {
  const bool passed = result == CUDA_SUCCESS;
  const std::string detail = DriverErrorDetail(result, driver);
  *output << (passed ? "[PASS] " : "[FAIL] ") << operation
          << ": code=" << static_cast<int>(result) << " (" << detail << ")\n";
  if (!passed) {
    failures->Record(operation, static_cast<int>(result), detail);
  }
  return passed;
}

void* LoadDriverSymbol(void* handle, const char* symbol,
                       FailureTracker* failures,
                       DiagnosticOutput* output) {
  dlerror();
  void* address = dlsym(handle, symbol);
  const char* error = dlerror();
  if (error != nullptr) {
    *output << "[FAIL] dlsym(" << symbol << "): " << error << '\n';
    failures->Record(std::string("dlsym(") + symbol + ')', -1, error);
    return nullptr;
  }
  *output << "[PASS] dlsym(" << symbol << ")\n";
  return address;
}

void ReportDriverLibraryPath(void* handle, void* driver_symbol,
                             DiagnosticOutput* output) {
  std::string loader_path;
  struct link_map* map = nullptr;
  if (dlinfo(handle, RTLD_DI_LINKMAP, &map) == 0 && map != nullptr &&
      map->l_name != nullptr && map->l_name[0] != '\0') {
    loader_path = map->l_name;
    *output << "libcuda dlinfo path: " << loader_path << '\n';
  } else {
    const char* error = dlerror();
    *output << "libcuda dlinfo path: unavailable"
            << (error != nullptr ? std::string(" (") + error + ')' : "") << '\n';
  }

  Dl_info info{};
  if (driver_symbol != nullptr && dladdr(driver_symbol, &info) != 0 &&
      info.dli_fname != nullptr) {
    loader_path = info.dli_fname;
    *output << "libcuda symbol path: " << loader_path << '\n';
  } else {
    *output << "libcuda symbol path: unavailable\n";
  }

  if (!loader_path.empty()) {
    char canonical_path[PATH_MAX];
    if (realpath(loader_path.c_str(), canonical_path) != nullptr) {
      *output << "libcuda canonical path: " << canonical_path << '\n';
    } else {
      *output << "libcuda canonical path: unavailable for " << loader_path << '\n';
    }
  }
}

__global__ void CompatibilityKernel(unsigned char* memory) {
  if (blockIdx.x == 0 && threadIdx.x == 0) {
    memory[0] = 0x5a;
  }
}

}  // namespace

int main() {
  DiagnosticOutput output(kLogPath);
  FailureTracker failures;

  output << "CUDA compatibility diagnostic\n"
         << "========================================\n"
         << "Compile-time CUDART_VERSION: " << CUDART_VERSION << " ("
         << VersionString(CUDART_VERSION) << ")\n";

  output << "\nCUDA Driver API and libcuda identity\n";
  void* driver_handle = dlopen("libcuda.so.1", RTLD_NOW | RTLD_LOCAL);
  DriverFunctions driver;
  if (driver_handle == nullptr) {
    const char* error = dlerror();
    const std::string detail = error != nullptr ? error : "unknown dlopen error";
    output << "[FAIL] dlopen(libcuda.so.1): " << detail << '\n';
    failures.Record("dlopen(libcuda.so.1)", -1, detail);
  } else {
    output << "[PASS] dlopen(libcuda.so.1)\n";
    void* init_symbol =
        LoadDriverSymbol(driver_handle, "cuInit", &failures, &output);
    driver.init = reinterpret_cast<CuInitFunction>(init_symbol);
    driver.get_version = reinterpret_cast<CuDriverGetVersionFunction>(
        LoadDriverSymbol(driver_handle, "cuDriverGetVersion", &failures, &output));
    driver.get_error_name = reinterpret_cast<CuGetErrorNameFunction>(
        LoadDriverSymbol(driver_handle, "cuGetErrorName", &failures, &output));
    driver.get_error_string = reinterpret_cast<CuGetErrorStringFunction>(
        LoadDriverSymbol(driver_handle, "cuGetErrorString", &failures, &output));

    ReportDriverLibraryPath(driver_handle, init_symbol, &output);

    if (driver.init != nullptr) {
      ReportDriverResult("cuInit(0)", driver.init(0), driver, &failures, &output);
    }

    if (driver.get_version != nullptr) {
      int direct_driver_version = 0;
      const CUresult result = driver.get_version(&direct_driver_version);
      if (ReportDriverResult("cuDriverGetVersion", result, driver, &failures,
                             &output)) {
        output << "Direct cuDriverGetVersion value: " << direct_driver_version
               << " (" << VersionString(direct_driver_version) << ")\n";
      }
    }
  }

  output << "\nCUDA Runtime API identity\n";
  int runtime_version = 0;
  cudaError_t runtime_result = cudaRuntimeGetVersion(&runtime_version);
  if (ReportRuntimeResult("cudaRuntimeGetVersion", runtime_result, &failures,
                          &output)) {
    output << "cudaRuntimeGetVersion value: " << runtime_version << " ("
           << VersionString(runtime_version) << ")\n";
  }

  int runtime_driver_version = 0;
  runtime_result = cudaDriverGetVersion(&runtime_driver_version);
  if (ReportRuntimeResult("cudaDriverGetVersion", runtime_result, &failures,
                          &output)) {
    output << "cudaDriverGetVersion value: " << runtime_driver_version << " ("
           << VersionString(runtime_driver_version) << ")\n";
    if (runtime_driver_version == 0) {
      output << "[FAIL] cudaDriverGetVersion reported no installed driver\n";
      failures.Record("cudaDriverGetVersion (no installed driver)", 0,
                      "cudaSuccess, but the reported driver version is 0");
    }
  }

  output << "\nCUDA device discovery\n";
  int device_count = 0;
  const cudaError_t count_result = cudaGetDeviceCount(&device_count);
  if (ReportRuntimeResult("cudaGetDeviceCount", count_result, &failures,
                          &output)) {
    output << "CUDA device count: " << device_count << '\n';
    if (device_count == 0) {
      output << "[FAIL] CUDA device discovery: no CUDA devices found\n";
      failures.Record("cudaGetDeviceCount (no devices)", 0,
                      "cudaSuccess, but no CUDA devices were found");
    }
  }

  if (count_result == cudaSuccess) {
    for (int device = 0; device < device_count; ++device) {
      cudaDeviceProp properties{};
      std::ostringstream properties_operation;
      properties_operation << "cudaGetDeviceProperties(" << device << ')';
      const cudaError_t properties_result =
          cudaGetDeviceProperties(&properties, device);
      if (!ReportRuntimeResult(properties_operation.str().c_str(), properties_result,
                               &failures, &output)) {
        continue;
      }

      output << "Device " << device << ": name=\"" << properties.name
             << "\", compute capability=" << properties.major << '.'
             << properties.minor << ", global memory="
             << properties.totalGlobalMem << " bytes ("
             << properties.totalGlobalMem / (1024 * 1024) << " MiB)\n";

      char pci_bus_id[32] = {};
      std::ostringstream pci_operation;
      pci_operation << "cudaDeviceGetPCIBusId(" << device << ')';
      const cudaError_t pci_result =
          cudaDeviceGetPCIBusId(pci_bus_id, sizeof(pci_bus_id), device);
      if (ReportRuntimeResult(pci_operation.str().c_str(), pci_result, &failures,
                              &output)) {
        output << "Device " << device << " PCI bus ID: " << pci_bus_id << '\n';
      }
    }
  }

  output << "\nGPU allocation and operation test (device 0, "
         << kAllocationBytes << " bytes)\n";
  if (count_result == cudaSuccess && device_count > 0) {
    const cudaError_t set_device_result = cudaSetDevice(0);
    if (ReportRuntimeResult("cudaSetDevice(0)", set_device_result, &failures,
                            &output)) {
      void* device_memory = nullptr;
      const cudaError_t malloc_result =
          cudaMalloc(&device_memory, kAllocationBytes);
      if (ReportRuntimeResult("cudaMalloc(1 MiB)", malloc_result, &failures,
                              &output)) {
        ReportRuntimeResult("cudaMemset(1 MiB)",
                            cudaMemset(device_memory, 0xa5, kAllocationBytes),
                            &failures, &output);

        CompatibilityKernel<<<1, 1>>>(
            static_cast<unsigned char*>(device_memory));
        ReportRuntimeResult("CompatibilityKernel launch", cudaGetLastError(),
                            &failures, &output);
        ReportRuntimeResult("cudaDeviceSynchronize", cudaDeviceSynchronize(),
                            &failures, &output);
        ReportRuntimeResult("cudaFree", cudaFree(device_memory), &failures,
                            &output);
      } else {
        output << "[SKIP] cudaMemset, kernel, synchronize, and free: allocation "
                  "did not succeed\n";
      }
    } else {
      output << "[SKIP] allocation and GPU operations: cudaSetDevice(0) did not "
                "succeed\n";
    }
  } else {
    output << "[SKIP] allocation and GPU operations: no usable CUDA device was "
              "discovered\n";
  }

  if (driver_handle != nullptr) {
    dlclose(driver_handle);
  }

  output << "\n========================================\n";
  if (!failures.failed) {
    output << "CUDA COMPATIBILITY CHECK: PASS\n"
           << "========================================\n";
    return 0;
  }

  output << "CUDA COMPATIBILITY CHECK: FAIL\n"
         << "First failing CUDA operation: " << failures.operation << '\n'
         << "Error: " << failures.code << ' ' << failures.detail << '\n'
         << "========================================\n";
  return 1;
}
