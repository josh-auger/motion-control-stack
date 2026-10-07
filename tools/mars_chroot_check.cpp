#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <dlfcn.h>
#include <link.h>
#include <sys/stat.h>
#include <sys/sysmacros.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr const char* kLogPath = "/tmp/share/mars_chroot_check.log";
constexpr unsigned int kDefaultHoldSeconds = 300;
constexpr size_t kMaximumEntriesPerLibraryRoot = 200000;

class DiagnosticOutput {
 public:
  explicit DiagnosticOutput(const char* log_path) {
    log_.open(log_path, std::ios::out | std::ios::trunc);
    if (!log_.is_open()) {
      std::cerr << "WARNING: unable to open MARS chroot diagnostic logfile "
                << log_path << "; continuing with console output only\n";
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

  bool LogIsOpen() const { return log_.is_open(); }

 private:
  std::ofstream log_;
};

std::string Lowercase(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char character) {
                   return static_cast<char>(std::tolower(character));
                 });
  return value;
}

std::string FileType(mode_t mode) {
  if (S_ISREG(mode)) return "regular file";
  if (S_ISDIR(mode)) return "directory";
  if (S_ISLNK(mode)) return "symbolic link";
  if (S_ISCHR(mode)) return "character device";
  if (S_ISBLK(mode)) return "block device";
  if (S_ISFIFO(mode)) return "FIFO";
  if (S_ISSOCK(mode)) return "socket";
  return "unknown";
}

std::string PermissionBits(mode_t mode) {
  std::string bits = "---------";
  constexpr mode_t masks[] = {S_IRUSR, S_IWUSR, S_IXUSR, S_IRGRP, S_IWGRP,
                              S_IXGRP, S_IROTH, S_IWOTH, S_IXOTH};
  constexpr char characters[] = {'r', 'w', 'x', 'r', 'w', 'x', 'r', 'w', 'x'};
  for (size_t index = 0; index < 9; ++index) {
    if ((mode & masks[index]) != 0) bits[index] = characters[index];
  }
  if ((mode & S_ISUID) != 0) bits[2] = (mode & S_IXUSR) != 0 ? 's' : 'S';
  if ((mode & S_ISGID) != 0) bits[5] = (mode & S_IXGRP) != 0 ? 's' : 'S';
  if ((mode & S_ISVTX) != 0) bits[8] = (mode & S_IXOTH) != 0 ? 't' : 'T';

  std::ostringstream result;
  result << '0' << std::oct << std::setw(4) << std::setfill('0')
         << (mode & 07777) << std::dec << " (" << bits << ')';
  return result.str();
}

bool ReadLink(const std::string& path, std::string* target, std::string* error) {
  std::vector<char> buffer(4096);
  const ssize_t length = readlink(path.c_str(), buffer.data(), buffer.size() - 1);
  if (length < 0) {
    *error = std::strerror(errno);
    return false;
  }
  buffer[static_cast<size_t>(length)] = '\0';
  *target = buffer.data();
  return true;
}

std::string CanonicalPath(const std::string& path) {
  char* resolved = realpath(path.c_str(), nullptr);
  if (resolved == nullptr) return "<unavailable: " + std::string(std::strerror(errno)) + '>';
  const std::string result(resolved);
  std::free(resolved);
  return result;
}

std::string PathDetails(const std::string& path) {
  struct stat information {};
  if (lstat(path.c_str(), &information) != 0) {
    const int error = errno;
    if (error == EACCES || error == EPERM) {
      return "PERMISSION DENIED (" + std::string(std::strerror(error)) + ')';
    }
    if (error == ENOENT || error == ENOTDIR) {
      return "MISSING (" + std::string(std::strerror(error)) + ')';
    }
    return "ERROR (" + std::string(std::strerror(error)) + ')';
  }

  std::ostringstream details;
  details << "PRESENT | type=" << FileType(information.st_mode)
          << " | permissions=" << PermissionBits(information.st_mode)
          << " | uid=" << information.st_uid << " | gid=" << information.st_gid
          << " | size=" << information.st_size;

  if (S_ISCHR(information.st_mode) || S_ISBLK(information.st_mode)) {
    details << " | device=" << major(information.st_rdev) << ':'
            << minor(information.st_rdev);
  }

  if (S_ISLNK(information.st_mode)) {
    std::string target;
    std::string error;
    if (ReadLink(path, &target, &error)) {
      details << " | symlink_target=" << target;
    } else {
      details << " | symlink_target=<unavailable: " << error << '>';
    }

    struct stat resolved_information {};
    if (stat(path.c_str(), &resolved_information) == 0) {
      details << " | resolved_size=" << resolved_information.st_size;
    }
  }

  details << " | canonical=" << CanonicalPath(path);
  return details.str();
}

std::string CurrentWorkingDirectory() {
  char* directory = getcwd(nullptr, 0);
  if (directory == nullptr) {
    return "<unavailable: " + std::string(std::strerror(errno)) + '>';
  }
  const std::string result(directory);
  std::free(directory);
  return result;
}

std::string Hostname() {
  std::vector<char> buffer(256, '\0');
  if (gethostname(buffer.data(), buffer.size() - 1) != 0) {
    return "<unavailable: " + std::string(std::strerror(errno)) + '>';
  }
  return buffer.data();
}

void ReportProcessIdentity(DiagnosticOutput* output) {
  *output << "\nProcess identity\n----------------\n"
          << "PID: " << getpid() << '\n'
          << "PPID: " << getppid() << '\n'
          << "UID/EUID: " << getuid() << '/' << geteuid() << '\n'
          << "GID/EGID: " << getgid() << '/' << getegid() << '\n'
          << "Hostname: " << Hostname() << '\n'
          << "Working directory: " << CurrentWorkingDirectory() << '\n';

  for (const char* path : {"/proc/self/exe", "/proc/self/root", "/proc/self/cwd"}) {
    std::string target;
    std::string error;
    if (ReadLink(path, &target, &error)) {
      *output << path << ": " << target << '\n';
    } else {
      *output << path << ": <unavailable: " << error << ">\n";
    }
  }
}

std::string ReportNamespaces(DiagnosticOutput* output) {
  *output << "\nNamespace identifiers\n---------------------\n";
  std::string mount_namespace = "<unavailable>";
  for (const char* name : {"mnt", "pid", "user", "net", "ipc", "uts"}) {
    const std::string path = std::string("/proc/self/ns/") + name;
    std::string target;
    std::string error;
    if (ReadLink(path, &target, &error)) {
      *output << path << ": " << target << '\n';
      if (std::string(name) == "mnt") mount_namespace = target;
    } else {
      *output << path << ": <unavailable: " << error << ">\n";
    }
  }
  return mount_namespace;
}

bool ContainsEnvironmentKeyword(const std::string& name) {
  const std::string lowercase_name = Lowercase(name);
  for (const char* keyword : {"cuda", "nvidia", "fire", "chroot", "mars"}) {
    if (lowercase_name.find(keyword) != std::string::npos) return true;
  }
  return false;
}

void ReportEnvironment(DiagnosticOutput* output) {
  std::map<std::string, std::string> environment;
  for (char** entry = environ; entry != nullptr && *entry != nullptr; ++entry) {
    const std::string value(*entry);
    const size_t separator = value.find('=');
    if (separator != std::string::npos) {
      environment[value.substr(0, separator)] = value.substr(separator + 1);
    }
  }

  const std::vector<std::string> required = {
      "PATH",          "LD_LIBRARY_PATH",          "CUDA_HOME",
      "CUDA_PATH",     "CUDA_VISIBLE_DEVICES",     "NVIDIA_VISIBLE_DEVICES",
      "NVIDIA_DRIVER_CAPABILITIES", "WORKDIR"};
  std::set<std::string> reported;

  *output << "\nRelevant environment\n--------------------\n";
  for (const std::string& name : required) {
    const auto value = environment.find(name);
    *output << name << '='
            << (value == environment.end() ? "<unset>" : value->second) << '\n';
    reported.insert(name);
  }
  for (const auto& variable : environment) {
    if (reported.count(variable.first) == 0 &&
        ContainsEnvironmentKeyword(variable.first)) {
      *output << variable.first << '=' << variable.second << '\n';
    }
  }
}

void ReportImportantPaths(DiagnosticOutput* output) {
  *output << "\nImportant filesystem paths\n--------------------------\n";
  const std::vector<std::string> paths = {
      "/",          "/lib",             "/lib64",       "/usr",
      "/usr/lib",   "/usr/lib64",       "/usr/local",   "/usr/local/cuda",
      "/usr/local/cuda-13", "/usr/local/cuda-13.2", "/tmp", "/tmp/share",
      "/dev",       "/proc",            "/sys"};
  for (const std::string& path : paths) {
    *output << path << ": " << PathDetails(path) << '\n';
  }
}

bool IsRelevantLibraryName(const std::string& name) {
  return name.rfind("libcuda.so", 0) == 0 ||
         name.rfind("libcudart.so", 0) == 0 ||
         (name.rfind("libnvidia-", 0) == 0 && name.find(".so") != std::string::npos);
}

std::vector<std::string> FindRelevantLibraries(DiagnosticOutput* output) {
  const std::vector<std::string> roots = {"/lib", "/lib64", "/usr/lib",
                                          "/usr/lib64", "/usr/local", "/opt"};
  std::set<std::string> libraries;

  for (const std::string& root : roots) {
    struct stat information {};
    if (stat(root.c_str(), &information) != 0) {
      *output << "Search root " << root << ": unavailable ("
              << std::strerror(errno) << ")\n";
      continue;
    }
    if (!S_ISDIR(information.st_mode)) {
      *output << "Search root " << root << ": not a directory\n";
      continue;
    }

    std::error_code error;
    std::filesystem::recursive_directory_iterator iterator(
        root, std::filesystem::directory_options::skip_permission_denied, error);
    const std::filesystem::recursive_directory_iterator end;
    size_t entries_examined = 0;
    while (iterator != end && entries_examined < kMaximumEntriesPerLibraryRoot) {
      if (!error) {
        const std::filesystem::path path = iterator->path();
        if (IsRelevantLibraryName(path.filename().string())) {
          libraries.insert(path.string());
        }
      }
      ++entries_examined;
      iterator.increment(error);
      if (error) {
        *output << "Library scan warning under " << root << ": "
                << error.message() << '\n';
        error.clear();
      }
    }
    if (entries_examined == kMaximumEntriesPerLibraryRoot && iterator != end) {
      *output << "Library scan under " << root << " stopped after "
              << kMaximumEntriesPerLibraryRoot << " entries\n";
    }
  }

  return std::vector<std::string>(libraries.begin(), libraries.end());
}

bool ReportLibraryInventory(DiagnosticOutput* output) {
  *output << "\nNVIDIA/CUDA library inventory\n-----------------------------\n"
          << "Expected host driver paths:\n";
  for (const char* path : {"/lib/x86_64-linux-gnu/libcuda.so",
                           "/lib/x86_64-linux-gnu/libcuda.so.1",
                           "/lib/x86_64-linux-gnu/libcuda.so.510.47.03"}) {
    *output << "  " << path << ": " << PathDetails(path) << '\n';
  }

  *output << "Discovered matching libraries:\n";
  const std::vector<std::string> libraries = FindRelevantLibraries(output);
  bool libcuda_visible = false;
  if (libraries.empty()) {
    *output << "  <none found>\n";
  }
  for (const std::string& path : libraries) {
    *output << "  " << path << ": " << PathDetails(path) << '\n';
    if (std::filesystem::path(path).filename() == "libcuda.so.1") {
      libcuda_visible = true;
    }
  }
  *output << "Matching library count: " << libraries.size() << '\n';
  return libcuda_visible;
}

void ReportTextFile(const std::string& path, DiagnosticOutput* output) {
  std::ifstream input(path);
  if (!input.is_open()) {
    *output << path << ": <unavailable: " << std::strerror(errno) << ">\n";
    return;
  }
  *output << path << ":\n";
  std::string line;
  bool had_lines = false;
  while (std::getline(input, line)) {
    *output << "  " << line << '\n';
    had_lines = true;
  }
  if (!had_lines) *output << "  <empty>\n";
}

void ReportLoaderConfiguration(DiagnosticOutput* output) {
  *output << "\nDynamic loader configuration\n----------------------------\n"
          << "/etc/ld.so.cache: " << PathDetails("/etc/ld.so.cache") << '\n';
  ReportTextFile("/etc/ld.so.conf", output);

  std::vector<std::string> configuration_files;
  std::error_code error;
  std::filesystem::directory_iterator iterator(
      "/etc/ld.so.conf.d", std::filesystem::directory_options::skip_permission_denied,
      error);
  const std::filesystem::directory_iterator end;
  while (!error && iterator != end) {
    if (iterator->path().extension() == ".conf") {
      configuration_files.push_back(iterator->path().string());
    }
    iterator.increment(error);
  }
  if (error) {
    *output << "/etc/ld.so.conf.d: <unavailable: " << error.message() << ">\n";
  }
  std::sort(configuration_files.begin(), configuration_files.end());
  for (const std::string& path : configuration_files) {
    ReportTextFile(path, output);
  }
}

struct LoaderResult {
  bool loaded = false;
  std::string resolved_path = "<not loaded>";
};

LoaderResult TestDynamicLoad(const char* library, DiagnosticOutput* output) {
  LoaderResult result;
  dlerror();
  void* handle = dlopen(library, RTLD_NOW | RTLD_LOCAL);
  if (handle == nullptr) {
    const char* error = dlerror();
    *output << "[FAIL] dlopen(\"" << library << "\"): "
            << (error != nullptr ? error : "unknown dlopen error") << '\n';
    return result;
  }

  result.loaded = true;
  *output << "[PASS] dlopen(\"" << library << "\")\n";
  struct link_map* map = nullptr;
  dlerror();
  if (dlinfo(handle, RTLD_DI_LINKMAP, &map) == 0 && map != nullptr &&
      map->l_name != nullptr && map->l_name[0] != '\0') {
    result.resolved_path = map->l_name;
    *output << "  dlinfo path: " << result.resolved_path << '\n'
            << "  canonical path: " << CanonicalPath(result.resolved_path) << '\n';
  } else {
    const char* error = dlerror();
    *output << "  resolved path unavailable"
            << (error != nullptr ? std::string(": ") + error : "") << '\n';
  }
  dlclose(handle);
  return result;
}

LoaderResult ReportDynamicLoaderTests(DiagnosticOutput* output) {
  *output << "\nDynamic libcuda resolution (no CUDA API calls)\n"
          << "----------------------------------------------\n";
  const LoaderResult primary = TestDynamicLoad("libcuda.so.1", output);
  TestDynamicLoad("libcuda.so", output);
  return primary;
}

enum class DevicePresence { kPresent, kMissing, kPermissionDenied, kError };

DevicePresence ReportDevicePath(const std::string& path, DiagnosticOutput* output) {
  struct stat information {};
  if (lstat(path.c_str(), &information) == 0) {
    *output << path << ": " << PathDetails(path) << '\n';
    return DevicePresence::kPresent;
  }
  if (errno == EACCES || errno == EPERM) {
    *output << path << ": PERMISSION DENIED (" << std::strerror(errno) << ")\n";
    return DevicePresence::kPermissionDenied;
  }
  if (errno == ENOENT || errno == ENOTDIR) {
    *output << path << ": MISSING\n";
    return DevicePresence::kMissing;
  }
  *output << path << ": ERROR (" << std::strerror(errno) << ")\n";
  return DevicePresence::kError;
}

struct DeviceInventory {
  std::map<std::string, DevicePresence> state;
  size_t primary_present = 0;
};

DeviceInventory ReportDeviceInventory(DiagnosticOutput* output) {
  *output << "\nNVIDIA device-node inventory\n----------------------------\n";
  DeviceInventory inventory;
  const std::vector<std::string> primary = {"/dev/nvidia0", "/dev/nvidiactl",
                                             "/dev/nvidia-uvm",
                                             "/dev/nvidia-uvm-tools"};
  for (const std::string& path : primary) {
    inventory.state[path] = ReportDevicePath(path, output);
    if (inventory.state[path] == DevicePresence::kPresent) {
      ++inventory.primary_present;
    }
  }

  inventory.state["/dev/nvidia-caps"] =
      ReportDevicePath("/dev/nvidia-caps", output);
  if (inventory.state["/dev/nvidia-caps"] == DevicePresence::kPresent) {
    *output << "/dev/nvidia-caps immediate contents:\n";
    std::vector<std::string> entries;
    std::error_code error;
    std::filesystem::directory_iterator iterator(
        "/dev/nvidia-caps",
        std::filesystem::directory_options::skip_permission_denied, error);
    const std::filesystem::directory_iterator end;
    while (!error && iterator != end) {
      entries.push_back(iterator->path().string());
      iterator.increment(error);
    }
    if (error) {
      *output << "  PERMISSION DENIED/ERROR: " << error.message() << '\n';
    }
    std::sort(entries.begin(), entries.end());
    if (entries.empty() && !error) *output << "  <empty>\n";
    for (const std::string& path : entries) {
      *output << "  " << path << ": " << PathDetails(path) << '\n';
    }
  }
  return inventory;
}

bool IsAtOrBelow(const std::string& path, const std::string& root) {
  if (root == "/") return path == "/";
  return path == root ||
         (path.size() > root.size() && path.compare(0, root.size(), root) == 0 &&
          path[root.size()] == '/');
}

bool IsRelevantMountLine(const std::string& line) {
  std::istringstream fields(line);
  std::string mount_id;
  std::string parent_id;
  std::string device;
  std::string root;
  std::string mount_point;
  fields >> mount_id >> parent_id >> device >> root >> mount_point;
  for (const char* relevant_root : {"/", "/dev", "/proc", "/sys", "/tmp",
                                    "/lib", "/usr", "/opt"}) {
    if (IsAtOrBelow(mount_point, relevant_root)) return true;
  }
  const std::string lowercase_line = Lowercase(line);
  for (const char* keyword : {"cuda", "nvidia", "fire", "chroot"}) {
    if (lowercase_line.find(keyword) != std::string::npos) return true;
  }
  return false;
}

void ReportMountInformation(DiagnosticOutput* output) {
  *output << "\nRelevant mount information (/proc/self/mountinfo)\n"
          << "-------------------------------------------------\n";
  std::ifstream input("/proc/self/mountinfo");
  if (!input.is_open()) {
    *output << "<unavailable: " << std::strerror(errno) << ">\n";
    return;
  }
  std::string line;
  size_t relevant_lines = 0;
  while (std::getline(input, line)) {
    if (IsRelevantMountLine(line)) {
      *output << line << '\n';
      ++relevant_lines;
    }
  }
  if (relevant_lines == 0) *output << "<no relevant mount entries found>\n";
}

void ReportProcStatus(DiagnosticOutput* output) {
  *output << "\nRelevant /proc/self/status fields\n"
          << "--------------------------------\n";
  const std::vector<std::string> required = {
      "Name", "Pid",  "PPid", "Uid",      "Gid",    "Groups",
      "CapInh", "CapPrm", "CapEff", "CapBnd", "NoNewPrivs", "Seccomp"};
  std::set<std::string> required_set(required.begin(), required.end());
  std::map<std::string, std::string> values;
  std::ifstream input("/proc/self/status");
  if (!input.is_open()) {
    *output << "<unavailable: " << std::strerror(errno) << ">\n";
    return;
  }
  std::string line;
  while (std::getline(input, line)) {
    const size_t separator = line.find(':');
    if (separator != std::string::npos) {
      const std::string name = line.substr(0, separator);
      if (required_set.count(name) != 0) values[name] = line;
    }
  }
  for (const std::string& name : required) {
    const auto value = values.find(name);
    *output << (value == values.end() ? name + ": <unavailable>" : value->second)
            << '\n';
  }
}

std::string PresenceText(DevicePresence presence) {
  switch (presence) {
    case DevicePresence::kPresent:
      return "PRESENT";
    case DevicePresence::kMissing:
      return "MISSING";
    case DevicePresence::kPermissionDenied:
      return "PERMISSION DENIED";
    case DevicePresence::kError:
      return "ERROR";
  }
  return "UNKNOWN";
}

void ReportSummary(bool libcuda_visible, const LoaderResult& loader,
                   const DeviceInventory& devices,
                   const std::string& mount_namespace, pid_t pid,
                   bool share_writable, DiagnosticOutput* output) {
  *output << "\n========================================\n"
          << "MARS/FIRE CHROOT DIAGNOSTIC SUMMARY\n"
          << "========================================\n"
          << "libcuda.so.1 visible: " << (libcuda_visible ? "YES" : "NO") << '\n'
          << "libcuda.so.1 dlopen: " << (loader.loaded ? "PASS" : "FAIL") << '\n'
          << "libcuda.so.1 resolved path: " << loader.resolved_path << '\n'
          << "NVIDIA primary device nodes visible: " << devices.primary_present
          << "/4\n";
  for (const char* path : {"/dev/nvidia0", "/dev/nvidiactl", "/dev/nvidia-uvm"}) {
    const auto state = devices.state.find(path);
    *output << path << ": "
            << (state == devices.state.end() ? "UNKNOWN"
                                             : PresenceText(state->second))
            << '\n';
  }
  *output << "/tmp/share writable: " << (share_writable ? "YES" : "NO") << '\n'
          << "Mount namespace: " << mount_namespace << '\n'
          << "PID: " << pid << '\n'
          << "========================================\n";
}

bool ParseHoldSeconds(int argc, char** argv, unsigned int* hold_seconds,
                      DiagnosticOutput* output) {
  *hold_seconds = kDefaultHoldSeconds;
  if (argc == 1) return true;
  if (argc != 3 || std::string(argv[1]) != "--hold-seconds") {
    *output << "Usage: " << argv[0] << " [--hold-seconds SECONDS]\n";
    return false;
  }

  errno = 0;
  char* end = nullptr;
  const long value = std::strtol(argv[2], &end, 10);
  if (errno != 0 || end == argv[2] || *end != '\0' || value < 0 || value > 86400) {
    *output << "Invalid hold duration: " << argv[2]
            << " (expected 0-86400 seconds)\n";
    return false;
  }
  *hold_seconds = static_cast<unsigned int>(value);
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  DiagnosticOutput output(kLogPath);
  unsigned int hold_seconds = kDefaultHoldSeconds;
  if (!ParseHoldSeconds(argc, argv, &hold_seconds, &output)) return 2;

  const pid_t pid = getpid();
  output << "MARS/FIRE chroot environment diagnostic\n"
         << "========================================\n";

  ReportProcessIdentity(&output);
  const std::string mount_namespace = ReportNamespaces(&output);
  ReportEnvironment(&output);
  ReportImportantPaths(&output);
  const bool libcuda_visible = ReportLibraryInventory(&output);
  ReportLoaderConfiguration(&output);
  const LoaderResult loader = ReportDynamicLoaderTests(&output);
  const DeviceInventory devices = ReportDeviceInventory(&output);
  ReportMountInformation(&output);
  ReportProcStatus(&output);
  ReportSummary(libcuda_visible, loader, devices, mount_namespace, pid,
                output.LogIsOpen(), &output);

  output << "\nDiagnostic collection complete.\n"
         << "Holding process for " << hold_seconds
         << " seconds for host-side inspection.\n"
         << "PID: " << pid << '\n';
  output.Flush();

  std::this_thread::sleep_for(std::chrono::seconds(hold_seconds));

  output << "Inspection hold complete. Exiting normally.\n";
  output.Flush();
  return 0;
}
