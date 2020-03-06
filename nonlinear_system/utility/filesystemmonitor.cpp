#include "filesystemmonitor.h"
#include "inotify.h"
#include <future>
#include <chrono>


static std::future<bool> _threadHandle;
static bool isInitialized = false;
static int sweepTime = 10; // 10ms correspond to 120fps

void FileSystemMonitor::Init(const char* path)
{
    try
    {
        Inotify::Init(path);
        //_threadHandle = std::async(std::launch::async,Inotify::Update);
        _threadHandle = std::async(std::launch::async, [] () {
                // Use sleep_for to wait specified time (or sleep_until).
                std::this_thread::sleep_for( std::chrono::milliseconds{sweepTime});
                // Do whatever you want.
                return Inotify::Update();
            } );
        isInitialized = true;
    }catch(...)
    {
        isInitialized = false;
    }
}

bool FileSystemMonitor::Update()
{
    if(!isInitialized) return false;
    bool result = false;
    if (_threadHandle.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
        result = _threadHandle.get();
        _threadHandle = std::async(std::launch::async, [] () {
                // Use sleep_for to wait specified time (or sleep_until).
                std::this_thread::sleep_for( std::chrono::milliseconds{sweepTime});
                // Do whatever you want.
                return Inotify::Update();
            } );
    }
    return result;
}
