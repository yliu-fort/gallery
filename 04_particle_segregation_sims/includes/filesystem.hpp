#pragma once
#include <vector>
#include <string>
#include <map>

namespace util
{

class IO
{
    using Dict = std::vector<std::map<std::string, std::string>>;
public:
    static void set_config_location(const char* val);
    static void set_output_location(const char* val);
    static void set_output_option_append(bool flag = false) { to_append = flag; }

    // return file handle for input
    std::map<std::string, std::string> query() { return record; }

    // read input "test_input.txt"
    bool scan_config(const char* config_file = nullptr);


    // return file handle for output (.lock file)
    // isfinite: check nan value, if present, refuse autosave and goto error handling
    template <typename T>
    void dump(const std::vector<T>& data, int stride = 1, const char* alt_fname = nullptr);

private:
    static void load_dict();
    static Dict dict;
    static bool loaded;
    static std::string db_loc;
    static std::string list_loc;
    static bool to_append;

    std::map<std::string, std::string> record;

};

}

