#include "filesystem.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>

using namespace std;
namespace fs =  std::experimental::filesystem;


namespace util
{

// todo: use static member might cause issue
// on feeding multiple configuration files
std::vector< std::map< std::string, std::string > > IO::dict;
bool IO::loaded = false;
bool IO::to_append = false;

std::string IO::db_loc = "./db/";
std::string IO::list_loc = "./resources/test_input.txt";


void IO::set_config_location(const char* val){ list_loc = val; };
void IO::set_output_location(const char* val){ db_loc = val; };

// read input "test_input.txt"
// add one uncompleted case into dict
bool IO::scan_config(const char* config_file)
{
    if(!fs::exists(db_loc))
    {
        fs::create_directory(db_loc);
    }

    // if config_file passed in && config file is not the current one
    // config file exists && regular file
    if(config_file && (fs::path(config_file) != fs::path(list_loc)))
    {
        fs::directory_entry entry(config_file);
        if((fs::exists(entry) && fs::is_regular_file(config_file)))
        {
            loaded = false;
            set_config_location(config_file);
        }
    }

    if(!loaded)
    {
        load_dict();
        loaded = true;
    }

    // pull a record
    if(dict.empty())
    {
        std::cout << "IO::dictionary is empty, return 0...\n ";
        return false;
    }
    record = dict.back();
    dict.pop_back();
    return true;

}

void IO::load_dict()
{
    // load dictionary
    string path = list_loc;
    ifstream infile(path);

    // Scan headers
    auto entrys = vector<string>();
    {
        string line;
        getline( infile, line );
        std::istringstream iss(line);
        string entry;
        while(iss >> entry)
        {
            //cout << entry << endl;
            entrys.push_back(entry);
        }
    }

    // scan records
    dict = Dict();
    for ( string record; getline( infile, record ); )
    {
        map<string, string> table;
        std::istringstream iss(record);
        string val;
        for(const auto i: entrys)
        {
            iss >> val;

            // bind to header
            table.insert(pair<string, string>(i, val));
        }

        // check if output file exists
        fs::directory_entry fid( fs::path(db_loc) /fs::path(table.at("outputFolder")) );
        if(fs::exists(fid) && fs::is_regular_file(fid))
        {
            std::cout << "IO::skip existing record: " << table.at("outputFolder") << "\n";
            continue;
        }

        // add record to local dict
        dict.push_back(table);
    }
}

// return file handle for output (.lock file)
// isfinite: check nan value, if present, refuse autosave and goto error handling
template <typename T>
void IO::dump(const vector<T>& data, int stride, const char* alt_fname)
{
    // Verbose
    std::cout << "Dumping output file...";


    // Compose file path
    fs::path filename = record.at("outputFolder");

    // if specfied filename, use it
    if(alt_fname)
        filename.replace_extension( fs::path(alt_fname) );

    fs::path path = db_loc + filename.string();
    fs::path tpath = path.string() + ".lock"; // prevent file damage

    // if append to file, no .lock file to create
    if(to_append)
        tpath = path.string();

    // Call ostream
    std::ofstream vtkstream;

    // open filestream
    if(to_append)
        vtkstream.open (tpath.string(), std::ios::out | std::ios::app); // std::ios::binary
    else
        vtkstream.open (tpath.string(), std::ios::out | std::ios::trunc); // std::ios::binary

    if (!vtkstream) std::cout<<"ERROR"<<std::endl;

    // Dump field
    int count = 1;
    for(auto& val: data)
    {
        if(!isfinite(val)) { std::cerr << "File dumping err: detected nan/inf values.\n"; }
        vtkstream << val;

        if(count++%stride == 0)
            vtkstream << "\n";
        else
            vtkstream << "\t";
    }

    // Test if ostream is still good
    if(vtkstream.bad() || vtkstream.fail()) { goto bad; }

    // Close filestream
    vtkstream.close();

    // Rename to overwrite old file
    fs::rename(tpath, path);


    cout << "IO::write output to: " << path << endl;

    return;

bad:
    vtkstream.close();
    cerr << "autosave interrupted because ostream failed." << endl;
    return;

}

#ifdef VERBOSE
#undef VERBOSE
#endif

// explicit initialization
template void IO::dump<float>(const vector<float>&, int, const char* );
template void IO::dump<double>(const vector<double>&, int, const char* );
}
