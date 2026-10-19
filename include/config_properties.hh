#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <unordered_set>
#include <array>

class config_properties
{
protected:
    // Stores the raw key-value pairs as they are successfully loaded
    std::vector<std::pair<std::string, std::string>> loaded_properties;

    // Helper function to strip leading/trailing whitespace
    static std::string trim(const std::string& str)
    {
        auto start = str.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        auto end = str.find_last_not_of(" \t\r\n");
        return str.substr(start, end - start + 1);
    }

    // Derived classes only need to define key parsing and optional validation
    virtual bool parseKey(const std::string& key, std::stringstream& ss) = 0;
    virtual std::vector<std::string> getRequiredKeys() const { return {}; }

public:
    virtual ~config_properties() = default;

    bool read(const std::string& filepath)
    {
        loaded_properties.clear();
        if (filepath.empty()) return false;

        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file: " << filepath << std::endl;
            return false;
        }

        std::unordered_set<std::string> parsed_keys;
        std::string line;

        while (std::getline(file, line))
        {
            std::string trimmedLine = trim(line);

            if (trimmedLine.empty() || trimmedLine[0] == '#') continue;

            size_t equalPos = trimmedLine.find('=');
            if (equalPos == std::string::npos) continue;

            std::string key = trim(trimmedLine.substr(0, equalPos));
            std::string value = trim(trimmedLine.substr(equalPos + 1));

            std::stringstream ss(value);

            // If the derived class successfully parses the key, store the pair for printing
            if (parseKey(key, ss))
            {
                parsed_keys.insert(key);
                loaded_properties.emplace_back(key, value);
            }
            else
            {
                std::cerr << "Warning: Unrecognized property '" << key 
                          << "' in file: " << filepath << "\n";
            }
        }

        file.close();

        // Validate critical keys
        for (const auto& req_key : getRequiredKeys())
        {
            if (parsed_keys.find(req_key) == parsed_keys.end())
            {
                std::cerr << "Warning: Missing critical property '" << req_key 
                          << "' in file: " << filepath 
                          << ". Defaulting to zero.\n";
            }
        }

        return true;
    }

    // Defined ONCE in the base class for ALL derived classes
    void print(std::ostream& os = std::cout) const
    {
        os << "--- Loaded Configuration Properties ---\n";
        for (const auto& [key, value] : loaded_properties)  os << key << " = " << value << "\n";
    }
};