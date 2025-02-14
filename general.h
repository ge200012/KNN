#ifndef GNNS_GENERAL_H
#define GNNS_GENERAL_H

#include <stdexcept>

//namespace gnns
//{
    class GnnsException : public std::runtime_error
      {
    public:
        GnnsException(const char* message) : std::runtime_error(message) { }
        GnnsException(const std::string& message) : std::runtime_error(message) { }
    };
//}

#endif //GNNS_GENERAL_H
