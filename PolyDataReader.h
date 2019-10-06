#if !defined(__POLYDATAREADER_H__)
#define __POLYDATAREADER_H__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace OPS
{

std::string &ltrim(std::string &str, const std::string &chars);
std::string &rtrim(std::string &str, const std::string &chars);
std::string &trim(std::string &str, const std::string &chars);
int read_polydata(std::string filename, std::vector<std::array<double, 3>> &coords,
                  std::vector<std::vector<int>> &cells);

} // namespace OPS
#endif //__POLYDATAREADER_H__