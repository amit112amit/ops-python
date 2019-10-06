#include "PolyDataReader.h"

namespace OPS
{

std::string &ltrim(std::string &str, const std::string &chars = "\t\n\v\f\r ")
{
	str.erase(0, str.find_first_not_of(chars));
	return str;
}

std::string &rtrim(std::string &str, const std::string &chars = "\t\n\v\f\r ")
{
	str.erase(str.find_last_not_of(chars) + 1);
	return str;
}

std::string &trim(std::string &str, const std::string &chars = "\t\n\v\f\r ")
{
	return ltrim(rtrim(str, chars), chars);
}

int read_polydata(std::string filename, std::vector<std::array<double, 3>> &coords,
				  std::vector<std::vector<int>> &cells)
{
	std::ifstream file(filename.c_str());
	if (!file.is_open())
	{
		std::cerr << "Failed to open file " << filename << "!\n";
		return -1;
	}

	const std::string coordsmarker = "POINTS";
	const std::string cellsmarker = "POLYGONS";

	coords.clear();
	cells.clear();

	// Main reading loop:
	//
	std::string line, ignorestring, numstring;
	int numpoints, numcells, count;
	double number;
	std::vector<double> flatcoords;
	while (std::getline(file, line))
	{
		if (line.find(coordsmarker) != std::string::npos)
		{
			line = trim(line);
			std::istringstream ss(line);
			// Ignore the `coordsmarker`
			ss >> ignorestring;
			// Read the number of POINTS
			ss >> numpoints;
			count = 3 * numpoints;
			// Keep reading lines till we get `count` floating point numbers
			while (count > 0 && (bool)std::getline(file, line))
			{
				line = trim(line);
				std::istringstream ss(line);
				std::array<double, 3> point;
				while (!ss.eof())
				{
					ss >> number;
					count--;
					flatcoords.push_back(number);
				}
			}
		}
		// Read the POLYGONS
		else if (line.find(cellsmarker) != std::string::npos)
		{
			line = trim(line);
			std::istringstream ss(line);
			// Ignore the `cellsmarker`
			ss >> ignorestring;
			// Read the number of lines of data cotaining cells
			ss >> numcells;
			// Read the number of floating points to read
			ss >> count;
			// Keep reading lines till we get `count` floating point numbers
			while (count > 0 && (bool)std::getline(file, line))
			{
				line = trim(line);
				std::istringstream ss(line);
				if (!ss.eof())
				{
					int numvertices;
					int vertexid;
					std::vector<int> cell;
					ss >> numvertices;
					count--;
					for (int i = 0; i < numvertices; ++i)
					{
						ss >> vertexid;
						count--;
						cell.push_back(vertexid);
					}
					cells.push_back(cell);
				}
			}
		}
	}

	// Reshape the flatcoords
	std::array<double, 3> point;
	int index = 0;
	for (const auto coord : flatcoords)
	{
		if (index < 3)
		{
			point[index++] = coord;
		}
		else
		{
			coords.push_back(point);
			index = 0;
			point[index++] = coord;
		}
	}
	coords.push_back(point);
	return 1;
}

} // namespace OPS