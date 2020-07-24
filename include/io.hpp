/*
  Very simple I/O
 */
#include <fstream>
#include <string>

namespace IO{
	template<typename T>
	void write(std::string filename, std::vector<T> vec){
		std::ofstream file;
		file.open(filename);
		for(int i=0; i<vec.size(); ++i){
			file << vec[i] << "\n";
		}
		file.close();
	};
}
