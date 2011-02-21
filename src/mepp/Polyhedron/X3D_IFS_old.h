#ifndef IFSDATA_H_
#define IFSDATA_H_

#include <vector>
#include <string>
#include <fstream>

#define NO_COLOR 0
#define COLOR_PER_VERTEX 1
#define COLOR_PER_FACE 2

class IFSData {

public:

	std::vector <std::vector <int> > face;
	std::vector <std::vector <float> > vertex;
	std::vector <std::vector <float> > color;
	std::vector <std::vector <int> > info;

	int colorMode;

	std::string faceFileName;
	std::string vertexFileName;
	std::string colorFileName;
	std::string infoFileName;

	std::fstream faceFile;
	std::fstream vertexFile;
	std::fstream colorFile;
	std::fstream infoFile;


	bool openFiles(std::string xmlFileName) {

		faceFileName.assign(xmlFileName + ".face");
		vertexFileName.assign(xmlFileName + ".vertex");
		colorFileName.assign(xmlFileName + ".color");
		infoFileName.assign(xmlFileName + ".info");

		faceFile.open(faceFileName.data(), std::ios::in | std::ios::out
				| std::ios::trunc);
		vertexFile.open(vertexFileName.data(), std::ios::in | std::ios::out
				| std::ios::trunc);
		colorFile.open(colorFileName.data(), std::ios::in | std::ios::out
				| std::ios::trunc);
		infoFile.open(infoFileName.data(), std::ios::in | std::ios::out
				| std::ios::trunc);

		if (!faceFile || !vertexFile || !colorFile || !infoFile) {
			std::cerr << "Failed to open data files."<< std::endl;
			closeFiles();
			deleteFiles();
			return (false);
		}

		return (true);

	}

	void seek_zero() {
		vertexFile.seekg(0);
		faceFile.seekg(0);
		colorFile.seekg(0);
		infoFile.seekg(0);
	}

	void closeFiles() {
		faceFile.close();
		vertexFile.close();
		colorFile.close();
		infoFile.close();
	}

	void deleteFiles() {
		remove(faceFileName.data());
		remove(vertexFileName.data());
		remove(colorFileName.data());
		remove(infoFileName.data());
	}

	void write_last_face() {
		faceFile << face.back().size();
		std::vector<int>::iterator i;
		for (i = face.back().begin(); i < face.back().end(); i++)
			faceFile << " "<< *i;
		faceFile << std::endl;
	}

	void write_last_vertex() {
		std::vector<float>::iterator i;
		for (i = vertex.back().begin(); i < vertex.back().end(); i++)
			vertexFile << " "<< *i;
		vertexFile << std::endl;
	}

	void write_last_color() {
		std::vector<float>::iterator i;
		for (i = color.back().begin(); i < color.back().end(); i++)
			colorFile << " "<< *i;
		colorFile << std::endl;
	}

};

#endif /*IFSDATA_H_*/
