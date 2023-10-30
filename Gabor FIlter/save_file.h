#pragma once
#include<string>
#include "opencv2/opencv.hpp"
#include <fstream>

class save_file
{
	public:
		save_file(std::string path, cv::Point2f tip, std::chrono::duration<double> time, double angle);
private: 

};

