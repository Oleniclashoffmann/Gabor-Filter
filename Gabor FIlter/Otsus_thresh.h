#pragma once
#include "opencv2/opencv.hpp"
#include <iostream>
#include <cmath>


class Otsus_thresh
{
public:
	Otsus_thresh(const cv::Mat src);
	~Otsus_thresh();
 
	int Hist(); 
	float tuning(int otsu);
	void drawHistogram(cv::Mat& hist); 

private:
	cv::Mat m_src, hist; 
	int histSize = 256; 
	std::vector<double> probabilities;

	int threshold();
	

};

