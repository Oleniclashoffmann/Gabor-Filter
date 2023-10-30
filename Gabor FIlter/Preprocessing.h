#pragma once
#include "opencv2/opencv.hpp"
#include <iostream>
#include "Otsus_thresh.h"

class Preprocessing
{
public: 

	Preprocessing(cv::Mat src, double angle); 
	~Preprocessing(); 

	void GaborFilter(); 
	cv::Mat retimage(); 
private: 

	cv::Mat m_src; 
	cv::Mat m_src_8bit; 
	cv::Mat dest; 

	double m_angle; 

	void connected_components(cv::Mat filteredImage); 
};

