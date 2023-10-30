#pragma once
#include "opencv2/opencv.hpp"
#include <iostream>
#include <ctime>
class MLESAC
{
public: 
	void initialize(cv::Mat src); 
	cv::Point2f algorithm();
	float ret_intercept() { return m_y_intercept; }; 
	float ret_slope() { return m_slope; }; 


private: 
	cv::Mat m_src; 
	cv::Mat original; 
	void create_data(); 
	cv::Point2f needleest();  
	std::vector<cv::Point2i> data;

	cv::Point2f punkt_1;
	cv::Point2f punkt_2;

	cv::Point2f tip_estimation; 

	float m_y_intercept; 
	float m_slope; 



};

