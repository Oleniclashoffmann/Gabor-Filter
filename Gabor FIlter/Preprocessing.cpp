#include "Preprocessing.h"

Preprocessing::Preprocessing(cv::Mat src, double angle)
{
	//both the image matrix and the rotation angle are saved to member variables of the class
	m_src = src.clone(); 
	m_src_8bit = src.clone(); 

	m_angle = angle; 
	

}

Preprocessing::~Preprocessing()
{

}

void Preprocessing::GaborFilter()
{
	
	//rotate Image
	cv::Point2f center(static_cast<float>(m_src.cols) / 2, static_cast<float>(m_src.rows) / 2); 
	cv::Mat rotateMatrix = cv::getRotationMatrix2D(center, m_angle, 1.0); 
	cv::warpAffine(m_src, m_src, rotateMatrix, m_src.size());

	cv::imshow("Original bild", m_src);
	cv::cvtColor(m_src, m_src, cv::COLOR_BGR2GRAY); //convert the color to gray scale image
	m_src.convertTo(m_src, CV_32F); //convert the image to a float image for later calculations

	

	double min, max;
	
	//define values for the Gabor Kernel and Gabor filter 
	int ksize = 51; //is the Kernel size in Bachelor thesis defined as N in section of Gabor filter

	double sigma = 4, theta = CV_PI / 2.0, lambd = 13, gamma = 0.7, psi = 0;
	/*
	- sigma is Standard deviation of gaussian envelope, 
	- theta is orientation of the normal to the parallel stripes of Gabor function, 
	- lambd is Wavelength of the sinusoidal factor, 
	- gamma is Spatial aspect ration, psi is phase offset*/

	cv::Mat dest; 
	cv::Mat kernel; 
	
	//Gabor Filter
	kernel = cv::getGaborKernel(cv::Size(ksize, ksize), sigma, theta, lambd, gamma, psi);
	cv::filter2D(m_src, dest, CV_32F, kernel);
	cv::Mat GaborKernel = dest.clone(); 

	cv::normalize(GaborKernel, GaborKernel, 0, 255, cv::NORM_MINMAX, CV_8U);

	dest = GaborKernel;  


	//Median Filter 
	cv::minMaxLoc(dest, &min, &max);
	dest.convertTo(dest, CV_8U, 255.0 / (max - min), -255.0 * min / (max - min));
	cv::medianBlur(dest, dest, 3);

	//Otsu's thresholding class
	int thresh;
	Otsus_thresh new_object(dest); 
	thresh = new_object.Hist(); 

	//this value is manual threhsold limit
	if (thresh <= 100)
		thresh = 110; 

	//threshold function 
	cv::threshold(dest, dest, thresh, 255, cv::THRESH_BINARY);
	cv::imshow("Otsu-thresh", dest); //shows the thresholded image
	

	//Morphological Structering algorithm to eliminate small objects in the image
	int kernelSize = 3; 
	cv::Mat closing; 
	cv::Mat elementClosing = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(kernelSize, kernelSize)); 
	cv::morphologyEx(dest, closing, cv::MORPH_CLOSE, elementClosing); 

	cv::Mat filteredImage; 
	cv::Mat elementOpening = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(kernelSize, kernelSize));
	cv::morphologyEx(closing, filteredImage, cv::MORPH_OPEN, elementOpening); 

	//Than Connected Components algorihtm is applied
	this->connected_components(filteredImage.clone()); 

	//Normalize float image to show it 
	cv::normalize(m_src, m_src, 0, 1, cv::NORM_MINMAX);

}

cv::Mat Preprocessing::retimage()
{
	return dest;
}

void Preprocessing::connected_components(cv::Mat filteredImage)
{
	//initialize minimum cluster size
	int minClusterSize = 300; 
	
	//perform connected component analysis
	cv::Mat labels, stats, centroids; 
	int imageType = filteredImage.type(); 

	//only works if image type is CV_8U (integer image)
	if (imageType != CV_8U)
		std::cout << "\nImage type ist not 8-bit " << std::endl; 
	else {

		//function to analyze image of the clusters (connected components) in the image
		int numLabels = cv::connectedComponentsWithStats(filteredImage, labels, stats, centroids);

		//Iterate through the connected components 
		for (int label = 1; label < numLabels; label++)
		{
			//get cluster size of iterated label 
			int area = stats.at<int>(label, cv::CC_STAT_AREA);

			//Check if the cluster is small and erase it
			if (area < minClusterSize)
			{
				for (int row = 0; row < filteredImage.rows; row++)
				{
					for (int col = 0; col < filteredImage.cols; col++)
					{
						if (labels.at<int>(row, col) == label)
						{
							filteredImage.at<uchar>(row, col) = 0;
						}
					}
				}
			}
		}
	}

	//update new image and display image
	dest = filteredImage.clone(); 
	cv::imshow("filtered Image in con", dest);
}
