#include "Otsus_thresh.h"

Otsus_thresh::Otsus_thresh(const cv::Mat src)
{
	//save image to memeber variable
	m_src = src; 
}

Otsus_thresh::~Otsus_thresh()
{

}

int Otsus_thresh::threshold()
{
	//total number of pixels
	int total_num_pixels = m_src.cols * m_src.rows;

	//calculate probability of each pixel intensity
	probabilities.resize(256, 0.0);
	for (int i = 0; i < histSize; i++)
	{
		probabilities[i] = hist.at<float>(i) / total_num_pixels;
		
	}

	// Calculate cumulative sums
	std::vector<double> cumulativeSum(histSize);
	cumulativeSum[0] = probabilities[0];
	for (int i = 1; i < histSize; i++) {
		cumulativeSum[i] = cumulativeSum[i - 1] + probabilities[i];
	}

	//Calculate cumulative mean
	std::vector<double> cumulativeMean(histSize);
	cumulativeMean[0] = 0; 
	for (int j = 1; j < histSize; j++)
	{
		cumulativeMean[j] = cumulativeMean[j - 1] + (j * probabilities[j]);
	}

	//calculate global mean
	double globalMean = cumulativeMean[histSize - 1]; 

	//Calculate between class variance and find maximum variance 
	double maxVariance = 0;
	int threshold = 0;
	for (int i = 0; i < histSize; i++) {
		double class0Prob = cumulativeSum[i];
		double class1Prob = 1.0 - class0Prob;

		if (class0Prob == 0 || class1Prob == 0) {
			continue;  // Skip if one of the classes is empty
		}

		double class0Mean = cumulativeMean[i] / class0Prob;
		double class1Mean = (globalMean - cumulativeMean[i]) / class1Prob;

		double variance = class0Prob * class1Prob * std::pow(class0Mean - class1Mean, 2);
		if (variance > maxVariance) {
			maxVariance = variance;
			threshold = i;
		}
	}
	return threshold; 
}

//main function in this class that is called in Preprocessing 
int Otsus_thresh::Hist()
{
	cv::Mat source; 

	source = m_src.clone(); 
	float range[] = { 0, 256 }; //range of intensity levels in integer image

	//calculation of histogram of the image
	bool uniform = true; 
	bool accumulate = false; 
	const float* histRange = { range }; 
	std::vector<cv::Mat> planes;
	cv::split(source, planes); 
	cv::calcHist(&planes[0], 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, uniform, accumulate); 

	//Otsu threshold function 
	int thresh = 0;
	thresh = this->threshold();
	
	//Otsu tuning function
	thresh=this->tuning(thresh);

	//draw the histogram
	//this->drawHistogram(hist);
	
	return thresh;
}

void Otsus_thresh::drawHistogram(cv::Mat& hist_1)
{
	int hist_w = m_src.cols;
	int hist_h = m_src.rows; 
    
    int bin_w = cvRound((double)hist_w / histSize);

    cv::Mat histImage(hist_h, hist_w, CV_8UC3, cv::Scalar(0, 0, 0));

    cv::normalize(hist, hist, 0, histImage.rows, cv::NORM_MINMAX, -1, cv::Mat());


    for (int i = 1; i < histSize; i++) {
        cv::line(
            histImage,
            cv::Point(bin_w * (i - 1), hist_h - cvRound(hist_1.at<float>(i - 1))),
            cv::Point(bin_w * (i), hist_h - cvRound(hist_1.at<float>(i))),
            cv::Scalar(255, 0, 0), 2, 8, 0);
    }
}

float Otsus_thresh::tuning(int otsu)
{
	//Entropy calculation
	double entropy= 0.0;

	for(int i = 0; i < histSize - 1; i++)
	{
		if (probabilities[i] > 0.0) {
			entropy = entropy + (probabilities[i] * log2(probabilities[i]));
		}
	}
	entropy = -1 * entropy; 
	 
	//Foreground entropy H_b
	float P_b = 0; 
	for (int i = 0; i <= otsu; i++)
	{
		P_b = probabilities[i] + P_b;
		 
	}
	std::vector<double> input1;
	input1 = probabilities; 

	float H_b = 0;
	for (int i = 0; i <= otsu; i++)
	{
		if (input1[i] > 0.0) {
			H_b = H_b + (input1[i] / P_b) * std::log((input1[i]) / P_b);
		}
	}
	H_b = -1*H_b;
 

	//Background Entropy H_f
	float P_f = 0; 
	float H_f = 0;
	for (int i = otsu+1; i <= histSize-1; i++)
	{
		P_f = probabilities[i] + P_f;
	}
	std::vector<double> input2;
	input2 = probabilities; 

	for (int k = otsu + 1; k <= histSize - 1; k++)
	{
		if (probabilities[k] > 0.0) {
			H_f = H_f + (probabilities[k] / P_f) * std::log((probabilities[k]) / P_f);
		}
	}
	H_f =-1* H_f ; 


	//Calculate W_ratio
	int W_ratio = 0;
	float W_ratio_final=0;
	uchar pixelvalue; 
	int pixel_val; 

	//W ration I(i, j) nur 1 oder Null addieren oder Intensitätswert ?
	for (int i = 0; i < m_src.rows; i++)
	{
		for (int k = 0; k < m_src.cols; k++)
		{
			pixelvalue = m_src.at<uchar>(i, k);
			pixel_val = static_cast<int>(pixelvalue); 

			if (pixel_val >= otsu)
			{
				W_ratio = W_ratio + 1;
			}
		}
	}
	W_ratio_final = static_cast<double>(W_ratio) / (m_src.cols * m_src.rows);


	//calculate final alpha tuning value 
	float alpha; 
	float T_tuned; 
	alpha =(entropy + (H_b / H_f) - W_ratio_final); 

	T_tuned = otsu * (alpha +1)/2; 
	T_tuned = T_tuned/2.2; 

	//in liver tuned /2.2
	//in gel 2.7
	//this value is the mentioned manual adjustment value in the thesis, since the calculated T_tuned is to big

	return T_tuned; 

}
