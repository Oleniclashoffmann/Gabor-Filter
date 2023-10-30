#include "MLESAC.h"

void MLESAC::initialize(cv::Mat src)
{
	//save the image
	m_src = src.clone();
	original = src.clone(); 
	cv::cvtColor(original, original, cv::COLOR_GRAY2BGR); 
	m_src.convertTo(m_src, CV_32F);
}

cv::Point2f MLESAC::algorithm()
{
	//initialize variables of the MLESAC algorithm
	int d = 2; // Number of random data points to be selected in each iteration
	float sigma = 0.01; //standard deviation of the inlier noise (algorithm assumes that inliers follow Gaussian distribution with standard deviation of sigma) 
	int n = 100; //iteration number 
	cv::Vec4f needleModel; // Will store the final model parameters for the needle



	//create dataset 
	this->create_data(); 
	if (data.size() != 0)
	{
		if (data.size() != 0)
		{
			std::cout << "Data not empty! " << std::endl;
		}
		//Seed the random number generator
		srand(static_cast<unsigned int> (time(0)));

		//MLESAC algorithm for needle location
		for (int i = 0; i < n; i++)
		{
			//Randomly select two pixels as MMS from data 
			std::vector<cv::Point2f> mms;
			for (int j = 0; j < d; j++)
			{
				int x_1 = rand() % data.size();
				int x_2 = rand() % data.size();
				mms.push_back(cv::Point2f(data[x_1].x, data[x_1].y));
				mms.push_back(cv::Point2f(data[x_2].x, data[x_2].y));
			}

			//Use the two pixels in MMS to estimate the needle axis
			cv::fitLine(mms, needleModel, cv::DIST_L2, 0, 0.01, 0.01);

			//Calculate the error for each pixel in image
			std::vector<float> errors;
			for (int y = 0; y < m_src.rows; y++)
			{
				for (int x = 0; x < m_src.cols; x++)
				{
					if (m_src.at<float>(y, x) != 0)
					{
						float A = needleModel[1];
						float B = -needleModel[0];
						float C = needleModel[0] * needleModel[3] - needleModel[1] * needleModel[2];
						float distance = std::abs(A * x + B * y + C) / std::sqrt(A * A + B * B);
						errors.push_back(distance);
					}
				}
			}

			// Estimate the mixing coefficient gamma using the EM (Expectation-Maximization) approach
			//liver 0.5
			float gamma = 0.5;
			for (int iter = 0; iter < 10; iter++)
			{
				// Compute inlier and outlier weights
				float sumInliers = 0;
				float sumOutliers = 0;
				int numInliers = 0;
				int numOutliers = 0;

				for (int j = 0; j < errors.size(); j++)
				{
					float inlierLikelihood = exp(-errors[j] * errors[j] / (2 * sigma * sigma));
					float outlierLikelihood = 1 - gamma;
					float weight = inlierLikelihood / (inlierLikelihood + outlierLikelihood);

					if (weight > 0.5)
					{
						sumInliers += errors[j];
						numInliers++;
					}
					else
					{
						sumOutliers += errors[j];
						numOutliers++;
					}
				}

				// Update gamma based on inlier and outlier counts
				float newGamma = static_cast<float>(numInliers) / (numInliers + numOutliers);
				if (std::abs(newGamma - gamma) < 0.001)
					break;

				gamma = newGamma;
			}

			// Compute the cost function and update the model if a better one is found
			float minCost = std::numeric_limits<float>::max();
			for (int j = 0; j < errors.size(); j++)
			{
				float inlierLikelihood = exp(-errors[j] * errors[j] / (2 * sigma * sigma));
				float outlierLikelihood = 1 - gamma;
				float cost = -log(gamma * inlierLikelihood + outlierLikelihood);

				
				if (cost < minCost)
				{
					minCost = cost;

					//if angle greater than 90° than needleModel[3] + errors[j] * needleModel[0]
					//if angle smaller than 90° than needleModel[3] - errors[j] * needleModel[0]

					// Update the model based on the cost
					needleModel = cv::Vec4f(needleModel[0], needleModel[1], needleModel[2] + errors[j], needleModel[3] - errors[j] * needleModel[0]);
				}

			}

		}


		// Visualization phase: Draw the needle line and display the result
		float slope = needleModel[1] / needleModel[0];
		float yIntercept = needleModel[3] / needleModel[0];
		m_slope = slope;
		m_y_intercept = yIntercept; 
		
		cv::Point2f p1(0, yIntercept);
		cv::Point2f p2(m_src.cols, slope * m_src.cols + yIntercept);
		cv::Mat result = m_src.clone();

		cv::cvtColor(result, result, cv::COLOR_GRAY2BGR);
		cv::line(result, p1, p2, cv::Scalar(0, 255, 0), 2);
		cv::imshow("MLESAC", result);

		punkt_1 = p1; 
		punkt_2 = p2; 

		//calculate the angle between the x-axis and the fitted line
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;

		// Calculate the angle in radians
		double angle_radians = atan2(dy, dx);

		// Convert the angle to degrees if needed
		double angle_degrees = angle_radians * 180.0 / CV_PI;
		 
		//estimate needle tip and return it
		cv::Point2f tip = this->needleest(); 

		return tip; 
	}
}

void MLESAC::create_data()
{
	//iterate over all pixels in the image and save x and y coordinates of the pixels that are white pixels
	for (int x = 0; x < m_src.cols; x++)
	{
		for (int y = 0; y < m_src.rows; y++)
		{
			if (m_src.at<float>(y, x) > 0)
			{
				cv::Point2i new_point(x, y);

				data.push_back(new_point);
			}
		}
	}
}

cv::Point2f MLESAC::needleest()
{
	//Initializing variables for intensity change
	float maxchange = 0.0f; 
	cv::Point maxChangePoint;

	//Create line iterator to iterate alongside the calculated needle axis
	cv::LineIterator it(punkt_1, punkt_2, 8); 
	std::vector<int> intensity;
	std::vector<int> intensity_aft; 

	//iterate alongside the needle axis
	for (int i = 1; i < it.count; ++i, ++it)
	{
		cv::Point currentPoint = it.pos(); //Current point in iteration 
		cv::Point copyPoint = it.pos(); //dublicate of current point

		// Collecting intensity values in a range around the current point to optimize the needle tip estimation 
		for (int y = -20; y < 20; y++)
		{
			copyPoint.y = currentPoint.y + y;
			if (copyPoint.x >= 0 && copyPoint.y >= 0 && copyPoint.x < m_src.cols && copyPoint.y < m_src.rows)
			{
				intensity.push_back(m_src.at<float>(copyPoint));
			}
		}

		// Setting up a second point (4 pixels to the right) for intensity comparison
		copyPoint = currentPoint; 
		copyPoint.x = copyPoint.x + 4; 

		// Collecting intensity values in a range around the second point
		for (int y_2 = -20; y_2 < 20; y_2++)
		{
			copyPoint.y = currentPoint.y + y_2;
			if (copyPoint.x >= 0 && copyPoint.y >= 0 && copyPoint.x < m_src.cols && copyPoint.y < m_src.rows)
			{
				
				intensity_aft.push_back(m_src.at<float>(copyPoint));
			}
		}

		std::vector <int> diff; // To store differences in intensity between the two points

		// Ensuring both intensity vectors are of the same size before processing
		if (intensity.size() == intensity_aft.size())
		{
			int sum = 0;
			for (int i = 1; i < intensity.size(); i++)
			{
				// Calculating the difference in intensity between the current and the second point
				int difference = intensity[i] - intensity_aft[i];
				diff.push_back(difference);
				sum += difference;
			}

			// Calculating the mean difference in intensity
			double mean_difference = static_cast<double>(sum) / diff.size();

			// Updating the maximum change in intensity and its corresponding point if the current mean difference is greater
			if (mean_difference > maxchange)
			{
				maxchange = mean_difference;
				maxChangePoint = currentPoint;
			}

			// Clearing the intensity vectors for the next iteration
			intensity.clear();
			intensity_aft.clear();
		}
	}

	// Returning the point with the maximum change in intensity, likely representing the needle tip
	return maxChangePoint; 
}