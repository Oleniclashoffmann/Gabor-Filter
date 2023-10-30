#include "opencv2/opencv.hpp"
#include <iostream>
#include "Preprocessing.h"
#include "MLESAC.h"
#include "save_file.h"


double angle = 90; //Rotation angle for the input image (in my case images where all rotated by 90° so 90° is 0°) 

int main()
{
    //Initialize variables
    cv::Mat image, dest; 
    int i = 0; // is the image number in specific folder (e.g. image name: \name3) 
    float est_angle; 
    double algo_angle = 0; 



    while (i < 100) //maximum is set to 100 in my case ince only 100 images neede to be read in 
    {
        auto start = std::chrono::high_resolution_clock::now(); //saves the time value to later calculate the computational time of the algorithm

        //following code is to open folder with the images and read the specific image from the iteration
        std::string path = "C:/Users/ohoff/Documents/Bachelorarbeit/Datensets/Pictures/Liver/liv01_pos01_oth01/";
        path += "name";
        path += std::to_string(i);
        path = path + ".jpg"; 
        image = cv::imread(path); 
      
        //if file empty error message occurs
        if (!image.data)
        {
            std::cout << "Could not open or find the image!" << std::endl; 
        }


        //first step is the preprocessing of the image that was read in
        ///////// PreProcessing /////////
        Preprocessing element1(image, angle);
        element1.GaborFilter(); 
        dest = element1.retimage();

        //if image after preprocessing is 0 error 
        if(dest.size != 0)
            cv::imshow(" Nach Gabor Filter ", dest); 

        //Second step is the MLESAC tracking algorithm and the tip algorithm in one class
        /////////MLESAC/////////
        MLESAC element_mlesac; 
        element_mlesac.initialize(dest); 
        cv::Point2f final_tip = element_mlesac.algorithm(); 

        //this class algorithm returns slope and intercept of the needle axis
        float slope = element_mlesac.ret_slope();
        float intercept = element_mlesac.ret_intercept();

        //this part is to rotate the image to the angle variable defined above in this case 90°
        cv::Point2f center(static_cast<float>(image.cols) / 2, static_cast<float>(image.rows) / 2);
        cv::Mat rotateMatrix = cv::getRotationMatrix2D(center, angle, 1.0);
        cv::warpAffine(image, image, rotateMatrix, image.size());

        //here the start and end point of the axis is defined based on the slope and intercept and the tip estimation
        cv::Point startPoint(0, intercept); 
        cv::Point endPoint(final_tip.x, slope * final_tip.x + intercept); 
        cv::line(image, startPoint, endPoint, cv::Scalar(255, 0, 0), 4); //this draws the line of the needle axis

        //the needle tip is drawn in the following code lines
        cv::circle(image, final_tip, 5, cv::Scalar(0, 0, 255), 2);
        cv::imshow("Circle", image);

        //as mentioned above this part is to calculate the computational time of the algorithm 
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> time = end - start;
        std::cout << "Computational time: " << time.count() << std::endl;

        //this part is to calculate the estimated angle of the algorithm 
        algo_angle = atan(slope); 
        algo_angle = algo_angle * (180.0 / CV_PI); 
        std::cout << "Estimated needle angle: " << algo_angle << std::endl;

        //the following code that is commented out saves all values in a txt file 
        //save_file object_file("name"+std::to_string(i), final_tip, time, algo_angle);

        std::cout << "-----------------------------------------------" << std::endl; 
        cv::waitKey(1);
        i++;
    }
}

