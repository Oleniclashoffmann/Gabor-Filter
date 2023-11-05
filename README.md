# Gabor FIlter
# Gabor Filter Algorithm for Needle Detection

## Overview
This repository is dedicated to the implementation of a needle detection algorithm using Gabor filtering and Otsu thresholding techniques. The project is designed to process images to enhance, detect, and estimate the axis of needle-like structures. This is particularly useful in medical imaging where precise detection is crucial.

## Project Structure
The codebase is organized with several key classes, each handling a specific aspect of the algorithm:
- `GaborFilter`: The main file of the program.
- `MLESAC`: Implements the MLESAC algorithm for robust axis estimation and needle tip estimation.
- `Otsu_thresh`: Applies Otsu's thresholding method as a preprocessing step.
- `Preprocessing`: Contains additional preprocessing routines to prepare images for analysis (Gabor filter).
- `save_file`: Provides functionality to save processed images or output data.

### Prerequisites
To ensure proper execution of the algorithm, you'll need the following:
- install OpenCV to run algorithm
- change the image repository in main file
- Code is commented for better understanding



