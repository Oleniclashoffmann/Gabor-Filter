# Gabor-Filter Needle Detection

> **Bachelor's Thesis Project** — 2 of 3 Real-time needle axis and tip estimation in ultrasound images using a ray-based statistical filtering algorithm.

---

## Problem Statement

Accurate needle tip localization during ultrasound-guided medical procedures (e.g. biopsies, regional anesthesia) is critical for patient safety. However, needles are often poorly visible in ultrasound images due to noise, speckle artifacts, and tissue interference. Manual tracking by clinicians is error-prone and limits procedural efficiency. This project addresses the challenge by implementing a **fully automated, multi-stage image processing pipeline** that enhances needle visibility, detects the needle trajectory, and precisely estimates the tip position — all in real time.

---

## Algorithm Overview

The detection pipeline combines classical computer-vision techniques in a multi-stage approach:

1. **Gabor Filtering** — A frequency- and orientation-selective filter enhances elongated, line-like structures (the needle) while suppressing irrelevant texture.
2. **Entropy-Tuned Otsu Thresholding** — A custom implementation of Otsu's method determines an optimal binarisation threshold, further refined by an entropy-based tuning factor that adapts to the foreground/background intensity distribution.
3. **Morphological Cleanup & Connected-Component Analysis** — Closing and opening operations remove small artefacts; connected-component labelling discards clusters below a minimum area.
4. **MLESAC Axis Estimation** — A Maximum-Likelihood variant of RANSAC iteratively fits a line model to the remaining white pixels, using Expectation-Maximisation (EM) to separate inliers from outliers and minimise a probabilistic cost function.
5. **Intensity-Gradient Tip Localisation** — The algorithm walks along the estimated axis and measures the mean intensity change in a local neighbourhood to detect the point of maximum gradient — the needle tip.

---

## Pipeline Steps

```
Input Image
    │
    ▼
┌──────────────────────┐
│  1. Rotate Image     │  (configurable angle correction)
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  2. Convert to       │  BGR → Grayscale → Float32
│     Grayscale        │
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  3. Gabor Filter     │  Kernel: 51×51, σ=4, θ=π/2, λ=13, γ=0.7
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  4. Median Filter    │  3×3 kernel — removes salt-and-pepper noise
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  5. Otsu Threshold   │  Custom implementation + entropy-based tuning
│     (Entropy-Tuned)  │
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  6. Morphological    │  Close → Open (3×3 rectangular element)
│     Operations       │
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  7. Connected-       │  Remove clusters < 300 px
│     Component Filter │
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  8. MLESAC           │  100 iterations, EM-based γ estimation
│     Axis Estimation  │
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  9. Needle Tip       │  Intensity-gradient analysis along axis
│     Localisation     │
└──────────┬───────────┘
           ▼
  Output: Axis line + Tip point + Angle + Computation time
```

---

## Tech Stack

| Layer | Technology |
|---|---|
| Language | C++17 |
| Computer Vision | OpenCV 4.x |
| Build System | Microsoft Visual Studio (`.vcxproj` / `.sln`) |
| Platform | Windows (readily portable to macOS / Linux with CMake) |

---

## Project Structure

```
Gabor-Filter/
├── Gabor FIlter.sln                 # Visual Studio solution file
├── README.md
└── Gabor FIlter/
    ├── Gabor FIlter.cpp             # Main entry point & processing loop
    ├── Preprocessing.h / .cpp       # Gabor filter, median filter, morphology, CC analysis
    ├── Otsus_thresh.h / .cpp        # Custom Otsu threshold with entropy-based tuning
    ├── MLESAC.h / .cpp              # MLESAC axis estimation & tip localisation
    ├── save_file.h / .cpp           # Result persistence (tip, angle, timing → .txt)
    ├── Gabor FIlter.vcxproj         # VS project file
    └── Gabor FIlter.vcxproj.filters # VS project filters
```

---

## Component Responsibilities

### `Gabor FIlter.cpp` — Main Loop
- Reads images sequentially from a configurable directory.
- Orchestrates the full pipeline: Preprocessing → MLESAC → visualisation.
- Measures and prints per-frame computation time and estimated needle angle.
- Draws the detected axis line and tip circle on the original (rotated) image.

### `Preprocessing` (`.h` / `.cpp`)
- **Image rotation** to a configurable angle.
- **Gabor filtering** with a 51×51 kernel (σ = 4, θ = π/2, λ = 13, γ = 0.7).
- **Median blur** (3×3) for noise reduction.
- **Otsu binarisation** via the `Otsus_thresh` class, with a safety floor of 110.
- **Morphological close + open** (3×3 rectangular structuring element).
- **Connected-component analysis** to discard small clusters (< 300 px).

### `Otsus_thresh` (`.h` / `.cpp`)
- Computes a 256-bin histogram and per-bin probabilities.
- Implements **Otsu's method** from scratch: cumulative sums, cumulative means, between-class variance maximisation.
- **Entropy-based tuning**: calculates image entropy, foreground/background entropies (H_b, H_f), and a pixel-ratio weight (W_ratio) to derive an adaptive correction factor α that shifts the threshold for improved segmentation.

### `MLESAC` (`.h` / `.cpp`)
- **Data extraction**: collects all non-zero pixel coordinates from the binary image.
- **MLESAC loop** (100 iterations, 2 random points per iteration):
  - Fits a line model via `cv::fitLine`.
  - Computes point-to-line distances for all data points.
  - Runs **EM** (up to 10 inner iterations) to estimate the inlier mixing coefficient γ.
  - Evaluates a negative-log-likelihood cost and retains the best model.
- **Needle tip estimation** (`needleest`): iterates along the estimated axis using `cv::LineIterator`, compares intensity sums in a ±20 px band between adjacent positions, and returns the point with the maximum mean intensity drop — the likely needle tip.

### `save_file` (`.h` / `.cpp`)
- Appends per-frame results (image name, tip coordinates, computation time, estimated angle) to text files for offline error analysis.

---

## Prerequisites

| Requirement | Details |
|---|---|
| **Compiler** | MSVC 2019+ (or any C++17-compatible compiler) |
| **OpenCV** | Version 4.x — core, imgproc, imgcodecs, highgui modules |
| **IDE** | Visual Studio 2019 / 2022 (solution file included) |

---

## Getting Started

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/Gabor-Filter.git

# 2. Open the solution in Visual Studio
#    File → Open → Project/Solution → "Gabor FIlter.sln"

# 3. Configure OpenCV
#    • Set OpenCV include path in Project → Properties → C/C++ → Additional Include Directories
#    • Set OpenCV lib path in Linker → Additional Library Directories
#    • Add required .lib files (e.g. opencv_world4xx.lib) to Linker → Input

# 4. Set the image directory
#    In "Gabor FIlter.cpp", update the `path` variable (line 26) to point to your image folder.

# 5. Build & Run
#    Press F5 (Debug) or Ctrl+F5 (Release) in Visual Studio.
```

> **Tip:** To port to macOS / Linux, create a `CMakeLists.txt` and link against the same OpenCV modules.

---

## Configuration Parameters

The following parameters can be tuned directly in the source code:

| Parameter | File | Default | Description |
|---|---|---|---|
| `angle` | `Gabor FIlter.cpp` | `90` | Global rotation applied to every input image (degrees) |
| `ksize` | `Preprocessing.cpp` | `51` | Gabor kernel size (N×N) |
| `sigma` | `Preprocessing.cpp` | `4.0` | Gaussian envelope standard deviation |
| `theta` | `Preprocessing.cpp` | `π/2` | Orientation of the Gabor stripes |
| `lambd` | `Preprocessing.cpp` | `13.0` | Wavelength of the sinusoidal component |
| `gamma` | `Preprocessing.cpp` | `0.7` | Spatial aspect ratio |
| `minClusterSize` | `Preprocessing.cpp` | `300` | Minimum connected-component area (px) |
| `n` (iterations) | `MLESAC.cpp` | `100` | Number of MLESAC iterations |
| `sigma` (MLESAC) | `MLESAC.cpp` | `0.01` | Inlier noise standard deviation |
| `gamma` (EM init) | `MLESAC.cpp` | `0.5` | Initial inlier mixing coefficient |
| `thresh` floor | `Preprocessing.cpp` | `110` | Manual safety floor for Otsu threshold |
| Tuning divisor | `Otsus_thresh.cpp` | `2.2` | Manual scaling of entropy-tuned threshold (tissue-dependent) |

---

## Context

This project was developed as the practical component of a **Bachelor's thesis** focused on **needle detection in medical ultrasound images**. The dataset consists of ultrasound frames of liver tissue with inserted biopsy needles at various insertion angles. The goal was to build a lightweight, real-time-capable detection pipeline that could eventually be integrated into an ultrasound guidance system. The algorithm was benchmarked on 100 sequential frames per dataset, measuring:

- **Needle angle accuracy** — estimated vs. ground-truth insertion angle.
- **Tip localisation error** — Euclidean distance between estimated and actual tip position.
- **Computational time** — per-frame processing duration.

---

## References

- Gabor, D. (1946). *Theory of communication.* Journal of the IEE, 93(26), 429–457.
- Otsu, N. (1979). *A Threshold Selection Method from Gray-Level Histograms.* IEEE Transactions on Systems, Man, and Cybernetics, 9(1), 62–66.
- Torr, P. H. S. & Zisserman, A. (2000). *MLESAC: A New Robust Estimator with Application to Estimating Image Geometry.* Computer Vision and Image Understanding, 78(1), 138–156.
- Kapur, J. N., Sahoo, P. K. & Wong, A. K. C. (1985). *A new method for gray-level picture thresholding using the entropy of the histogram.* Computer Vision, Graphics, and Image Processing, 29(3), 273–285.

---

## License

This project is provided for **educational and portfolio purposes**. No explicit license is applied. If you wish to use or adapt this code, please contact the author.
