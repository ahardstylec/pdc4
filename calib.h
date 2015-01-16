#ifndef CALIB_H
#define CALIB_H
#include <opencv2/calib3d/calib3d.hpp>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <iostream>

void calibCams(int camL, int camR,int CAPW, int CAPH);

#endif // CALIB_H

