#include "opencv2/highgui/highgui.hpp"
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//feature detect
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include "calib.h"
#include "stopwatch.h"

#include <sstream>
#include <iostream>
#include <cstdio>
#include <queue>
#include <omp.h>
#include <iostream>
#include <cassert>
#include <chrono>

#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <cmath> // for fabs()
#include <iomanip>
#include <sstream>
#include <thread>
#ifdef __APPLE__
#define CAPW 1440/3 //screen resolution x / 3
#define CAPH 810/3  //screen resolution y / 3
#else
#define CAPW 1920/3 //screen resolution x / 3
#define CAPH 1080/3  //screen resolution y / 3
#endif

#define DEBUG true
#define USECALIBDATA true
#define SGBM true

#if DEBUG
#include <chrono>
using namespace std::chrono;
#endif

using namespace cv;
using namespace std;
void blur(Mat& image, int i, int k)
{
	try{
		Mat subImg = image(Range(((image.rows/4)*(i)) , ((image.rows/4)*(i+1))),  Range(((image.cols/4)*(k)), ((image.cols/4)*(k+1)) ));
		medianBlur(subImg,subImg, 97);	
		//bitwise_not(subImg, subImg);	
	}catch(...){};
	
}

int main(int argc, char *argv[]) {

	int threadnum = atoi(argv[1]);
	// const int threadnum = 8;
	int real_thread_num= 0;
	vector<thread>threads;
	chrono::system_clock::time_point startTime = chrono::system_clock::now();
	Mat image = imread(string(argv[2]));
	// cout << image.cols << "x" << image.rows << endl;
	omp_set_num_threads(threadnum);
	int i;
#pragma omp parallel for private(i) lastprivate(real_thread_num)
	for(int i=0; i < threadnum;i++)
	{
		real_thread_num = omp_get_num_threads();
		for(int k=0;k < threadnum;k++){
	//		cout  << "WIDTH" << i << "="<< ((image.cols/4)*(k)) << "-" << ((image.cols/4)*(k+1)) << endl;
	//		cout << "HEIGHT" << i << "="<< ((image.rows/4)*(i)) << "-" << ((image.rows/4)*(i+1)) << endl;
	//		cout << "___________________________________" << endl;
			blur(image, i, k);
			//threads.push_back(thread([&](){blur(image, i, k);}));
		}
	}
//	for (thread& t : threads) {
//		t.join();
//	}
	
	chrono::system_clock::time_point endTime = chrono::system_clock::now();
	chrono::microseconds microRunTime =
			chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	double runTime = microRunTime.count() / 1000000.0;

    cout << setprecision(8) << "time " << runTime << " seconds." << endl
		 << flush;
	cout << "There were " << real_thread_num << " threads." << endl;

	// imshow("Bild for Filteranwendung", image);
	// waitKey(100000);
	exit(EXIT_SUCCESS);
}