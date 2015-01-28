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

int main(int argc, char *argv[]) {

    int threadnum = atoi(argv[3]);
    int do_times=1, tiles=0;
    do_times = atoi(argv[1]);
    tiles = atoi(argv[2]);
    cout << "do times "<< do_times<< endl;
//    const int threadnum = 1;
	int real_thread_num= 0;
	vector<thread>threads;
	chrono::system_clock::time_point startTime = chrono::system_clock::now();
    Mat image = imread(string(argv[4]));
	// cout << image.cols << "x" << image.rows << endl;
	omp_set_num_threads(threadnum);
    int i,offset;
     for(int k=0; k < do_times;k++){
#pragma omp parallel for private(offset, i) lastprivate(real_thread_num)
        for(int i=0; i < tiles;i++)
        {

                real_thread_num= omp_get_num_threads();
        //		cout  << "WIDTH" << i << "="<< ((image.cols/4)*(k)) << "-" << ((image.cols/4)*(k+1)) << endl;
        //		cout << "HEIGHT" << i << "="<< ((image.rows/4)*(i)) << "-" << ((image.rows/4)*(i+1)) << endl;
        //		cout << "___________________________________" << endl;
                offset = i>0 ? -2 : 0;
                Mat subImg = image(Range(((image.rows/tiles)*(i))+offset , ((image.rows/tiles)*(i+1))),
                                   Range(0 , image.cols));
                //try{
                    erode(subImg, subImg, Mat(), Point(-1, -1), 1, 0, 0);
                //}catch(...){};
        }
    }
	
	chrono::system_clock::time_point endTime = chrono::system_clock::now();
	chrono::microseconds microRunTime =
			chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	double runTime = microRunTime.count() / 1000000.0;

    cout << setprecision(8) << "time " << runTime << " seconds." << endl
		 << flush;
	cout << "There were " << real_thread_num << " threads." << endl;


    imwrite( "image.jpg", image );
    waitKey(100000);
	exit(EXIT_SUCCESS);
}
