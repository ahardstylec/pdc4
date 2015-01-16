/*Modified by Jakob Österling and Fred Fluegge to fit our needs*/


/* This is a code by Team SAHE, India */

/* *************** Stereo Camera Calibration in Real Time **************
 This code can be used to calibrate stereo camera or two cameras to get the intrinsic
 and extrinsic files.
 This code also generated rectified image, and also shows RMS Error and Reprojection error
 to find the accuracy of calibration.
 This code captures stereo images from two different cameras (or stereo camera), whose
 You can set no of stereo pairs you want to use bby editing 'noOfStereoPairs' global
 variable.
 Cheers
 Abhishek Upperwal
 */

/* FF License Text was annoyingly long, look it up on google, its BSD, all changes lisenced unter WTFPL, no warranty what so ever*/

#include "calib.h" //interface to main program
#include "stopwatch.h" //performance messurement

#define timeGap 3000000000U

bool firstRun=true;

using namespace cv;
using namespace std;

//globals
enum Modes { CAPTURING, CALIBRATING};
Modes mode = CAPTURING;
const int noOfStereoPairs = 14;
int stereoPairIndex = 0, cornerImageIndex=0;
int goIn = 1;
Mat _leftOri, _rightOri;
int64 prevTickCount;
vector<Point2f> cornersLeft, cornersRight;
vector<vector<Point2f> > cameraImagePoints[2];

stopwatch probe1,probe2,probe3;

//proto
Mat displayCapturedImageIndex(Mat);
bool findChessboardCornersAndDraw(Mat, Mat, Size);
void displayImages();
void saveImages(Mat, Mat, int);
void calibrateStereoCamera(Size, Size);

//functions
Mat displayCapturedImageIndex(Mat img) {
	std::ostringstream imageIndex;
	imageIndex<<stereoPairIndex<<"/"<<noOfStereoPairs;
	putText(img, imageIndex.str().c_str(), Point(50, 70), FONT_HERSHEY_PLAIN, 0.9, Scalar(0,0,255), 2);
	return img;
}


bool findChessboardCornersAndDraw(Mat inputLeft, Mat inputRight, Size boardSize) {
	_leftOri = inputLeft;
	_rightOri = inputRight;
	bool foundLeft = false, foundRight = false;

	//convert to gray.. (fast)
	cvtColor(inputLeft, inputLeft, COLOR_BGR2GRAY);
	cvtColor(inputRight, inputRight, COLOR_BGR2GRAY);

	//find corners (super slow)
	//you realy need to play with the flags, highly unreliable
	probe1.startClock();
	foundLeft = findChessboardCorners(inputLeft, boardSize, cornersLeft, CALIB_CB_ADAPTIVE_THRESH| CALIB_CB_NORMALIZE_IMAGE| CALIB_CB_ADAPTIVE_THRESH);
	//   | CALIB_CB_FAST_CHECK
	probe1.stopClock();
	probe2.startClock();
	foundRight = findChessboardCorners(inputRight, boardSize, cornersRight, CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_FAST_CHECK);
	probe2.stopClock();

	//draw corners (fast)
	drawChessboardCorners(_leftOri, boardSize, cornersLeft, foundLeft);
	drawChessboardCorners(_rightOri, boardSize, cornersRight, foundRight);

	_rightOri = displayCapturedImageIndex(_rightOri);

	return foundLeft&foundRight;
}

void displayImages(int CAPW) {
	imshow("Left Image", _leftOri);
	imshow("Right Image", _rightOri);
	if (firstRun) {
		moveWindow("Left Image", CAPW*2, 23);
		moveWindow("Right Image", CAPW, 23);
		firstRun=false;
	}
}

void saveImages(Mat leftImage, Mat rightImage, int pairIndex) {
	cameraImagePoints[0].push_back(cornersLeft);
	cameraImagePoints[1].push_back(cornersRight);
	cvtColor(leftImage, leftImage, COLOR_BGR2GRAY);
	cvtColor(rightImage, rightImage, COLOR_BGR2GRAY);
	std::ostringstream leftString, rightString;
	leftString<<"left"<<pairIndex<<".jpg";
	rightString<<"right"<<pairIndex<<".jpg";
	imwrite(leftString.str().c_str(), leftImage);
	imwrite(rightString.str().c_str(), rightImage);
}

void calibrateStereoCamera(Size boardSize, Size imageSize) {
	vector<vector<Point3f> > objectPoints;
//	Mat* objectPoints = (Mat*)new vector<vector<Point3f>>;
	objectPoints.resize(noOfStereoPairs);
	cerr << "foobar" << endl << flush;
	for (int i=0; i<noOfStereoPairs; i++) {
		for (int j=0; j<boardSize.height; j++) {
			for (int k=0; k<boardSize.width; k++) {
				objectPoints[i].push_back(Point3f(float(j),float(k),0.0));
			}
		}
	}
	Mat cameraMatrix[2], distCoeffs[2];
	cameraMatrix[0] = Mat::eye(3, 3, CV_64F);
	cameraMatrix[1] = Mat::eye(3, 3, CV_64F);
	Mat R, T, E, F;
//	CvMat CvObjectPoints = *objectPoints;
//	CvMat CvCameraImagePoints0 = Mat(cameraImagePoints[0]);
//	CvMat CvCameraImagePoints1 = Mat(cameraImagePoints[1]);
//	CvMat CvCameraMatrix0 = cameraMatrix[0];
//	CvMat CvDistCoeffs0 = distCoeffs[0];
//	CvMat CvCameraMatrix1 = cameraMatrix[1];
//	CvMat CvDistCoeffs1 = distCoeffs[1];
//	Size Sizee = imageSize;
//	CvMat tmp_npoints;////////////Needs to be something useful
//	CvMat CvR = R;
//	CvMat CvT = T;
//	CvMat CvE =	E;
//	CvMat CvF	 = F;
CvTermCriteria CvTermCriteriaa;
	CvTermCriteriaa.type = TermCriteria::COUNT + TermCriteria::EPS;
	CvTermCriteriaa.max_iter = 100;
	CvTermCriteriaa.epsilon = 1e-5;
	double rms = stereoCalibrate(objectPoints,cameraImagePoints[0],cameraImagePoints[1],cameraMatrix[0],distCoeffs[0],cameraMatrix[1],distCoeffs[1],imageSize,R,T,E,F);
	//double rms = cvStereoCalibrate(&CvObjectPoints, &CvCameraImagePoints0, &CvCameraImagePoints1, &tmp_npoints,
//								 &CvCameraMatrix0, &CvDistCoeffs0,
//								 &CvCameraMatrix1, &CvDistCoeffs1,
//								 Sizee, &CvR, &CvT, &CvE, &CvF,CvTermCriteriaa,
//								 CV_CALIB_FIX_ASPECT_RATIO +
//								 CV_CALIB_ZERO_TANGENT_DIST +
//								 CV_CALIB_SAME_FOCAL_LENGTH +
//								 CV_CALIB_RATIONAL_MODEL +
//								 CV_CALIB_FIX_K3 +CV_CALIB_FIX_K4 + CV_CALIB_FIX_K5
//								  );
	cout<<"RMS Error: "<<rms<<"\n";
	double err = 0;
	int npoints = 0;
	vector<Vec3f> lines[2];
	for(int i = 0; i < noOfStereoPairs; i++ )
	{
		int npt = (int)cameraImagePoints[0][i].size();
		Mat imgpt[2];
		for(int k = 0; k < 2; k++ )
		{
			imgpt[k] = Mat(cameraImagePoints[k][i]);
			undistortPoints(imgpt[k], imgpt[k], cameraMatrix[k], distCoeffs[k], Mat(), cameraMatrix[k]);
			computeCorrespondEpilines(imgpt[k], k+1, F, lines[k]);
			cout << "imgpt[k] = " << imgpt[k] << endl;
		}
		for(int j = 0; j < npt; j++ )
		{
			double errij = fabs(cameraImagePoints[0][i][j].x*lines[1][j][0] +
								cameraImagePoints[0][i][j].y*lines[1][j][1] + lines[1][j][2]) +
			fabs(cameraImagePoints[1][i][j].x*lines[0][j][0] +
				 cameraImagePoints[1][i][j].y*lines[0][j][1] + lines[0][j][2]);
			err += errij;
			cout << err << "= Err" << endl << flush;
		}
		npoints += npt;
		cout << "i = "<< i << endl << flush;
	}
	cout << "Average Reprojection Error: " <<  err/npoints << endl;
	FileStorage fs("intrinsics.yml", FileStorage::WRITE);
	if (fs.isOpened()) {
		fs << "M1" << cameraMatrix[0] << "D1" << distCoeffs[0] <<
		"M2" << cameraMatrix[1] << "D2" << distCoeffs[1];
		fs.release();
	}
	else
		cout<<"Error: Could not open intrinsics file.";
	cout << "File initristic.yml should be writen" << endl << flush;
	Mat R1, R2, P1, P2, Q;
	Rect validROI[2];
	stereoRectify(cameraMatrix[0], distCoeffs[0], cameraMatrix[1], distCoeffs[1], imageSize, R, T, R1, R2, P1, P2, Q, CALIB_ZERO_DISPARITY, 1, imageSize, &validROI[0], &validROI[1]);
	fs.open("extrinsics.yml", FileStorage::WRITE);
	if (fs.isOpened()) {
		fs << "R" << R << "T" << T << "R1" << R1 << "R2" << R2 << "P1" << P1 << "P2" << P2 << "Q" << Q;
		fs.release();
	}
	else
		cout<<"Error: Could not open extrinsics file";
	cout << "File extristic.yml should be writen" << endl << flush;
	bool isVerticalStereo = fabs(P2.at<double>(1, 3)) > fabs(P2.at<double>(0, 3));
	Mat rmap[2][2];
	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, P1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, P2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
	Mat canvas;
	double sf;
	int w, h;
	if (!isVerticalStereo) {
		sf = 600./MAX(imageSize.width, imageSize.height);
		w = cvRound(imageSize.width*sf);
		h = cvRound(imageSize.height*sf);
		canvas.create(h, w*2, CV_8UC3);
	}
	else {
		sf = 300./MAX(imageSize.width, imageSize.height);
		w = cvRound(imageSize.width*sf);
		h = cvRound(imageSize.height*sf);
		canvas.create(h*2, w, CV_8UC3);
	}
	String file;
	namedWindow("rectified");
	for (int i=0; i < noOfStereoPairs; i++) {
		for (int j=0; j < 2; j++) {
			if (j==0) {
				file = "left";
			}
			else if (j==1) {
				file = "right";
			}
			ostringstream st;
			st<<file<<i+1<<".jpg";
			Mat img = imread(st.str().c_str()), rimg, cimg;
			remap(img, rimg, rmap[j][0], rmap[j][1], INTER_LINEAR);
			cimg=rimg;
			Mat canvasPart = !isVerticalStereo ? canvas(Rect(w*j, 0, w, h)) : canvas(Rect(0, h*j, w, h));
			resize(cimg, canvasPart, canvasPart.size(), 0, 0, INTER_AREA);
			Rect vroi(cvRound(validROI[j].x*sf), cvRound(validROI[j].y*sf),
					  cvRound(validROI[j].width*sf), cvRound(validROI[j].height*sf));
			rectangle(canvasPart, vroi, Scalar(0,0,255), 3, 8);
		}
		if( !isVerticalStereo )
			for(int j = 0; j < canvas.rows; j += 16 )
				line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
		else
			for(int j = 0; j < canvas.cols; j += 16 )
				line(canvas, Point(j, 0), Point(j, canvas.rows), Scalar(0, 255, 0), 1, 8);
		imshow("rectified", canvas);
	}
}

void calibCams(int camL, int camR,int CAPW, int CAPH){
	VideoCapture camLeft(camL), camRight(camR);

	probe1.init("probe1");
	probe2.init("probe2");
	probe3.init("probe3");

	//downsize cam stream
	camLeft.set(CV_CAP_PROP_FRAME_WIDTH, CAPW);
	camLeft.set(CV_CAP_PROP_FRAME_HEIGHT, CAPH);
	camRight.set(CV_CAP_PROP_FRAME_WIDTH, CAPW);
	camRight.set(CV_CAP_PROP_FRAME_HEIGHT, CAPH);

	//fail if no cams
	if (!camLeft.isOpened() || !camRight.isOpened()){
		cout << "No cameras found" << endl;
		exit(EXIT_FAILURE);
	}

	//chessboard search pattern
	Size boardSize = Size(3,3);
	Mat inputLeft, inputRight, copyImageLeft, copyImageRight;
	bool foundCornersInBothImage = false;
	namedWindow("Left Image");
	namedWindow("Right Image");

	// DA LOOP
	while (true) {
		//capture and flip (2ms)
		camLeft>>inputLeft;camRight>>inputRight;
		flip(inputLeft, inputLeft,-1);

		//0.2ms
		inputLeft.copyTo(copyImageLeft);
		inputRight.copyTo(copyImageRight);

		//superslow 800ms
		foundCornersInBothImage = findChessboardCornersAndDraw(inputLeft, inputRight, boardSize);

		//superfast if idle 200ns
		//used to be capturing
		if (foundCornersInBothImage && stereoPairIndex<14) {
			int64 thisTick = getTickCount();
			int64 diff = thisTick - prevTickCount;
			if (goIn==1 || diff >= timeGap) {
				goIn=0;
				saveImages(copyImageLeft, copyImageRight, ++stereoPairIndex);
				prevTickCount = getTickCount();
			}
		}

		//fast 5ms
		displayImages(CAPW);
		waitKey(1);

		//used to be p mode
		if (stereoPairIndex==14) {
			calibrateStereoCamera(boardSize, Size(inputLeft.size()));
			waitKey();
			exit(EXIT_SUCCESS);
		}
	}
}
