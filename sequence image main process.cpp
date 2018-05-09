//============================================================================
// Name        : CRDIR.cpp
// Author      : Yiyou Chen
// Version     : 14bit image integers input
// Copyright   : Rochester Insitute of Technology Imaging Science Eye Tracking Lab
// Description : Cosmic Ray Damaged Images Repair in C++
//============================================================================

//headfiles
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core/cuda.hpp"
#include <bits/stdc++.h>
#include <fstream>
using namespace cv;
using namespace std;
//end headfiles

int ans;
//functions clarification
void Read_Raw(Mat& RawImage, string Name);
void Find_Isolated_High_Values(Mat& RawImage, Mat& Mask, int MinimumDirections, double MinimumAdjRatio, int Least_Operation_Value);
void Apply_Filters_On_Hot_Pixels(Mat& RawImage, Mat& Mask);
void CannyEdgeDetection_Contour(Mat& NewImage, Mat& Drawing, int CannyParameter);
void Information_Restoration(Mat& FinalImage, Mat& TiffImage, Mat& Drawing);
long long Find_Suitable_Parameter(Mat& NewImage);
void Find_SuperSaturated_Pixels(Mat& Drawing, Mat& SatMask, Mat& FinalImage, int SaturatedRatio);
void Remove_Saturated_Noise(Mat& FinalImage, Mat& SatMask);
void Normalize(Mat& HotPixelImage2);
void Create_Resulting_Image(Mat& ResultingImage, int LeastNumber);
//end function clarification

//main
int main()
{
	cin.tie(0);
	ios_base::sync_with_stdio(false);
	int n(1);
	Mat ComparisonImage(3292, 4992, CV_16UC3, Scalar::all(0)), ResultingImage;
	char ch1, ch2, ch3;
	ch1 = '0';
	for(int p = 663; p <= 726; ++p) {
		string Name = "125H7";
		ch1 = '0' + (p / 100);
		Name = Name + ch1;
		ch2 = '0' + ((p % 100) / 10);
		Name = Name + ch2;
		ch3 = '0' +(p % 10);
		Name = Name + ch3;
		Name = Name + ".NEF.txt.txt";
		Mat HotPixelImage1(3292, 4992, CV_16UC1, Scalar::all(0)), HotPixelImage2(3292, 4992, CV_16UC3, Scalar::all(0));
		Mat RawImage(3292, 4992, CV_16UC1, Scalar::all(0)), Mask(3292, 4992, CV_16UC1, Scalar::all(0)), NewImage, TiffImage, FinalImage, SatMask(3292, 4992, CV_16UC3, Scalar::all(0));
			int MinimumDirections(3), Least_Operation_Value, CannyParameter(300), SaturatedRatio(3), ImageBit(12);
			Least_Operation_Value = 600 * pow(2, ImageBit - 14);
			double MinimumAdjRatio(0.7);
			Read_Raw(RawImage, Name);
			demosaicing(RawImage, TiffImage, COLOR_BayerBG2BGR, 0);
				imwrite("demos125I1430.NEF.tif", TiffImage * 50);
				TiffImage.convertTo(NewImage, CV_8UC3, (double)pow(2.0, 8 - ImageBit));
				Mat Drawing = Mat::zeros(NewImage.size(), CV_8UC3);
				CannyEdgeDetection_Contour(NewImage, Drawing, CannyParameter);
				Find_Isolated_High_Values(RawImage, Mask, MinimumDirections, MinimumAdjRatio, Least_Operation_Value);
				Apply_Filters_On_Hot_Pixels(RawImage, Mask);
				demosaicing(RawImage, FinalImage, COLOR_BayerBG2BGR, 0);
				//imwrite("inter125I1430.NEF.tif", FinalImage * 16);
				Information_Restoration(FinalImage, TiffImage, Drawing);
				//imwrite("drawing125I1430.NEF.tif", Drawing);
				Find_SuperSaturated_Pixels(Drawing, SatMask, FinalImage, SaturatedRatio);
				Remove_Saturated_Noise(FinalImage, SatMask);
				//imwrite("out125I1430.NEF.tif", FinalImage * 50);
		//Normalize(HotPixelImage2);
		for(int i = 0; i < FinalImage.rows; ++i) {
			ushort* FinalImageRowPointer = FinalImage.ptr<ushort>(i);
			ushort* TiffImageRowPointer = TiffImage.ptr<ushort>(i);
			ushort* HotPixelImageRowPointer = HotPixelImage2.ptr<ushort>(i);
			for(int j = 0; j < FinalImage.cols * FinalImage.channels(); ++j) {
				if(FinalImageRowPointer[j] != TiffImageRowPointer[j]) {
					HotPixelImageRowPointer[j] = (ushort)1;
				}
				else HotPixelImageRowPointer[j] = (ushort)0;
			}
		}
		string Name2 = "125I14";
				Name2 = Name2 + ch1;
				Name2 = Name2 + ch2;
				Name2 = Name2 + ch3;
		imwrite(Name2 + "out.tif", HotPixelImage2 * 1000);
		ComparisonImage = ComparisonImage + HotPixelImage2;
	}
	freopen("out.out", "w", stdout);
	//cout << ComparisonImage;
	//imwrite("out.tif", ComparisonImage * 1000);
	ComparisonImage.copyTo(ResultingImage);
	//cout << ComparisonImage;
	//cout << ResultingImage;
	Create_Resulting_Image(ResultingImage, (int)50);
	//cout << ResultingImage;
	//
	imwrite("out.tif", ResultingImage * 5000);
	waitKey(0);
}
//end main

//functions
//Read_Raw: Read a raw image from a text file that contains the information of the raw image;
//functions
//Read_Raw: Read a raw image from a text file that contains the information of the raw image;
void Read_Raw(Mat& RawImage, string Name) {
	ushort RawData;
	ushort* ReadingRowPointer = RawImage.ptr<ushort>(0);
	ifstream File;
	File.open(Name);
	int CurrentCol = 0, CurrentRow = 0;
	//double FloatData;
	while(File >> RawData) {
		//RawData = (ushort)FloatData;
		if(CurrentCol >= RawImage.cols * RawImage.channels()) {
			CurrentCol = 0;
			CurrentRow++;
			ReadingRowPointer = RawImage.ptr<ushort>(CurrentRow);
		}
		ReadingRowPointer[CurrentCol] = RawData;
		CurrentCol++;
	}
	File.close();
}
//end Read_Raw

//Find_Isolated_High_Values: Find hot pixels
//Find_Isolated_High_Values: Find hot pixels
void Find_Isolated_High_Values(Mat& RawImage, Mat& Mask, int MinimumDirections, double MinimumAdjRatio, int Least_Operation_Value) {
	int DirectionX[8] = {-2, -2, -2, 0, 0, 2, 2, 2};
	int DirectionY[8] = {-2 * RawImage.channels(), 0, 2 * RawImage.channels(), -2 * RawImage.channels(), 2 * RawImage.channels(), -2 * RawImage.channels(), 0, 2 * RawImage.channels()};
	ushort* RawCenterRowPointer;
	ushort* MaskCenterRowPointer;
	ushort* MaskTestRowPointer;
	ushort* RawTestRowPointer;
	for(int i = 2; i < RawImage.rows - 2; ++i) {
		RawCenterRowPointer = RawImage.ptr<ushort>(i);
		MaskCenterRowPointer = Mask.ptr<ushort>(i);
		for(int j = 2 * RawImage.channels();  j < (RawImage.cols - 2) * RawImage.channels(); j += RawImage.channels()) {
			bool flag = 0;
			for(int Position = j; Position < j + RawImage.channels(); ++Position) {
				if(RawCenterRowPointer[Position] < Least_Operation_Value) continue;
				int count = 0;
				for(int k = 0; k < 8; ++k) {
					RawTestRowPointer = RawImage.ptr<ushort>(i + DirectionX[k]);
					MaskTestRowPointer = Mask.ptr<ushort>(i + DirectionX[k]);
					if(MaskTestRowPointer[Position + DirectionY[k]] == 65535) continue;
					if((double)RawTestRowPointer[Position + DirectionY[k]] > (double)MinimumAdjRatio * RawCenterRowPointer[Position] || RawCenterRowPointer[Position] - RawTestRowPointer[Position + DirectionY[k]] < 150) {
						count++;
					}
				}
				if(count < 3) {
					flag = 1;
					break;
				}
			}
			if(flag) {
				for(int Position = j; Position < j + RawImage.channels(); ++Position) {
					MaskCenterRowPointer[Position] = 65535;
				}
			}
		}
	}
}
//end Find_Isolated_High_Values

//Apply Average Filter
void Apply_Filters_On_Hot_Pixels(Mat& RawImage, Mat& Mask) {
	int DirectionX[8] = {-2, -2, -2, 0, 0, 2, 2, 2};
	int DirectionY[8] = {-2 * RawImage.channels(), 0, 2 * RawImage.channels(), -2 * RawImage.channels(), 2 * RawImage.channels(), -2 * RawImage.channels(), 0, 2 * RawImage.channels()};
	for(int i = 2; i < Mask.rows - 2; ++i) {
		ushort* RawCenterRowPointer;
		ushort* MaskCenterRowPointer;
		ushort* MaskTestRowPointer;
		ushort* RawTestRowPointer;
		MaskCenterRowPointer = Mask.ptr<ushort>(i);
		RawCenterRowPointer = RawImage.ptr<ushort>(i);
		for(int j = 2 * Mask.channels(); j < (Mask.cols - 2) * Mask.channels(); j += Mask.channels()) {
			for(int Position = j; Position < j + Mask.channels(); ++Position) {
				if(MaskCenterRowPointer[j] == 65535) {
					int sum = 0, num = 0;
					for(int k = 0; k < 8; ++k) {
						MaskTestRowPointer = Mask.ptr<ushort>(i + DirectionX[k]);
						RawTestRowPointer = RawImage.ptr<ushort>(i + DirectionX[k]);
						if(MaskTestRowPointer[Position + DirectionY[k]] != 65535) {
							sum += RawTestRowPointer[Position + DirectionY[k]];
							num++;
						}
					}
					RawCenterRowPointer[Position] = (ushort)sum / num;
					MaskCenterRowPointer[Position] = 0;
				}
			}
		}
	}
}
//End Apply_Filters

//start canny edge detection
void CannyEdgeDetection_Contour(Mat& NewImage, Mat& Drawing, int CannyParameter) {
	GaussianBlur(NewImage, NewImage, Size(0, 0), 0.8 , 0, BORDER_DEFAULT);
	cvtColor(NewImage, NewImage, CV_BGR2GRAY );
	Canny( NewImage, NewImage, (int)0, (int)CannyParameter, 3);
	vector <vector<Point>> contours;
	vector <Vec4i> hierarchy;
	findContours( NewImage, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
	for(int i = 0; i < contours.size(); i++) {
		Scalar color = Scalar(255, 255, 255);
		drawContours(Drawing, contours, i, color, 2, 8, hierarchy, 0, Point());
	}
}
//end canny edge detection


//start information restoration
void Information_Restoration(Mat& FinalImage, Mat& TiffImage, Mat& Drawing) {
	for(int i = 0; i < FinalImage.rows; ++i) {
		ushort* FinalImageRowPointer = FinalImage.ptr<ushort>(i);
		ushort* TiffImageRowPointer = TiffImage.ptr<ushort>(i);
		uchar* DrawingRowPointer = Drawing.ptr<uchar>(i);
		for(int j = 0; j < FinalImage.cols * FinalImage.channels(); j += FinalImage.channels()) {
			bool flag = 0;
			for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
				if(DrawingRowPointer[Position] != (uchar) 255) {
					flag = 1;
					break;
				}
			}
			if(!flag) {
				for(int Position = j; Position < j + FinalImage.channels(); ++Position)
					FinalImageRowPointer[Position] = TiffImageRowPointer[Position];
			}
		}
	}
}
//End information restoration

//Start find suitable parameter
long long Find_Suitable_Parameter(Mat& NewImage) {
	long long result = 0;
	for(int i = 0; i < NewImage.rows; ++i) {
		uchar* NewImageRowPointer = NewImage.ptr<uchar>(i);
		for(int j = 0; j < NewImage.cols * NewImage.channels(); j += NewImage.channels()) {
			for(int Position = j; Position < j + NewImage.channels(); ++Position) {
				result += (int) NewImageRowPointer[Position];
			}
		}
	}
	return result / (NewImage.rows * NewImage.cols * NewImage.channels());
}
//end find suitable parameter

//start get mask
void Find_SuperSaturated_Pixels(Mat& Drawing, Mat& SatMask, Mat& FinalImage, int SaturatedRatio) {
	Mat Mask(FinalImage.rows, FinalImage.cols, CV_16UC3, Scalar::all(0));
	for(int i = 1; i < FinalImage.rows - 1; ++i) {
		ushort* FinalImageRowPointer = FinalImage.ptr<ushort>(i);
		uchar* DrawingRowPointer = Drawing.ptr<uchar>(i);
		for(int j = FinalImage.channels(); j < (FinalImage.cols - 1) * FinalImage.channels(); j += FinalImage.channels()) {
			int MaximumChannel = 0, SumChannels = 0;
			for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
				MaximumChannel = max(MaximumChannel, (int)FinalImageRowPointer[Position]);
				SumChannels += FinalImageRowPointer[Position];
			}
			ushort* Sat_MaskRowPointer = SatMask.ptr<ushort>(i);
			if(MaximumChannel > SaturatedRatio * (SumChannels - MaximumChannel) / (FinalImage.channels() - 1) && MaximumChannel > 75) {
				bool flag = 0;
				for(int Position = j; Position < j + Drawing.channels(); ++Position) {
					if(DrawingRowPointer[Position] != (uchar)0) {
						flag = 1;
						break;
					}
				}
				if(flag) {
					for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
						Sat_MaskRowPointer[Position] = (ushort)1;
					}
				}
				else {
					for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
						Sat_MaskRowPointer[Position] = (ushort)0;
					}
				}
			}
			else {
				for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
					Sat_MaskRowPointer[Position] = (ushort)0;
				}
			}
		}
	}
}
//end get mask

//Start remove_Saturated_noise
void Remove_Saturated_Noise(Mat& FinalImage, Mat& SatMask) {
	Mat Mask(FinalImage.rows, FinalImage.cols, CV_16UC3, Scalar::all(0));
	int DirectionX[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
	int DirectionY[8] = {-FinalImage.channels(), 0, FinalImage.channels(), -FinalImage.channels(), FinalImage.channels(), -FinalImage.channels(), 0, FinalImage.channels()};
	for(int i = 1; i < FinalImage.rows - 1; ++i) {
		ushort* FinalImageRowPointer = FinalImage.ptr<ushort>(i);
		ushort* SatMaskRowPointer = SatMask.ptr<ushort>(i);
		for(int j = FinalImage.channels(); j < (FinalImage.cols - 1) * FinalImage.channels(); j += FinalImage.channels()) {
			for(int Position = j; Position < j + FinalImage.channels(); ++Position) {
				if(SatMaskRowPointer[Position]) {
					int sum = 0, num = 0;
					for(int k = 0; k < 8; ++k) {
						ushort* SatMaskTestRowPointer = SatMask.ptr<ushort>(i + DirectionX[k]);
						ushort* FinalImageTestRowPointer = FinalImage.ptr<ushort>(i + DirectionX[k]);
						if(!SatMaskTestRowPointer[Position + DirectionY[k]]) {
							num++;
							sum += FinalImageTestRowPointer[Position + DirectionY[k]];
						}
					}
					FinalImageRowPointer[Position] = (ushort)sum / num;
					SatMaskRowPointer[Position] = (ushort)0;
				}
			}
		}
	}
}


void Normalize(Mat& HotPixelImage2) {
	for(int i = 0; i < HotPixelImage2.rows; ++i) {
		ushort* HotPixelImageRowPointer = HotPixelImage2.ptr<ushort>(i);
		for(int j = 0; j < HotPixelImage2.cols * HotPixelImage2.channels(); ++j) {
			if(HotPixelImageRowPointer[j] != (ushort)0) HotPixelImageRowPointer[j] = (ushort)1;
		}
	}
}
void Create_Resulting_Image(Mat& ResultingImage, int LeastNumber) {
	for(int i = 0; i < ResultingImage.rows; ++i) {
		ushort* ResultingImageRowPointer = ResultingImage.ptr<ushort>(i);
		for(int j = 0; j < ResultingImage.cols * ResultingImage.channels(); ++j) {
			if((int)ResultingImageRowPointer[j] < (int)LeastNumber) ResultingImageRowPointer[j] = (ushort)0;
		}
	}
}
