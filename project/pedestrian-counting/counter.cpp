/* count_things.cpp
 * a program to count pedestrians and bikers moving through a video
 * hopefully a live stream too

 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * Some of this code is based on helpful tutorials available here:
 * http://docs.opencv.org/3.1.0/d5/d07/tutorial_multitracker.html#gsc.tab=0
 * https://www.youtube.com/user/khounslow/featured

 * If by some miracle you find this software useful, thanks are accepted in
 * the form of chocolate, coffee, or introductions to potential employers.

 */

#include <opencv/cv.hpp>
#include "opencv2/core.hpp"
// #include <raspicam/raspicam_cv.h>
// #include "opencv2/background_segm.hpp"
// #include <opencv/highgui.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "object.h"
#include "useful_functions.h"
#include "image_input.h"
#include "image_output.h"
#include "trackers.h"

#define REMOTE 1 == 1

using namespace std;
using namespace cv;

//our sensitivity value to be used in the threshold() function
static double MAX_DIST_SQD = 6000000; // maximum distance between to centers to consider it one object
static int SENSITIVITY_VALUE_1 = 200; // values for cleaning noise out of difference images
static int SENSITIVITY_VALUE_2 = 50; 
//size of blur used to smooth the image to remove possible noise and
//increase the size of the object we are trying to track. (Much like dilate and erode)
static int BLUR_SIZE_1 = 200;
static int BLUR_SIZE_2 = 200;
static double MIN_OBJ_AREA = 1000;

//TODO don't do this
ImageOutput* video_out; 

void track_with_adaptive_BS(ImageInput* capture, Mat& grayBackground, bool use_static_back,
														double& next_id, int& count_LR, int& count_RL);
void do_adaptive_BS(Ptr<BackgroundSubtractorMOG2> subtractor, Mat &image, bool debugMode, Mat &thresholdImage);
void search_for_movement(Mat &thresholdImage, Mat &display, 
												bool loop_switch, double &next_id, int &count_LR, int &count_RL,
												vector<Object> &objects_0, vector<Object> &objects_1);

void update_object(Object &prev_obj, Object &curr_obj, double mid_row, int &count_LR, int &count_RL);
char is_center_crossed(const Point2d &a, const Point2d &b, double middle);
char is_center_crossed(const Object &obj_a, const Object &obj_b, double middle);

void get_settings_file(int argc, char** argv, string& vid_name, string& back_name, char& bs_type);
void interpret_input(char c, bool &debugMode, bool &trackingEnabled, bool &pause);
void draw_rectangles(vector<Rect2d> &obj_rects, Mat &display);
void draw_centers(vector<Object> &objects, Mat &display);

void show_help();


int main(int argc, char** argv){

	//TODO combine static back and bs type?
	//TODO make this thing a freaking class
	bool use_static_back = false;		  // use a static image for background
	bool background_is_video = true; 	// obtain static back from video
	bool success = false;					    // boolean set when image capture works
	char bs_type = 'N';						    // back subtraction algo 'M' for MOG2, non-adaptive is default
	double next_id = 0;					    	// the next id to use
	int count_LR = 0, count_RL = 0;		// counts of objects
	string vid_name;									// name of video file to use
	string back_name;									// optional name for background file

	Mat grayBackground;
	vector<Object> objects_0, objects_1;
	Mat thresholdImage;
	ImageInput* capture;

	int num_videos = 4;
	char** name_list   = new char*[4]; 
  name_list[0] = (char*)"tracking_video";
  name_list[1] = (char*)"difference_image";
  name_list[2] = (char*)"threshold_image";
  name_list[3] = (char*)"final_threshold";

	if(argc == 2) {
		get_settings_file(argc, argv, vid_name, back_name, bs_type);
	} else {
		show_help();
	}

	//TODO move to get settings?
	Trackers::set_max_dist_sqd(MAX_DIST_SQD);

	if(vid_name == "RASPICAM") {
		capture = new ImageInput();
		cout << "using live stream" << endl;
	} else {
		capture = new ImageInput(vid_name);
		cout << "using video: " << vid_name << endl;
	}
	
	if(use_static_back)
		set_background(back_name, background_is_video, grayBackground, use_static_back);

	success = capture->open();
	if(!success){
			cout << "ERROR ACQUIRING VIDEO FEED named \"" << vid_name << "\"\n";
			getchar();
			exit(1);
	}

	Size S =  Size((int) capture->get(CV_CAP_PROP_FRAME_WIDTH), (int) capture->get(CV_CAP_PROP_FRAME_HEIGHT));
  video_out = new ImageOutput();
  if(!video_out->setup(REMOTE, name_list, S, num_videos))
  	exit(1);

	//TODO we won't need this loop for live streaming
	while(1){

		if(!success){
			cout<<"ERROR ACQUIRING VIDEO FEED\n";
			getchar();
			exit(1);
		} 

		if (bs_type == 'M') {
			cout << endl <<  "Using adaptive (MOG2) Background subtraction" << endl;
			track_with_adaptive_BS(capture, grayBackground, use_static_back, next_id, count_LR, count_RL);
		} else {
			cout << endl <<  "Using non adaptive (Naive) Background subtraction" << endl;
			track_with_non_adaptive_BS(capture, grayBackground, use_static_back, next_id, count_LR, count_RL);
		}

		//release the capture before re-opening and looping again.
		capture->release();
		cout << "Next id for object (total id'd): " << next_id << endl;
		cout << "objects moving Left to Right:    " << count_LR << endl;
		cout << "objects moving Right to Left:    " << count_RL << endl;

		if(REMOTE)
			break;
	
		success = capture->open();

	} // outer while loop (infinite)

	delete[] name_list;
	return 0;
} //main


//track objects through video using GMM background subtraction
//TODO do we actually want gray images for this version?
void track_with_adaptive_BS(ImageInput* capture, Mat& grayBackground, bool use_static_back,
														double& next_id, int& count_LR, int& count_RL) {
	bool debugMode = false;
	bool trackingEnabled = false;
	bool pause = false;
	bool success = false;
	bool loop_switch = true;

	int frames = 1;
	double tot_time = 0;
	clock_t start_t;

	if (REMOTE) {
		trackingEnabled =  true;
		debugMode = true;
	}

	cout << endl;
	cout << "Tracking: " << ((trackingEnabled) ? "Enabled" : "Disabled") << endl;
	cout << "Debug:    " << ((debugMode) ? "Enabled" : "Disabled") << endl;
	cout << endl;

	Ptr<BackgroundSubtractorMOG2> subtractor = createBackgroundSubtractorMOG2();
	Mat frame, image;
	Mat thresholdImage;
	vector<Object> objects_0, objects_1;

	success = capture->read(frame);
	if(!success){
		cout << endl << "ERROR: frame failed to be read" << endl;
		getchar();
		exit(1);
	}

	while( success ) {
		start_t = clock();
		image = frame.clone();
		do_adaptive_BS(subtractor, image, debugMode, thresholdImage);

		if(trackingEnabled) {
			search_for_movement( thresholdImage, frame, loop_switch, next_id, count_LR, count_RL, objects_0, objects_1); 
		}

		char c = video_out->output_track_frame(frame);
		interpret_input(c, debugMode, trackingEnabled, pause);

		success = capture->read(frame);
		loop_switch = !loop_switch;
		tot_time += double(clock() - start_t ) /  CLOCKS_PER_SEC;
		frames++;
	} 
	cout << "Time    = " << tot_time << endl;
	cout << "Frames  = " << frames << endl;
	cout << "t per f = " << tot_time / (double)frames << endl;
}

//@finds movement blobs based on GMM background subtraction
//	also displays the stages if requested
void do_adaptive_BS(Ptr<BackgroundSubtractorMOG2> subtractor, Mat &image, bool debugMode, Mat &thresholdImage) {
	Mat mat_list[3];
	Mat differenceImage, blurImage, firstThreshold;
	int diff = 0, thresh = 1, final = 2; // TODO maybe define these more universally

	subtractor->apply(image, differenceImage);
	//blur(differenceImage, blurImage, Size(BLUR_SIZE_1, BLUR_SIZE_1));
	threshold(differenceImage, firstThreshold, SENSITIVITY_VALUE_1, 255, THRESH_BINARY);

	// TODO determine if this is useful
	blur(firstThreshold, blurImage, Size(BLUR_SIZE_2, BLUR_SIZE_2));
	//TODO paramertize %
	//TODO figure out if this is actually better
	//dynamic_threshold(blurImage, thresholdImage, 0.5, debugMode);
	threshold(blurImage, thresholdImage, SENSITIVITY_VALUE_2, 255, THRESH_BINARY);

	//TODO maybe add blur image
	if(debugMode){
		mat_list[diff]   = differenceImage;
		mat_list[thresh] = firstThreshold;
		mat_list[final]  = thresholdImage;
		video_out->output_debug_frames(mat_list);
		// video_out->output_one_frame_to_file(differenceImage, diff+1);
		// video_out->output_one_frame_to_file(firstThreshold, thresh+1);
		// video_out->output_one_frame_to_file(thresholdImage, final+1);
	}
	else {
		video_out->close_debug_frames();
	}
}

//@identifies objects based on threshold image and previous objects
//@
void search_for_movement(Mat &thresholdImage, Mat &display, 
												bool loop_switch, double &next_id, int &count_LR, int &count_RL,
												vector<Object> &objects_0, vector<Object> &objects_1){

	int obj_count = 0, i = 0;
	double mid_row = (double)(thresholdImage.cols >> 1); // half way across the screen
	double obj_area = 0;
	vector< vector<Point> > contours;
	Mat temp;
	Rect2d temp_rect;
	vector<Rect2d> obj_rects;
	vector<Vec4i> hierarchy;
	Object *prev_obj = NULL;

	thresholdImage.copyTo(temp);

	//TODO clean this up so it makes sense. maybe make some functions
	if(loop_switch){
		findContours(temp, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
		if(contours.size() > 0){
			for(vector< vector<Point> >::iterator it_0 = contours.begin(); it_0 != contours.end(); it_0++) {
				temp_rect = boundingRect(*it_0);
				obj_area = temp_rect.area();

				if(obj_area >= MIN_OBJ_AREA){
					obj_count++;
					obj_rects.push_back(Rect2d(temp_rect));
					objects_0.push_back(Object(*it_0));
					prev_obj = NULL;
					if(objects_1.size() > 0) {
						prev_obj = Trackers::find_previous_object(objects_1, *objects_0.rbegin());
					}
					if(prev_obj == NULL) {
						objects_0.rbegin()->set_id(next_id++);
						objects_0.rbegin()->set_is_counted(false);
					} else {
						update_object(*prev_obj, *objects_0.rbegin(), mid_row, count_LR, count_RL);
					}
				}//if obj_area >= MIN_OBJ_AREA
			draw_centers(objects_0, display);
			draw_centers(objects_1, display);
			objects_1.clear();
			}//outer for
		} //if contour_1 > 0
	} else {
		findContours(temp, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
		if(contours.size() > 0){
			for(vector< vector<Point> >::iterator it_0 = contours.begin(); it_0 != contours.end(); it_0++) {
				temp_rect = boundingRect(*it_0);
				obj_area = temp_rect.area();

				if(obj_area >= MIN_OBJ_AREA){
					obj_count++;
					obj_rects.push_back(Rect2d(temp_rect));
					objects_1.push_back(Object(*it_0));
					prev_obj = NULL;
					if(objects_0.size() > 0) {
						prev_obj = Trackers::find_previous_object(objects_0, *objects_1.rbegin());
					}
					if(prev_obj == NULL) {
						objects_1.rbegin()->set_id(next_id++);
						objects_1.rbegin()->set_is_counted(false);
					} else {
						update_object(*prev_obj, *objects_1.rbegin(), mid_row, count_LR, count_RL);
					}
				}
			draw_centers(objects_0, display);
			draw_centers(objects_1, display);
			objects_0.clear();
			}
		}
	}

	draw_rectangles(obj_rects, display);
	line(display, Point(mid_row, 0), Point(mid_row, display.cols), Scalar( 0, 255, 0 ), 2, 1);

} //search for movement

//@searches through list of old object to find match for the new one
//@returns pointer to the old one
void update_object(Object &prev_obj, Object &curr_obj, double mid_row, int &count_LR, int &count_RL) {
	curr_obj.set_id(prev_obj.get_id());
	curr_obj.set_is_counted(prev_obj.get_is_counted());
	if( !(curr_obj.get_is_counted()) ) {
		switch(is_center_crossed(prev_obj, curr_obj, mid_row)) {
		case 'R':
			curr_obj.set_is_counted();
			count_LR++;
			cout << "object moving Left to Right id = " << curr_obj.get_id() << endl;
			break;
		case 'L':
			curr_obj.set_is_counted();
			count_RL++;
			cout << "object moving Right to Left id = " << curr_obj.get_id() << endl;
			break;
		//default case (implied) - if it isn't R or L dont do anything
		} 
	} 
}

//@checks if the center is crossed
//@returns N- no, L- right to left, R- left to right
char is_center_crossed(const Point2d &a, const Point2d &b, double middle) {
	if( (a.x < middle) && (b.x >= middle) )
		return 'R';
	else if( (b.x < middle) && (a.x >= middle) )
		return 'L';
	else
		return 'N';
}

//@checks if the center is crossed
//@returns N- no, L- right to left, R- left to right
char is_center_crossed(const Object &obj_a, const Object &obj_b, double middle) {
	Point2d a, b;
	obj_a.get_center(a);
	obj_b.get_center(b);
	if( (a.x < middle) && (b.x >= middle) )
		return 'R';
	else if( (b.x < middle) && (a.x >= middle) )
		return 'L';
	else
		return 'N';
}



/*****************************************************************************\

																IO FUNCTIONS


\*****************************************************************************/

//@read file to get  proper settings and file names
//TODO dynamic threshold setting (maybe not)
//TODO tracking algo
void get_settings_file(int argc, char** argv, string& vid_name, string& back_name, char& bs_type) {
	string next_line;
	int input_cnt = 0;
	bool done = false;
	ifstream file;
	file.open(argv[1]);

	if (file.is_open()) {
		while ( getline(file, next_line) && !done ) {
			if(next_line[0] != '#') {
				switch (input_cnt) {
					case 0:
						//TODO handle live stream
						vid_name = next_line.c_str();
						break;
					case 1:
						back_name = next_line.c_str();
						break;
					case 2:
						MAX_DIST_SQD = str_to_int(next_line);
						break;
					case 3:
						SENSITIVITY_VALUE_1 = str_to_int(next_line);
						break;
					case 4:
						SENSITIVITY_VALUE_2 = str_to_int(next_line);
						break;
					case 5:
						BLUR_SIZE_1 = str_to_int(next_line);
						break;
					case 6:
						BLUR_SIZE_2 = str_to_int(next_line);
						break;
					case 7:
						MIN_OBJ_AREA = str_to_int(next_line);
						break;
					case 8:
						bs_type = next_line[0];
						break;
					case 9:
						Trackers::set_algo(next_line[0]);
						break;
				} //switch
				input_cnt++;
			} // if not comment
		} //while
		file.close();
	} else {
		cout << "ERROR: Could nt open configuration file." << endl;
		getchar();
		exit(1);
	}
}

//@interpret keyboard input for runtime options
void interpret_input(char c, bool &debugMode, bool &trackingEnabled, bool &pause) {
	char c2 = 'x';
	bool wait = pause;
	//TODO set defines or somthing for these numbers
	switch(c){
	// case 1048603:
	case 27: //'esc' key has been pressed, exit program.
		cout << "Have a nice day! :)" << endl;
		exit(0);
	// case 1048692:
	case 116: //'t' has been pressed. this will toggle tracking
		trackingEnabled = !trackingEnabled;
		if(trackingEnabled == false) cout << "Tracking disabled." << endl;
		else cout << "Tracking enabled." << endl;
		break;
	// case 1048676:
	case 100: //'d' has been pressed. this will debug mode
		debugMode = !debugMode;
		if(debugMode == false) cout << "Debug mode disabled." << endl;
		else cout << "Debug mode enabled." << endl;
		break;
	// case 1048688:
	case 112: //'p' has been pressed. this will pause/resume the code.
		pause = !pause;
		wait = pause;
		cout << "Code paused, press 'p' again to resume, 's' to step" << endl;
	}

	if(pause == true){
		while (wait){
			//stay in this loop until
			c2 = waitKey(10);
			switch (c2){
			case 112: // p is for unpause
				pause = false;
				wait  = false;
				cout << "Code resumed." << endl;
				break;
			case 115: // s is for step
				pause = true;
				wait  = false;
				break;
			}
		}
	}

}

//@draws the rectagles
void draw_rectangles(vector<Rect2d> &obj_rects, Mat &display) {
	for(unsigned j = 0; j < obj_rects.size(); j++) {
	  rectangle( display, obj_rects[j], Scalar( 255, 0, 0 ), 2, 1 ); // draw rectangle around object
	  // int mid_x = obj_rects[j].x + (obj_rects[j].width / 2);  // was this important?
	  // int mid_y = obj_rects[j].y - (obj_rects[j].height / 2);
	}
}

//@draws the rectagles
void draw_centers(vector<Object> &objects, Mat &display) {
	Point2d temp_pt;
	for(unsigned j = 0; j < objects.size(); j++) {
		objects[j].get_center(temp_pt);
	  circle( display, temp_pt, 5, Scalar( 0, 0, 255 ), 2, 1 );
	  //circle( display, temp_pt, MAX_DIST_SQD, Scalar( 0, 255, 255 ), 2, 1 );
	  putText(display,"Object: " + int_to_str(objects[j].get_id()), temp_pt, 1, 1, Scalar(255,0,0), 2);
	}
}

//@print instructions to standard output and crash program
void show_help() {
  cout << endl <<
  " Usage: ./counter.out <configuration file>\n"
  " example:\n"
  " ./counter.out config_example.txt\n"
  << endl << endl;
  exit(1);
}
