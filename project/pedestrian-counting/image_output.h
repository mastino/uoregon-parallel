/* image_output.h
 * a class to make it easier to switch between saving opencv output
 * and outputing to the screen

 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * If by some miracle you find this software useful, thanks are accepted in
 * the form of chocolate or introductions to potential employers.

 */


#ifndef IMAGE_OUTPUT_H
#define IMAGE_OUTPUT_H

#include <opencv/cv.hpp>
#include "opencv2/core.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "useful_functions.h"


using namespace cv;
using namespace std;

//TODO stop being lazy and don't hardcode this stuff
#define WIN_HIEGHT 384
#define WIN_LENGTH 512

class ImageOutput {
public:
  //@constructor with default values
  ImageOutput();
  //@constructor with specific values
  ImageOutput(bool new_to_file, char** new_name_list, Size new_size, int new_video_count);
  //@destructor
  ~ImageOutput();

  //@updates data members and setup video files if neceessary
  //@params
  //@return true if success false for failure
  bool setup(bool new_to_file, char** new_name_list, Size new_size, int new_video_count);

  //@outputs the requested frames
  //@track returns the result of a wait key
  //@usage output_trackframe should be called shortly after any others
  char output_track_frame(Mat &frame);
  bool output_debug_frames(Mat frame[]);
  void output_one_frame_to_file(Mat &frame, int i);

  //@closes output
  void close_track_frame();
  void close_debug_frames();

  //@getters (no setters only initial setup)
  bool get_to_file() const; 
  int  get_video_count() const;

private:
  bool to_file;  // whether the output should be written to file or screen
                 // note: cannot be changed after initial setup

  VideoWriter** video_list;  // list of video files for output
  char**        name_list;   // list of video names; if file, .h264 will be added
                             // note: first in list is main (tracking) video
  Size          frame_size;  // size for output frames, only needed if to_file
  int           video_count; // how many videos there are
};

#endif