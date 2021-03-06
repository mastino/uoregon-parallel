#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <opencv2/opencv.hpp>

using namespace cv;

struct pixel {
  double red;
  double green;
  double blue;
  
  pixel(double r, double g, double b) : red(r), green(g), blue(b) {};
};


/*
 * The Prewitt kernels can be applied after a blur to help highlight edges
 * The input image must be gray scale/intensities:
 *     double intensity = (in[in_offset].red + in[in_offset].green + in[in_offset].blue)/3.0;
 * Each kernel must be applied to the blured images separately and then composed:
 *     blurred[i] with prewittX -> Xedges[i]
 *     blurred[i] with prewittY -> Yedges[i]
 *     outIntensity[i] = sqrt(Xedges[i]*Xedges[i] + Yedges[i]*Yedges[i])
 * To turn the out intensity to an out color set each color to the intensity
 *     out[i].red = outIntensity[i]
 *     out[i].green = outIntensity[i]
 *     out[i].blue = outIntensity[i]
 *
 * For more on the Prewitt kernels and edge detection:
 *     http://en.wikipedia.org/wiki/Prewitt_operator
 */
// void prewittX_kernel(const int rows, const int cols, double * const kernel) {
// 	if(rows != 3 || cols !=3) {
// 		std::cerr << "Bad Prewitt kernel matrix\n";
// 		return;
// 	}
//   kernel[0] = 1;
//   kernel[1] = 0;
//   kernel[2] = -1;
//   kernel[3] = 1;
//   kernel[4] = 0;
//   kernel[5] = -1;
//   kernel[6] = 1;
//   kernel[7] = 0;
//   kernel[8] = -1;
// }

// void prewittY_kernel(const int rows, const int cols, double * const kernel) {
//   if(rows != 3 || cols !=3) {
//     std::cerr << "Bad Prewitt kernel matrix\n";
//     return;
//   }
//   kernel[0] = 1;
//   kernel[1] = 1;
//   kernel[2] = 1;
//   kernel[3] = 0;
//   kernel[4] = 0;
//   kernel[5] = 0;
//   kernel[6] = -1;
//   kernel[7] = -1;
//   kernel[8] = -1;
// }
void prewittX_kernel(const int rows, const int cols, double * const kernel) {
  if(rows != 3 || cols !=3) {
    std::cerr << "Bad Prewitt kernel matrix\n";
    return;
  }
  for(int i=0;i<3;i++) {
    kernel[0 + (i*rows)] = -1.0;
    kernel[1 + (i*rows)] = 0.0;
    kernel[2 + (i*rows)] = 1.0;
  }
}

void prewittY_kernel(const int rows, const int cols, double * const kernel) {
  if(rows != 3 || cols !=3) {
    std::cerr << "Bad Prewitt kernel matrix\n";
    return;
  }
  for(int i=0;i<3;i++) {
    kernel[i + (0*rows)] = 1.0;
    kernel[i + (1*rows)] = 0.0;
    kernel[i + (2*rows)] = -1.0;
  }
}

/*
 * The gaussian kernel provides a stencil for blurring images based on a 
 * normal distribution
 */
void gaussian_kernel(const int rows, const int cols, const double stddev, double * const kernel) {
  const double denom = 2.0 * stddev * stddev;
  const double g_denom = M_PI * denom;
  const double g_denom_recip = (1.0/g_denom);
  double sum = 0.0;
  for(int i = 0; i < rows; ++i) {
    const double row_dist = i - (rows/2);
    for(int j = 0; j < cols; ++j) {
      const double col_dist = j - (cols/2);
      const double dist_sq = (row_dist * row_dist) + (col_dist * col_dist);
      const double value = g_denom_recip * exp((-dist_sq)/denom);
      kernel[i + (j*rows)] = value;
      sum += value;
    }
  }
  // Normalize
  const double recip_sum = 1.0 / sum;
  for(int i = 0; i <= rows-2; i+=2) {
    #pragma ivdep
    for(int j = 0; j <= cols-5; j+=5) {
      kernel[(i) + ((j)*rows)] *= recip_sum;
      kernel[(i) + ((j+1)*rows)] *= recip_sum;
      kernel[(i) + ((j+2)*rows)] *= recip_sum;
      kernel[(i) + ((j+3)*rows)] *= recip_sum;
      kernel[(i) + ((j+4)*rows)] *= recip_sum;
    }   
    #pragma ivdep
    for(int j = 0; j <= cols-5; j+=5) {
      kernel[(i+1) + ((j)*rows)] *= recip_sum;
      kernel[(i+1) + ((j+1)*rows)] *= recip_sum;
      kernel[(i+1) + ((j+2)*rows)] *= recip_sum;
      kernel[(i+1) + ((j+3)*rows)] *= recip_sum;
      kernel[(i+1) + ((j+4)*rows)] *= recip_sum;
    }    
  }
}

void apply_stencil(const int radius, const double stddev,
       const int rows, const int cols, pixel * const in, pixel * const out) {
  const int dim = radius*2+1;
  double kernel[dim*dim];
  gaussian_kernel(dim, dim, stddev, kernel);
  #pragma omp parallel for
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < cols; ++j) {
      const int out_offset = i + (j*rows);
      // For each pixel, do the stencil
      for(int x = i - radius, kx = 0; x <= i + radius; ++x, ++kx) {
        for(int y = j - radius, ky = 0; y <= j + radius; ++y, ++ky) {
          if(x >= 0 && x < rows && y >= 0 && y < cols) {
            const int in_offset = x + (y*rows);
            const int k_offset = kx + (ky*dim);
            out[out_offset].red   += kernel[k_offset] * in[in_offset].red;
            out[out_offset].green += kernel[k_offset] * in[in_offset].green;
            out[out_offset].blue  += kernel[k_offset] * in[in_offset].blue;
          }
        }
      }
    }
  } //main for
}

void apply_prewitt(const int rows, const int cols, pixel * const in, double *  out) {

  int ind;
  const int radius = 1;
  const int dim = radius*2+1;
  double intensity[rows*cols];
  double Xedges[rows*cols];
  double Yedges[rows*cols];
  double Xkernel[9];
  double Ykernel[9];

  prewittX_kernel(dim, dim, Xkernel);
  prewittY_kernel(dim, dim, Ykernel);

  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < cols; ++j) {
      ind = i + (j*rows);
      intensity[ind] = (in[ind].red + in[ind].green + in[ind].blue)/3.0;
      Xedges[ind] = 0;
      Yedges[ind] = 0;
    }
  }

  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < cols; ++j) {
      const int out_offset = i + (j*rows);
      for(int x = i - radius, kx = 0; x <= i + radius; ++x, ++kx) {
        for(int y = j - radius, ky = 0; y <= j + radius; ++y, ++ky) {
          if(x >= 0 && x < rows && y >= 0 && y < cols) {
            const int in_offset = x + (y*rows);
            const int k_offset = kx + (ky*dim);
            Xedges[out_offset] += Xkernel[k_offset] * intensity[in_offset];
            Yedges[out_offset] += Ykernel[k_offset] * intensity[in_offset];
          }
        }
      }
    }
  } 

  for(int i = 0; i <= rows-1; i+=1) {
    for(int j = 0; j <= cols-1; j+=1) {
      ind = i + (j*rows);
      out[ind] = sqrt(Xedges[ind]*Xedges[ind] + Yedges[ind]*Yedges[ind]);
    }
  } 

}

int main( int argc, char* argv[] ) {
  double start, runtime;


  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " imageName\n";
    return 1;
  }

  // Read image
  Mat image;
  image = imread(argv[1], CV_LOAD_IMAGE_COLOR);
  if(!image.data ) {
    std::cout <<  "Error opening " << argv[1] << std::endl;
    return -1;
  }
	
  // Get image into C array of doubles for processing
  const int rows = image.rows;
  const int cols = image.cols;
  pixel * imagePixels = (pixel *) malloc(rows * cols * sizeof(pixel));

  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < cols; ++j) {
      Vec3b p = image.at<Vec3b>(i, j);
      imagePixels[i + (j*rows)] = pixel(p[0]/255.0,p[1]/255.0,p[2]/255.0);
    }
  }
	
  // Create output array
  pixel * outPixels = (pixel *) malloc(rows * cols * sizeof(pixel));
  for(int i = 0; i < rows * cols; ++i) {
    outPixels[i].red = 0.0;
    outPixels[i].green = 0.0;
    outPixels[i].blue = 0.0;
  }

  double * out = (double *) malloc(rows * cols * sizeof(double));
	
  // Do the stencil
  start = omp_get_wtime(); 
  apply_stencil(3, 32.0, rows, cols, imagePixels, outPixels);
  apply_prewitt(rows, cols, outPixels, out);
  runtime = omp_get_wtime() - start;

  // Create an output image (same size as input)
  Mat dest(rows, cols, CV_8UC1);
  // Copy C array back into image for output
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < cols; ++j) {
      dest.at<char>(i, j) = (char)floor(out[i + (j*rows)] * 255.0);
    }
  }
	
  imwrite("out.jpg", dest);
  printf("  runtime for apply_stencil was %f sec.\n", runtime);
	
  free(imagePixels);
  free(outPixels);
  free(out);
  return 0;
}
