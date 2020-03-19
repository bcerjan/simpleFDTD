#ifndef IMAGE_FUNC
#define IMAGE_FUNC

#include <armadillo>
#include <string>
#include "CImg.h"

void imageOutput(arma::mat &field, cimg_library::CImg<unsigned char> &WallImg, arma::uvec src_pt, std::string filePath);
cimg_library::CImg<unsigned char> findEdge(cimg_library::CImg<unsigned char> inputImg);

#endif
