#include "SDL.h"
#include "fdtd_proto.h"
#include "sdl_funcs.h"
#include <stdio.h>
#include <stdlib.h>

/* Constants for SDL Rendering, set by imageInit */
static SDL_Window *MainWindow;
static SDL_Renderer *renderer;

void imageInit(struct Grid *g) {

  const unsigned int WINDOW_WIDTH = xSize ;
  const unsigned int WINDOW_HEIGHT = ySize ;
  // Initialize SDL:
  SDL_Init(SDL_INIT_VIDEO);

  // Don't eat all keyboard inputs:
  SDL_SetHint(SDL_HINT_EMSCRIPTEN_KEYBOARD_ELEMENT, "#canvas");

  // Window:
  MainWindow = SDL_CreateWindow("simpleFDTD",
                               SDL_WINDOWPOS_CENTERED,
                               SDL_WINDOWPOS_CENTERED,
                               WINDOW_WIDTH, WINDOW_HEIGHT,
                               SDL_WINDOW_SHOWN
                               );

  // Renderer
  renderer = SDL_CreateRenderer(MainWindow, -1, SDL_RENDERER_PRESENTVSYNC);
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
  SDL_RenderClear(renderer); // Make scene Black

  return;
}

void imageFree() {
  SDL_DestroyRenderer(renderer);
  return;
}

/** Function to find walls using arma mat instead of image: **/
void findMatEdge(struct Grid *g) {
  int i,j,n,k;

  //double **out = AllocateMemory(sizeX, sizeY, 0.0);
  double Gx;
  double Gy;
  double gx[3][3];
  double gy[3][3];

  // Specify Sobel Operators:
  gx[0][0] = -1.0;
  gx[0][2] = 1.0;
  gx[1][0] = -2.0;
  gx[1][2] = 2.0;
  gx[2][0] = -1.0;
  gx[2][2] = 1.0;

  gy[0][0] = 1.0;
  gy[0][1] = 2.0;
  gy[0][2] = 1.0;
  gy[2][0] = -1.0;
  gy[2][1] = -2.0;
  gy[2][2] = -1.0;

/*
  //Sobel Operators:
  mat gx = { {-1, 0, 1},
             {-2, 0, 2},
             {-1, 0, 1} };

  mat gy = { {1, 2, 1},
             {0, 0, 0},
             {-1, -2, -1} };
*/
  for (i = 1; i < xSize-3; i++) {
    for (j = 1; j < ySize-3; j++) {
      Gx = 0.0;
      Gy = 0.0;
      // Manual 2-D convolution:
      for (n = -1; n < 2; n++) { // CHECK ME. MIGHT BE OFF-BY-ONE!!
        for (k = -1; k < 2; k++) {
          Gx += gx[n+1][k+1] * object_locs[i+n][j+k];
          Gy += gy[n+1][k+1] * object_locs[i+n][j+k];
        } /* kForLoop */
      } /* nForLoop */

      edgeMat[i+1][j+1] = pow(Gx*Gx + Gy*Gy, 0.5);
      //printf("edgeMat[%i][%i]: %f\n",i+1,j+1,edgeMat[i+1][j+2]);
    } /* jForLoop */
  } /* iForLoop */

  return;

}

/** Function to plot E^2 Field on our grid: **/
void PlotField (struct Grid *g, double maximumValue, double minimumValue) {
  int i,j;
  double scaleValue;

  /***********************************************************************/
  //     Plot fields
  /***********************************************************************/
  // Send E^2 field to our plotting routine:
  scaleValue = 256.0 / (maximumValue - minimumValue);
  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      e2Field[i][j] = scaleValue * ( ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j] );
    } /* xForLoop */
  } /* yForLoop */

  imageShow(g);

}


/** Function to visualize standard C array using SDL **/
/** Scales the size of the array so that 1 pixel in the array -> 4 pixels
    in the image **/

void imageShow(struct Grid *g) {
  // Initialization:
  int pitch;
  int scaleFactor = 2; // scaleFactor is squared to get number of output pixels
  void *pixels;
  Uint8 *base;
  SDL_Texture *texture = NULL;

  //const unsigned int WINDOW_WIDTH = xSize * scaleFactor;
  //const unsigned int WINDOW_HEIGHT = ySize * scaleFactor;
  //printf("width: %i\n", WINDOW_WIDTH);
  //printf("height: %i\n", WINDOW_HEIGHT);
  const unsigned int WINDOW_WIDTH = xSize ;
  const unsigned int WINDOW_HEIGHT = ySize ;

  const double CENTER_X = WINDOW_WIDTH/2.0;
  const double CENTER_Y = WINDOW_HEIGHT/2.0;

  float z,b,r;

  unsigned int xc, yc;



  texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                              SDL_TEXTUREACCESS_STREAMING,
                              WINDOW_WIDTH, WINDOW_HEIGHT);

  SDL_LockTexture(texture, NULL, &pixels, &pitch);

  for( int i = 0; i<WINDOW_WIDTH; i+=1 ) {
    for( int j = 0; j<WINDOW_HEIGHT; j+=1 ) {
      //printf( "[i][j]: [%i][%i]\n", i,j );
      xc = CENTER_X - j;
      yc = CENTER_Y - i;
      z = (float )e2Field[i][j]; // Get value at this pixel (scaled)
      b = 0.0;
      r = 0.0;
      if (edgeMat[i][j] > 4.0) {
        b = 160.0;
        z = 160.0;
        r = 160.0;
      }
      //z = 200.0* (float )object_locs[i][j]; // Get value at this pixel (scaled)
      //z = 125.0;
      //float [r,g,b] = colormap(z)
      // or base[0],base[1],base[2] = colormap(z) or something...
      /*for( int k = 0; k < scaleFactor; k++) {
        for (int l = 0; l < scaleFactor; l++) {
          if(j == 10 && k == 0 && l == 0) {
            printf("fixed j i,j: %i,%i\n",i,j,k,l );
          }
          if(i == 10 && k == 0 && l == 0) {
            printf("fixed i i,j: %i,%i\n",i,j,k,l );
          }*/
          //base = ((Uint8 *)pixels) + (4*((j+l)*WINDOW_WIDTH + (i+k)));
          base = ((Uint8 *)pixels) + (4*((j)*WINDOW_WIDTH + (i)));
          base[0] = b; // Blue
          base[1] = z; // Green
          base[2] = r; // Red
          base[3] = 255; // Opacity? 255 => totally opaque
        //} /* lForLoop */
      //} /* kForLoop */
    } /* jForLoop */
  } /* iForLoop */

  SDL_UnlockTexture(texture);
  SDL_RenderCopy(renderer, texture, NULL, NULL);

  SDL_RenderPresent(renderer);

  // Memory Management:
  //SDL_DestroyRenderer(renderer);
  //SDL_DestroyTexture(texture); // Not needed -- above line should free the texture too
  //SDL_DestroyWindow(MainWindow);
  return;

}
