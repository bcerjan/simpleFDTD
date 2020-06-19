# simpleFDTD
An embeddable 2D electromagnetic FDTD simulation for approximate results in the visible region.
Uses the leapfrog alternating-direction implicit FDTD (LADI-FDTD) method. 

The goal of this project is to make a simple, user-friendly, and embeddable 2D FDTD simulation using emscripten to compile to [WebAssembly](https://webassembly.org/). This allows the computation to be run in the user's browser, rather than on the server hosting the page allowing the server to run smoothly and serve simulations even under high load.

To see an example of the generated web page, see [here](https://simplefdtd.waldocorp.com).

## Embedding On Your Site:
To use on your site (subject to the [license](https://github.com/bcerjan/simpleFDTD/blob/master/LICENSE.txt)) just download the three files in the `website` directory and then adjust them as you see fit for your site. It is strongly recommended that you only edit the `.html` file directly unless you know what you are doing.

## Compiling:
To compile, a (perhaps inelegant) combination of `make` and `cmake` is used. `CMake` is used for compiling to executables to run for testing / debugging, while `make` is used to generate the final `.wasm, .html,` and `.js` files.

Specifically, clone the git repo (`git clone https://github.com/bcerjan/simpleFDTD`) then `cd simpleFDTD/c_code`, `mkdir build`, `cd build`, and `cmake ..`. From there, you can now use `make empty` or `make structure` to create executables for testing without compiling all the way to `.wasm`. Note that if you change any parameters, you need to first run `./empty` before compiling `make structure` as `./structure` needs the header files (re)generated by running `./empty`.

To make the `.wasm` (and associated ephemera) enter the `c_code` directory, and run `make`. This compiles to WebAssembly. When you edit things, it can be helpful to
  1. use `make clean_html` to remove the previous version of the `.wasm, .html,` and `.js` files
  1. before running `make`.
