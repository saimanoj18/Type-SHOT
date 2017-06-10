# Type-SHOT
Low Bitrate 3D Feature Descriptors with Lattice Quantization

Source code for paper titled **Low Bit-rate 3D Feature Descriptors for Depth Data from Kinect-style Sensors**

PCL needs to installed in the system!

commands to run the package are
1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./type_shot


**Please look at type_shot.cpp for comments and to know the work flow**



**_Please choose the parameters based on the experimental results reported in the paper_** 

    int m, n, grids, bins;
    m = 22;// this is the length of the subvector from the descriptor that is being compressed!
    n = 3;// this is the number of quatization levels as discussed in the paper
    grids = 16;// grids = dimension of the descriptor/bins, i.e., 352/22 = 16 ( 352 is the dimensionality of SHOT descriptor)
    bins = 22;// this should be equal to `m`, i.e., bins  = m;
