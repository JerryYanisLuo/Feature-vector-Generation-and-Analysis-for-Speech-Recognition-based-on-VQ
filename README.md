# Feature-vector-Generation-and-Analysis-for-Speech-Recognition-based-on-VQ
**1.introduction**

This project is implemented by C programming, including feature extraction, feature classification and result evaluation. The input can be any mono wav file with 16 bps. The program will process the wav file to yield feature vectors, and generate codebook by using LBG algorithm. The evaluation part is realized by using PCA to produce the two-dimensional coordinates of the feature vectors and codebook. Additionally, a re-arranged test wav file will be produced as well in which the signal segments from the same cluster are gathered together. All the results are recorded as txt files.
  
**2.Operating environment**

gcc complier with version 8.2.0 or higher.

**3.Operating method**

This project mainly contains 3 header files and 4 c files,
* wave.h
* LBG.h
* PCA.h
* wave.c
* lbg.c
* pca.c
* main.c

The input wav files are all stored in the folder called "wav". In the file main.c, there are several parameters which can be adjusted, which has been explained in the comments. After set suitable values for each parameter, the joint compilation should be used as follows, gcc wave.c lbg.c pca.c main.c -o output.exe

After operating the executable file, the output results are stored in different txt files. By using Scilab or any other related sofeware, the data can be plotted in two-dimensional graph.

**4.Change log**

2019/5/31 ver1.0.0 First Version
