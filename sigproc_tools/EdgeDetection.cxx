#ifndef __SIGPROC_TOOLS_EDGEDETECTION_CXX__
#define __SIGPROC_TOOLS_EDGEDETECTION_CXX__

#include "EdgeDetection.h"


long sigproc_tools::EdgeDetection::CantorEnum(const int &x, const int &y) const
{
  int n = ((x + y) * (x + y + 1) / 2) + y;
  return n;
}


void sigproc_tools::EdgeDetection::Convolve2D(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const Array2D<float>& kernel) const
{
  // Input kernel must be normalized. 

  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  int kernelX = kernel.size() / 2;
  int kernelY = kernel.at(0).size() / 2;

  int kernelWidth = kernel.size();
  int kernelHeight = kernel.at(0).size();

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {

      for (int m=0; m<kernelWidth; ++m) {
        for (int n=0; n<kernelHeight; ++n) {
          
          int offsetX = m - kernelX;
          int offsetY = n - kernelY;

          int ix = i + offsetX;
          int iy = j + offsetY;

          if (0 <= ix && ix < numChannels && 0 <= iy && iy < numTicks) {
            output2D[i][j] += kernel[m][n] * input2D[ix][iy];
          }
        }
      }
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::SobelX(
  const Array2D<float>& input2D,
  Array2D<float>& gradient) const
{
  Array2D<float> kernel(3);
  for (auto& v : kernel) {
    v.resize(3);
  }
  // Define SobelX
  kernel[0][0] = 1.0;
  kernel[0][1] = 0.0;
  kernel[0][2] = -1.0;
  kernel[1][0] = 2.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = -2.0;
  kernel[2][0] = 1.0;
  kernel[2][1] = 0.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);
  return;
}


void sigproc_tools::EdgeDetection::SobelY(
  const Array2D<float>& input2D,
  Array2D<float>& gradient) const
{
  Array2D<float> kernel(3);
  for (auto& v : kernel) {
    v.resize(3);
  }

  // Define SobelY
  kernel[0][0] = 1.0;
  kernel[0][1] = 2.0;
  kernel[0][2] = 1.0;
  kernel[1][0] = 0.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = 0.0;
  kernel[2][0] = -1.0;
  kernel[2][1] = -2.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);
  return;
}


void sigproc_tools::EdgeDetection::Sobel(
  const Array2D<float>& input2D,
  Array2D<float>& sobelX,
  Array2D<float>& sobelY,
  Array2D<float>& gradient,
  Array2D<float>& direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  SobelXFast(input2D, sobelX);
  SobelYFast(input2D, sobelY);

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      float g = sqrt(sobelX[i][j]*sobelX[i][j] + sobelY[i][j]*sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::SobelXFast(
  const Array2D<float> &input2D,
  Array2D<float> &gradient) const
{
  const size_t numChannels = input2D.size();
  const size_t numTicks = input2D.at(0).size();

  assert(gradient.size() == numChannels);
  assert(gradient.at(0).size() == numTicks);

  Array2D<float> buffer(numChannels, std::vector<float>(numTicks));

  for (size_t i=0; i<numChannels; ++i) {
    SobelXFastRow(input2D[i], buffer[i]);
  }

  for (size_t j=0; j<numTicks; ++j) {
    std::vector<float> columnIn(numChannels);
    std::vector<float> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      columnIn[i] = buffer[i][j];
    }
    SobelXFastCol(columnIn, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      gradient[i][j] = columnOut[i];
    }
  }
  return;
}

void sigproc_tools::EdgeDetection::SobelXFastRow(
  const std::vector<float>& inputRow,
  std::vector<float>& outputRow) const
{
  const size_t N = inputRow.size();
  assert(outputRow.size() == N);

  outputRow[0] = -inputRow[1];
  outputRow[N-1] = inputRow[N-2];

  for (size_t i=1; i<N-1; ++i) {
    outputRow[i] = inputRow[i-1] - inputRow[i+1];
  }
  return;
}

void sigproc_tools::EdgeDetection::SobelXFastCol(
  const std::vector<float>& inputRow,
  std::vector<float>& outputRow) const
{
  const size_t N = inputRow.size();
  assert(outputRow.size() == N);

  outputRow[0] = 2.0 * inputRow[0] + inputRow[1];
  outputRow[N-1] = 2.0 * inputRow[N-1] + inputRow[N-2];

  for (size_t i=1; i<N-1; ++i) {
    outputRow[i] = inputRow[i-1] + inputRow[i+1] + 2.0 * inputRow[i];
  }
  return;
}

void sigproc_tools::EdgeDetection::SobelYFast(
  const Array2D<float> &input2D,
  Array2D<float> &gradient) const
{
  const size_t numChannels = input2D.size();
  const size_t numTicks = input2D.at(0).size();

  assert(gradient.size() == numChannels);
  assert(gradient.at(0).size() == numTicks);

  Array2D<float> buffer(numChannels, std::vector<float>(numTicks));

  for (size_t i=0; i<numChannels; ++i) {
    SobelXFastCol(input2D[i], buffer[i]);
  }

  for (size_t j=0; j<numTicks; ++j) {
    std::vector<float> columnIn(numChannels);
    std::vector<float> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      columnIn[i] = buffer[i][j];
    }
    SobelXFastRow(columnIn, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      gradient[i][j] = columnOut[i];
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::LSDGradX(
  const Array2D<float>& input2D,
  Array2D<float>& output2D) const
{
  // Input kernel must be normalized. 

  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  for (int i=0; i<numChannels-1; ++i) {
    for (int j=0; j<numTicks-1; ++j) {
      output2D[i][j] = (input2D[i+1][j] + input2D[i+1][j+1] - \
                        input2D[i][j] - input2D[i][j+1]) / 2.0;
    }
  }

  // Use Periodic Boundary
  for (int j=0; j<numTicks-1; ++j) {
    output2D[numChannels-1][j] = (input2D[0][j] + input2D[0][j+1] - \
                                  input2D[numChannels-1][j] - input2D[numChannels-1][j+1]) / 2.0;
  }
  for (int i=0; i<numChannels-1; ++i) {
    output2D[i][numTicks-1] = (input2D[i+1][numTicks-1] + input2D[i+1][0] - \
                      input2D[i][numTicks-1] - input2D[i][0]) / 2.0;
  }
  output2D[numChannels-1][numTicks-1] = (input2D[0][numTicks-1] + input2D[0][0] - \
                                         input2D[numChannels-1][numTicks-1] - input2D[numChannels-1][0]) / 2.0;
  return;
}

void sigproc_tools::EdgeDetection::LSDGradY(
  const Array2D<float>& input2D,
  Array2D<float>& output2D) const
{
  // Input kernel must be normalized. 

  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  for (int i=0; i<numChannels-1; ++i) {
    for (int j=0; j<numTicks-1; ++j) {
      output2D[i][j] = (input2D[i][j+1] + input2D[i+1][j+1] - \
                        input2D[i][j] - input2D[i+1][j]) / 2.0;
    }
  }

  // Use Periodic Boundary
  for (int j=0; j<numTicks-1; ++j) {
    output2D[numChannels-1][j] = (input2D[numChannels-1][j+1] + input2D[0][j+1] - \
                                  input2D[numChannels-1][j] - input2D[0][j]) / 2.0;
  }
  for (int i=0; i<numChannels-1; ++i) {
    output2D[i][numTicks-1] = (input2D[i][0] + input2D[i+1][0] - \
                      input2D[i][numTicks-1] - input2D[i+1][numTicks-1]) / 2.0;
  }
  output2D[numChannels-1][numTicks-1] = (input2D[numChannels-1][0] + input2D[0][0] - \
                                         input2D[numChannels-1][numTicks-1] - input2D[0][numTicks-1]) / 2.0;
  return;
}

void sigproc_tools::EdgeDetection::LSDGrad(
  const Array2D<float>& input2D,
  Array2D<float>& gradX,
  Array2D<float>& gradY,
  Array2D<float>& gradient,
  Array2D<float>& direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D.at(0).size();

  LSDGradX(input2D, gradX);
  LSDGradY(input2D, gradY);

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      float g = sqrt(gradX[i][j]*gradX[i][j] + gradY[i][j]*gradY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(gradX[i][j], -gradY[i][j]);
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::EdgeNMS(
  const Array2D<float>& gradient2D,
  const Array2D<float>& degrees2D,
  Array2D<float>& output2D
) const
{
  int numChannels = gradient2D.size();
  int numTicks = gradient2D.at(0).size();

  // Compensate for boundary
  for (int i=1; i<numChannels-1; ++i) {
    for (int j=1; j<numTicks-1; ++j) {
      if ((112.5 <= degrees2D[i][j] && degrees2D[i][j] < 157.5)
        || (-67.5 <= degrees2D[i][j] && degrees2D[i][j] < -22.5)) {
        // Upper Diag Left - Lower Diag Right
        float v1 = gradient2D[i-1][j+1];
        float v2 = gradient2D[i+1][j-1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else if ((22.5 <= degrees2D[i][j] && degrees2D[i][j] < 67.5) 
        || (-157.5 <= degrees2D[i][j] && degrees2D[i][j] < -112.5)) {
        // Upper Diag Right - Lower Diag Left
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i-1][j-1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else if ((67.5 <= degrees2D[i][j] && degrees2D[i][j] < 112.5)
        || (-112.5 <= degrees2D[i][j] && degrees2D[i][j] < -67.5)) {
        // Up-Down
        float v1 = gradient2D[i][j-1];
        float v2 = gradient2D[i][j+1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
      else{
        // Left-Right
        float v1 = gradient2D[i+1][j];
        float v2 = gradient2D[i-1][j];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2) {
          output2D[i][j] = gradient2D[i][j];
        } 
        else {
          output2D[i][j] = 0.0;
        }
      }
    }
  }

  return;
}


void sigproc_tools::EdgeDetection::EdgeNMSInterpolation(
  const Array2D<float>& gradient2D,
  const Array2D<float>& gradX,
  const Array2D<float>& gradY,
  const Array2D<float>& degrees2D,
  Array2D<float>& output2D) const
{
  int numChannels = gradient2D.size();
  int numTicks = gradient2D.at(0).size();

  // Implementation was partly adapted from:
  // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

  // Compensate for boundary
  for (int i=1; i<numChannels-1; ++i) {
    for (int j=1; j<numTicks-1; ++j) {
      if ((0 <= degrees2D[i][j] && degrees2D[i][j] < 45) 
        || (-180 <= degrees2D[i][j] && -135 > degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i+1][j];
        float v3 = gradient2D[i-1][j];
        float v4 = gradient2D[i-1][j-1];
        float G = abs(gradY[i][j] / gradX[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = gradient2D[i][j];
        } else {
          output2D[i][j] = 0;
        }
      }
      else if ((-45 <= degrees2D[i][j] && degrees2D[i][j] < 0) 
        || (135 <= degrees2D[i][j] && 180 >= degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j-1];
        float v2 = gradient2D[i+1][j];
        float v3 = gradient2D[i-1][j];
        float v4 = gradient2D[i-1][j+1];
        float G = abs(gradY[i][j] / gradX[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = gradient2D[i][j];
        } else {
          output2D[i][j] = 0;
        }
      }
      else if ((45 <= degrees2D[i][j] && degrees2D[i][j] < 90) 
        || (-135 <= degrees2D[i][j] && -90 > degrees2D[i][j])) {
        float v1 = gradient2D[i+1][j+1];
        float v2 = gradient2D[i][j+1];
        float v3 = gradient2D[i][j-1];
        float v4 = gradient2D[i-1][j-1];
        float G = abs(gradX[i][j] / gradY[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = gradient2D[i][j];
        } else {
          output2D[i][j] = 0;
        }
      }
      else {
        float v1 = gradient2D[i-1][j+1];
        float v2 = gradient2D[i][j+1];
        float v3 = gradient2D[i][j-1];
        float v4 = gradient2D[i+1][j-1];
        float G = abs(gradX[i][j] / gradY[i][j]);
        if ((gradient2D[i][j] >= v1 + (v2 - v1)*G) && 
            (gradient2D[i][j] >= v4 + G*(v3 - v4))) {
          output2D[i][j] = gradient2D[i][j];
        } else {
          output2D[i][j] = 0;
        }
      }
    }
  }

  return;
}


void sigproc_tools::EdgeDetection::DoubleThresholding(
  const Array2D<float>& doneNMS2D,
  Array2D<bool>& binary2D,
  std::vector<int>& strongEdgeRows,
  std::vector<int>& strongEdgeCols,
  std::vector<int>& weakEdgeRows,
  std::vector<int>& weakEdgeCols,
  float lowThreshold,
  float highThreshold) const
{
  int numChannels = doneNMS2D.size();
  int numTicks = doneNMS2D.at(0).size();

  // Implementation was partly adapted from:
  // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

  strongEdgeRows.reserve(numChannels * numTicks);
  strongEdgeCols.reserve(numChannels * numTicks);
  weakEdgeRows.reserve(numChannels * numTicks);
  weakEdgeCols.reserve(numChannels * numTicks);
  // Enumerator for strong / weak edges

  // Compensate for boundary
  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (doneNMS2D[i][j] >= highThreshold) {
        binary2D[i][j] = true;
        strongEdgeRows.push_back(i);
        strongEdgeCols.push_back(j);
      }
      else if (doneNMS2D[i][j] < highThreshold && doneNMS2D[i][j] >= lowThreshold) {
        weakEdgeRows.push_back(i);
        weakEdgeCols.push_back(j);
      }
      else {
        binary2D[i][j] = false;
      }
    }
  }
  return;
}

void sigproc_tools::EdgeDetection::HysteresisThresholdingFast(
  const Array2D<float>& doneNMS2D,
  float lowThreshold,
  float highThreshold,
  Array2D<bool>& outputROI) const
{
  const int numChannels = doneNMS2D.size();
  const int numTicks = doneNMS2D.at(0).size();

  const int forestSize = numChannels * numTicks;

  DisjointSetForest forest(forestSize, forestSize);

  forest.MakeSet();
  // std::cout << "size = " << forest.size << std::endl;

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      int flatIndex = i * numTicks + j;
      if (doneNMS2D[i][j] >= highThreshold) forest.parent[flatIndex] = forestSize;
    }
  }

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      int flatIndex = i * numTicks + j;
      int lowerBoundx = std::max(i-1, 0);
      int upperBoundx = std::min(i+2, (int) numChannels);
      int lowerBoundy = std::max(j-1, 0);
      int upperBoundy = std::min(j+2, (int) numTicks);
      // Process strong edges and its neighbors
      if (doneNMS2D[i][j] >= highThreshold) {
        if (forest.Find(flatIndex) == flatIndex) {
          forest.parent[flatIndex] = forestSize;
        }
        // std::cout << "Strong Edge: " << i << ", " << j << "    " << "Parent: " << forest.parent[flatIndex] << std::endl;
        for (int k=lowerBoundx; k<upperBoundx; ++k) {
          for (int l=lowerBoundy; l<upperBoundy; ++l) {
            int flatIndexNeigh = k * numTicks + l;
            const float &grad = doneNMS2D[k][l];
            if (grad >= lowThreshold) {
              // std::cout << "    Neighbor Index: " << k << ", " << l << std::endl;
              forest.Union(flatIndexNeigh, flatIndex);
            }
          }
        }
      }
      // Process weak edges
      else if ( (doneNMS2D[i][j] < highThreshold) && (doneNMS2D[i][j] >= lowThreshold)) {
        // std::cout << "Weak Edge: " << i << ", " << j << std::endl;
        for (int k=lowerBoundx; k<upperBoundx; ++k) {
          for (int l=lowerBoundy; l<upperBoundy; ++l) {
            int flatIndexNeigh = k * numTicks + l;
            const float &grad = doneNMS2D[k][l];
            // std::cout << "    Weak Edge: " << k << ", " << l << ", grad = " << grad << std::endl;
            if (grad >= lowThreshold) {
              forest.Union(flatIndexNeigh, flatIndex);
            }
          }
        }
      }
      else continue;
    }
  }
  // std::cout << "-----------------------------------------------" << std::endl;
  for (int flatIdx=0; flatIdx<forestSize; ++flatIdx) {
    int rep = forest.Find(flatIdx);
    int row = (flatIdx / numTicks);
    int col = (flatIdx % numTicks);
    if (rep == forestSize) outputROI[row][col] = true;
    else outputROI[row][col] = false;
  }
  return;
}


void sigproc_tools::EdgeDetection::HTFastLowMem(
  const Array2D<float>& doneNMS2D,
  float lowThreshold,
  float highThreshold,
  Array2D<bool>& outputROI) const
{
  const int numChannels = doneNMS2D.size();
  const int numTicks = doneNMS2D.at(0).size();

  // const int forestSize = numChannels * numTicks;

  // DisjointSetForest forest(forestSize);

  std::unordered_map<long, EdgeCandidate> edges;
  // std::vector<long> keyValues;

  // forest.MakeSet();
  // std::cout << "size = " << forest.size << std::endl;
  int count_id = 1;
  // int numWeakEdges = 0;
  // int numStrongEdges = 0;
  // First Pass
  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (doneNMS2D[i][j] >= highThreshold) {
        // if edge is a strong edge, assign to root node 
        // ID of 0 is reserved for root note reference
        EdgeCandidate strongEdge(i, j, 0, true);
        long key = CantorEnum(i, j);
        // keyValues.push_back(key);
        // count_id++;
        edges.emplace(std::make_pair(key, strongEdge));
        count_id++;
        // numStrongEdges++;
      }
      else if ( (doneNMS2D[i][j] < highThreshold) && 
                (doneNMS2D[i][j] >= lowThreshold)) {
        EdgeCandidate weakEdge(i, j, count_id, false);
        long key = CantorEnum(i, j);
        // keyValues.push_back(key);
        edges.emplace(std::make_pair(key, weakEdge));
        count_id++;
        // numWeakEdges++;
      }
      else continue;
    }
  }

  // std::cout << "Final Num IDs = " << count_id << std::endl;

  const int forestSize = edges.size();

  // std::cout << "forestSize = " << edges.size() << std::endl;

  // std::cout << "numStrongEdges = " << numStrongEdges << std::endl;

  // std::cout << "numWeakEdges = " << numWeakEdges << std::endl;

  // std::sort(keyValues.begin(), keyValues.end());
  // int uniqueCount = std::unique(keyValues.begin(), keyValues.end()) - keyValues.begin();

  // std::cout << "uniqueCount = " << uniqueCount << std::endl;

  DisjointSetForest forest(forestSize);
  forest.MakeSet();

  // Assign all strong edge to root node.
  for (auto& node : edges) {

    EdgeCandidate &edge = node.second;
    int i = edge.row;
    int j = edge.col;

    // std::cout << "Edge Index: i = " << i << ", j = " << j << std::endl;

    int lowerBoundx = std::max(i-1, 0);
    int upperBoundx = std::min(i+2, (int) numChannels);
    int lowerBoundy = std::max(j-1, 0);
    int upperBoundy = std::min(j+2, (int) numTicks);

    if (edge.edgeType) {
      // std::cout << "    Strong EdgeCandidate: row(" << edge.row << ") col(" << edge.col << ") type(" << edge.edgeType << ") id(" << edge.id << ")" << std::endl;
      if (forest.Find(edge.id) == edge.id) {
        // Assign strong edge to root node "0"
        forest.parent[edge.id] = 0;
      }
      // Handle neighbors
      for (int k=lowerBoundx; k<upperBoundx; ++k) {
        for (int l=lowerBoundy; l<upperBoundy; ++l) {
          const float &grad = doneNMS2D[k][l];
          if (grad >= lowThreshold) {
            long key = CantorEnum(k, l);
            // EdgeCandidate &neighborEdge = edges[key];
            // std::cout << "        Neighbor EdgeCandidate: row(" << edges[key].row << ") col(" << edges[key].col << ") type(" << edges[key].edgeType << ") id(" << edges[key].id << ")" << std::endl;
            forest.Union(edges[key].id, edge.id);
          }
        }
      }
    }
    // Process Weak Edges
    else {
      // std::cout << "    Weak EdgeCandidate: row(" << edge.row << ") col(" << edge.col << ") type(" << edge.edgeType << ") id(" << edge.id << ")" << std::endl;
      for (int k=lowerBoundx; k<upperBoundx; ++k) {
        for (int l=lowerBoundy; l<upperBoundy; ++l) {
          const float &grad = doneNMS2D[k][l];
          if (grad >= lowThreshold) {
            long key = CantorEnum(k, l);
            // std::cout << "        Neighbor EdgeCandidate: row(" << edges[key].row << ") col(" << edges[key].col << ") type(" << edges[key].edgeType << ") id(" << edges[key].id << ")" << std::endl;
            // EdgeCandidate &neighborEdge = edges[key];
            forest.Union(edges[key].id, edge.id);
          }
        }
      }
    }
  }

  for (auto& node: edges) {
    EdgeCandidate &edge = node.second;
    int rep = forest.Find(edge.id);
    const int &row = edge.row;
    const int &col = edge.col;

    // std::cout << "row(" << edge.row << ") col(" << edge.col << ") Rep = " << rep << ", EdgeType = " << edge.edgeType << std::endl;

    if (rep == 0) outputROI[row][col] = true;
  }

  return;
}


void sigproc_tools::EdgeDetection::HysteresisThresholding(
  const Array2D<bool>& binary2D,
  const std::vector<int>& strongEdgeRows,
  const std::vector<int>& strongEdgeCols,
  const std::vector<int>& weakEdgeRows,
  const std::vector<int>& weakEdgeCols,
  Array2D<bool>& output2D) const
{
  int numChannels = binary2D.size();
  int numTicks = binary2D.at(0).size();

  bool converged = false;

  sigproc_tools::Morph2DFast morph2d;

  // Initialize fixed point iteration buffer
  Array2D<bool> tempBuffer(numChannels);
  for (auto& v : tempBuffer) {
    v.resize(numTicks);
  }

  for (size_t i=0; i<strongEdgeRows.size(); ++i) {
    tempBuffer[strongEdgeRows[i]][strongEdgeCols[i]] = true;
    output2D[strongEdgeRows[i]][strongEdgeCols[i]] = true;
  }

  assert(strongEdgeRows.size() == strongEdgeCols.size());
  assert(weakEdgeRows.size() == weakEdgeCols.size());

  // Construct strong edge / weak edge 2D binary arrays from indices
  // REMARK: binary2D already contains strong edges. 
  // Array2D<bool> strongEdges2D(numChannels);
  // for (auto& v : strongEdge2D) {
  //   v.resize(numTicks);
  // }
  Array2D<bool> weakEdges2D(numChannels);
  for (auto& v : weakEdges2D) {
    v.resize(numTicks);
  }

  for (size_t i=0; i<weakEdgeRows.size(); ++i) {
    weakEdges2D[weakEdgeRows[i]][weakEdgeCols[i]] = true;
  }
  // This fixed point iteration always converges. 
  while (!converged) {
    // Dilation + compute overlap
    morph2d.getDilation(output2D, 3, 3, tempBuffer);
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        // Pixels have to be true in both dilated binary image and weak edges
        tempBuffer[i][j] = output2D[i][j] || (tempBuffer[i][j] && weakEdges2D[i][j]);
      }
    }
    converged = true;
    // Check array equality
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        if (tempBuffer[i][j] != output2D[i][j]) {
          converged = false;
          goto nextIter; // goto is used to avoid outer loop break
        }
      }
    }
    nextIter:
    // Update buffer
    for (int i=0; i<numChannels; ++i) {
      for (int j=0; j<numTicks; ++j) {
        output2D[i][j] = tempBuffer[i][j];
      }
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::Canny(
  const Array2D<float>& waveLessCoherent,
  Array2D<bool>& output2D,
  const unsigned int sx,
  const unsigned int sy,
  const float sigma_x,
  const float sigma_y,
  const float sigma_r,
  const float lowThreshold,
  const float highThreshold,
  const char mode) const
{
  /*
  This implementation of canny edge detection replaces gaussian smoothing with
  edge-preserving bilateral filtering.
  */
  int numChannels = waveLessCoherent.size();
  int numTicks = waveLessCoherent.at(0).size();

  sigproc_tools::BilateralFilters filter;

  sigproc_tools::Morph2DFast morph2d;

  Array2D<float> temp(numChannels);
  for (auto& v : temp) {
    v.resize(numTicks);
  }

  Array2D<bool> boolTemp(numChannels);
  for (auto& v : boolTemp) {
    v.resize(numTicks);
  }

  Array2D<float> gradient(numChannels);
  for (auto& v : gradient) {
    v.resize(numTicks);
  }

  Array2D<float> direction(numChannels);
  for (auto& v : direction) {
    v.resize(numTicks);
  }

  Array2D<float> sobelX(numChannels);
  for (auto& v : sobelX) {
    v.resize(numTicks);
  }

  Array2D<float> sobelY(numChannels);
  for (auto& v : sobelY) {
    v.resize(numTicks);
  }

  Array2D<float> morphed2D(numChannels);
  for (auto& v : morphed2D) {
    v.resize(numTicks);
  }


  Sobel(waveLessCoherent, sobelX, sobelY, gradient, direction);

  // 1. Perform Edge-Preserving Smoothing
  filter.directional(waveLessCoherent, direction, temp, sx, sy, sigma_x, sigma_y, sigma_r, 360);

  // 2. Run Morphological FIlter
  if (mode == 'e') {
    morph2d.getErosion(temp, sx, sy, morphed2D);
  }
  else if (mode == 'd') {
    morph2d.getDilation(temp, sx, sy, morphed2D);
  }
  else {
    morph2d.getGradient(temp, sx, sy, morphed2D);
  }

  // 2. Run Sobel Edge Filtering
  Sobel(morphed2D, sobelX, sobelY, gradient, direction);

  // 3. NMS Edge with Interpolation
  EdgeNMSInterpolation(gradient, sobelX, sobelY, direction, temp);

  // 4. Run double thresholding

  std::vector<int> strongEdgeRows;
  std::vector<int> strongEdgeCols;
  std::vector<int> weakEdgeRows;
  std::vector<int> weakEdgeCols;

  DoubleThresholding(temp, boolTemp, 
    strongEdgeRows, strongEdgeCols, 
    weakEdgeRows, weakEdgeCols, lowThreshold, highThreshold);

  HysteresisThresholding(boolTemp, strongEdgeRows, strongEdgeCols, weakEdgeRows, weakEdgeCols, output2D);

  return;
}


void sigproc_tools::EdgeDetection::gradientRegionGrow(
  const Array2D<float>& direction,
  const int anchorX,
  const int anchorY,
  const int regionID,
  const float tolerance,
  const unsigned int windowX, 
  const unsigned int windowY, 
  Array2D<int>& partitions) const
{
  assert(tolerance > 0);
  assert(regionID > 0);
  int numChannels = direction.size();
  int numTicks = direction.at(0).size();

  float Sx = 0.0;
  float Sy = 0.0;
  float meanAngle = direction[anchorX][anchorY];
  std::queue<int> anchorsX;
  std::queue<int> anchorsY;

  anchorsX.push(anchorX);
  anchorsY.push(anchorY);

  unsigned int count = 1;

  while (!anchorsX.empty()) {
    // std::cout << anchorsX.size() << std::endl;
    // std::cout << anchorsY.size() << std::endl;
    // Neighbors of current pixel
    int ix = anchorsX.front();
    int iy = anchorsY.front();
    anchorsX.pop();
    anchorsY.pop();

    int lbx = ix - (int) windowX;
    int ubx = ix + (int) windowX;
    int lby = iy - (int) windowY;
    int uby = iy + (int) windowY;
    int lowerBoundx = std::max(lbx, 0);
    int upperBoundx = std::min(ubx, (int) numChannels);
    int lowerBoundy = std::max(lby, 0);
    int upperBoundy = std::min(uby, (int) numTicks);

    for (int n=lowerBoundx; n<upperBoundx; ++n) {
      for (int m=lowerBoundy; m<upperBoundy; ++m) {
        float diffAngle = (float) abs(meanAngle - direction[n][m]);
        // std::cout << "Mean Angle = " << meanAngle << std::endl;
        // std::cout << "Diff Angle = " << diffAngle << std::endl;
        if (partitions[n][m] == 0 && diffAngle < tolerance) {
          anchorsX.push(n);
          anchorsY.push(m);
          count++;
          partitions[n][m] = regionID;
          Sx = Sx + cos(direction[n][m] * M_PI / 180.0);
          Sy = Sy + sin(direction[n][m] * M_PI / 180.0);
          meanAngle = atan2(Sy, Sx) * 180.0 / M_PI;
        }
      }
    }
  }
  return;
}


void sigproc_tools::EdgeDetection::regionGrow2D(
  const Array2D<float>& direction,
  const std::vector<int>& anchorsX,
  const std::vector<int>& anchorsY,
  const float tolerance,
  const unsigned int windowX, 
  const unsigned int windowY, 
  Array2D<int>& partitions
) const
{

  assert(anchorsX.size() == anchorsY.size());

  int groupID = 1;
  for (size_t i=0; i<anchorsX.size(); ++i) {
    // std::cout << "Label = " << groupID << std::endl;
    int ix = anchorsX[i];
    int iy = anchorsY[i];
    gradientRegionGrow(direction, ix, iy, groupID, tolerance, windowX, windowY, partitions);
    groupID++;
  }
  return;
}

#endif
