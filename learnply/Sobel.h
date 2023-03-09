#pragma once
#include <string>
#include <vector>

#include "Polyline.h"
#include "icMatrix.h"
#include "glError.h"
#include "polyhedron.h"

void initSobel();
void displaySobel();
void sobelFilter(const std::string& fname);

void createEdgeFieldFromSobel();
void drawstreamlines();

void initImage();
void displayImage();
void imageFilter(const std::string& fname);

icVector3 getVector(Quad* q, const icVector3& p);