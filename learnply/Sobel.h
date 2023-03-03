#pragma once
#include <string>

void initSobel();
void displaySobel();
void sobelFilter(const std::string& fname);

void createEdgeFieldFromSobel();

void initImage();
void displayImage();
void imageFilter(const std::string& fname);