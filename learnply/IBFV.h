#pragma once
#include <string>


void initIBFV();
void display_IBFV();
void displayIBFV();
void displaySobel();
void makePatternsImg(const std::string& fname);
void makePatternsImgNoise(const std::string& fname, float w);
void sobelFilter(const std::string& fname);