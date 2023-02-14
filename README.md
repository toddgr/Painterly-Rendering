# Thesis

Current progress on Honors undergraduate thesis project. Aiming to create an application that allows a "painterly rendering" of an image or video as input.

We present a tool that converts any image or video into a painterly rendered one, or one that looks as though it was hand-painted. To do this, we implemented the Sobel filter to detect the directional fields that are used to create the direction of the brush strokes. The system allows the user to adjust a few stylistic parameters regarding the brush strokes (width, color, texture, etc.) to provide the user with additional control over the rendered output. These brush strokes are then composited to create the final painterly rendering. This system will be valuated by a handful of users, and their feedback will ultimately be used to simplify the user interface.

## Sobel Filter

The first aspect of this project is the Sobel filter, which will be used to generate the streamlines used for brush strokes. The Sobel filter takes an image as input and creates a new image with emphasized edges. In this process, the gradient vector field of the input image is approximately computed in order to determine the intensity for each pixel. The image is then convolved with the filter operator in both the horizontal and vertical directions. In terms of image processing, convolution refers to a 3x3 matrix of surrounding pixels, which are used to compare the change of intensity in the selected pixel in any given direction. This intensity change is used to compute the gradient vector for that pixel, which points in the direction of the greatest change in intensity. For example, these images from Spongebob will present the following results when the Sobel operator is applied:

![Spongebob edge fields](spongebobresults.JPG)

## Streamline Creation

This part has not been implemented yet.

## Brushstroke Creation & stylization

This part has not been implemented yet.

Source code provided by mentor Dr. Eugene Zhang and Jinta Jeng
