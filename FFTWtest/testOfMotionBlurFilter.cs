using ImageReadCS;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Drawing;
using System.Drawing.Imaging;



namespace RidgeDetection
{
    class testOfMotionBlurFilter
    {

        public static GrayscaleFloatImage rotateImage(GrayscaleFloatImage image, bool isCW, float angleDeg)
        {
            float angle = angleDeg * (float)Math.PI / 180;
            int x = 0;
            int y = 0;
            int iCentreX = image.Width / 2;
            int iCentreY = image.Height / 2;
            //int deltaPlus=500;
            int bilinearRotationWidth = (int)Math.Round(image.Width * Math.Abs(Math.Cos(angle)) + image.Height * Math.Abs(Math.Sin(angle)));//image.Width + deltaPlus;
            int bilinearRotationHeight = (int)Math.Round(image.Width * Math.Abs(Math.Sin(angle)) + image.Height * Math.Abs(Math.Cos(angle)));//image.Height + deltaPlus;
            GrayscaleFloatImage bilinearInterpolatationForRotation = new GrayscaleFloatImage(bilinearRotationWidth, bilinearRotationHeight);
            int radius = 0;//1;
            GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);

            if (angleDeg % 90 == 0)
            {

                if (angleDeg % 360 == 0)
                {
                    //ImageIO.ImageToFile(image, fileOutName);
                    return image;
                }
                else if (angleDeg % 180 == 0)
                {
                    for (int i = 0; i < image.Width; i++)
                    {
                        for (int j = 0; j < image.Height; j++)
                        {
                            bilinearInterpolatationForRotation[image.Width - i - 1, image.Height - j - 1] = image[i, j];
                        }
                    }
                    //ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                    return bilinearInterpolatationForRotation;
                }
                else if (angleDeg % 270 == 0)
                {
                    int isPlus270 = (int)angleDeg / 270;
                    if (((isPlus270 > 0) && (isPlus270 % 2 == 0)) || ((isPlus270 < 0) && (isPlus270 % 2 != 0)))//+90
                    {
                        for (int i = 0; i < image.Width; i++)
                        {
                            for (int j = 0; j < image.Height; j++)
                            {
                                bilinearInterpolatationForRotation[image.Height - j - 1, i] = image[i, j];
                            }
                        }

                    }
                    else
                    {
                        for (int i = 0; i < image.Width; i++)
                        {
                            for (int j = 0; j < image.Height; j++)
                            {
                                bilinearInterpolatationForRotation[j, image.Width - i - 1] = image[i, j];
                            }
                        }
                    }
                    //ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                    return bilinearInterpolatationForRotation;
                }

                //+90
                int isPlus90 = (int)angleDeg / 90;

                if (((isPlus90 > 0) && (isPlus90 % 2 != 0)) || ((isPlus90 < 0) && (isPlus90 % 2 == 0)))//+90
                {
                    for (int i = 0; i < image.Width; i++)
                    {
                        for (int j = 0; j < image.Height; j++)
                        {
                            bilinearInterpolatationForRotation[image.Height - j - 1, i] = image[i, j];
                        }
                    }

                }
                else
                {
                    for (int i = 0; i < image.Width; i++)
                    {
                        for (int j = 0; j < image.Height; j++)
                        {
                            bilinearInterpolatationForRotation[j, image.Width - i - 1] = image[i, j];
                        }
                    }
                }

                //ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                return bilinearInterpolatationForRotation;
            }

            //mirrow edges
            //copy center
            //копируем в середину
            for (int k = 0; k < image.Width; k++)
            {
                for (int l = 0; l < image.Height; l++)
                {
                    imageCopyExtended[k + radius, l + radius] = image[k, l];
                }
            }
            //копируем края
            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                }
            }

            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int k = 0; k < radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                }
            }

            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                }
            }
            //ImageIO.ImageToFile(imageCopyExtended, fileOutName);



            for (int i = 0; i < bilinearRotationHeight; i++)
            {
                for (int j = 0; j < bilinearRotationWidth; j++)
                {

                    // convert raster to Cartesian
                    x = j - (int)bilinearRotationWidth / 2;
                    y = (int)bilinearRotationHeight / 2 - i;

                    // convert Cartesian to polar
                    float fDistance = (float)Math.Sqrt(x * x + y * y);
                    float fPolarAngle = 0.0f;
                    if (x == 0)
                    {
                        if (y == 0)
                        {
                            // centre of image, no rotation needed
                            //bilinearInterpolatationForRotation[j, i] = image[x + image.Width / 2, y + image.Width / 2];
                            bilinearInterpolatationForRotation[j, i] = imageCopyExtended[x + image.Width / 2 + radius, y + image.Width / 2 + radius];
                            continue;
                        }
                        else if (y < 0)
                        {
                            fPolarAngle = (float)(1.5f * Math.PI);
                        }
                        else
                        {
                            fPolarAngle = (float)(0.5 * Math.PI);
                        }
                    }
                    else
                    {
                        fPolarAngle = (float)Math.Atan2((double)y, (double)x);
                    }

                    // the crucial rotation part
                    // "reverse" rotate, so minus instead of plus
                    fPolarAngle += angle;

                    // convert polar to Cartesian
                    float fTrueX = (float)(fDistance * Math.Cos(fPolarAngle));
                    float fTrueY = (float)(fDistance * Math.Sin(fPolarAngle));

                    // convert Cartesian to raster
                    fTrueX = (float)(fTrueX + (double)iCentreX);
                    fTrueY = (float)((double)iCentreY - fTrueY);

                    int iFloorX = (int)(Math.Floor(fTrueX));
                    int iFloorY = (int)(Math.Floor(fTrueY));
                    int iCeilingX = (int)(Math.Ceiling(fTrueX));
                    int iCeilingY = (int)(Math.Ceiling(fTrueY));

                    //check bounds
                    if (iFloorX < 0 || iCeilingX < 0 || iFloorX >= image.Width ||
                        iCeilingX >= image.Width || iFloorY < 0 || iCeilingY < 0 ||
                        iFloorY >= image.Height || iCeilingY >= image.Height) continue;

                    //double fDeltaX = fTrueX - (double)iFloorX;
                    //double fDeltaY = fTrueY - (double)iFloorY;

                    //ColorFloatPixel clrTopLeft = image[iFloorX, iFloorY];
                    //ColorFloatPixel clrTopRight = image[iCeilingX, iFloorY];
                    //ColorFloatPixel clrBottomLeft = image[iFloorX, iCeilingY];
                    //ColorFloatPixel clrBottomRight = image[iCeilingX, iCeilingY];

                    //check bounds
                    //if (iFloorX < -1 || iCeilingX < -1 || iFloorX >= (image.Width+1) ||
                    //    iCeilingX >= (image.Width+1) || iFloorY < -1 || iCeilingY < -1 ||
                    //    iFloorY >= (image.Height+1) || iCeilingY >= (image.Height+1)) continue;

                    //double fDeltaX = fTrueX - (double)iFloorX;
                    //double fDeltaY = fTrueY - (double)iFloorY;

                    //ColorFloatPixel clrTopLeft = imageCopyExtended[iFloorX+radius, iFloorY+radius];
                    //ColorFloatPixel clrTopRight = imageCopyExtended[iCeilingX + radius, iFloorY + radius];
                    //ColorFloatPixel clrBottomLeft = imageCopyExtended[iFloorX+radius, iCeilingY+radius];
                    //ColorFloatPixel clrBottomRight = imageCopyExtended[iCeilingX + radius, iCeilingY + radius];
                    // check bounds
                    //if (iFloorX < -2 && iCeilingX < -2 || iFloorX >= image.Width+2 ||
                    //   iCeilingX >= image.Width+2 || iFloorY < -2 || iCeilingY < -2 ||
                    //   iFloorY >= image.Height+2 || iCeilingY >= image.Height + 2) continue;


                    //if (iFloorX < 0) iFloorX = 0;
                    //if (iFloorY < 0) iFloorY = 0;
                    //if (iCeilingX < 0) iCeilingX = 0;
                    //if (iCeilingY < 0) iCeilingY = 0;
                    //if (iFloorX >= image.Width) iFloorX = image.Width-1;
                    //if (iFloorY >= image.Height) iFloorY = image.Height-1;
                    //if (iCeilingX >= image.Width) iCeilingX = image.Width-1;
                    //if (iCeilingY >= image.Height) iCeilingY = image.Height-1;



                    double fDeltaX = fTrueX - (double)iFloorX;
                    double fDeltaY = fTrueY - (double)iFloorY;

                    float clrTopLeft = image[iFloorX, iFloorY];
                    float clrTopRight = image[iCeilingX, iFloorY];
                    float clrBottomLeft = image[iFloorX, iCeilingY];
                    float clrBottomRight = image[iCeilingX, iCeilingY];

                    // linearly interpolate horizontally between top neighbours
                    float fTopRed = (float)((1 - fDeltaX) * clrTopLeft + fDeltaX * clrTopRight);
                    //float fTopGreen = (float)((1 - fDeltaX) * clrTopLeft.g + fDeltaX * clrTopRight.g);
                    //float fTopBlue = (float)((1 - fDeltaX) * clrTopLeft.b + fDeltaX * clrTopRight.b);

                    // linearly interpolate horizontally between bottom neighbours
                    float fBottomRed = (float)((1 - fDeltaX) * clrBottomLeft + fDeltaX * clrBottomRight);
                    //float fBottomGreen = (float)((1 - fDeltaX) * clrBottomLeft.g + fDeltaX * clrBottomRight.g);
                    //float fBottomBlue = (float)((1 - fDeltaX) * clrBottomLeft.b + fDeltaX * clrBottomRight.b);

                    // linearly interpolate vertically between top and bottom interpolated results
                    float iRed = (int)(Math.Round((1 - fDeltaY) * fTopRed + fDeltaY * fBottomRed));
                    //float iGreen = (int)(Math.Round((1 - fDeltaY) * fTopGreen + fDeltaY * fBottomGreen));
                    //float iBlue = (int)(Math.Round((1 - fDeltaY) * fTopBlue + fDeltaY * fBottomBlue));

                    // make sure colour values are valid
                    if (iRed < 0) iRed = 0;
                    if (iRed > 255) iRed = 255;
                    //if (iGreen < 0) iGreen = 0;
                    //if (iGreen > 255) iGreen = 255;
                    //if (iBlue < 0) iBlue = 0;
                    //if (iBlue > 255) iBlue = 255;

                    float newPixelInterpolated = iRed;//new ColorFloatPixel();
                    //newPixelInterpolated.r = iRed;
                    //newPixelInterpolated.g = iGreen;
                    //newPixelInterpolated.b = iBlue;

                    bilinearInterpolatationForRotation[j, i] = newPixelInterpolated;



                }
            }
            return bilinearInterpolatationForRotation;
            
        }

        static ColorFloatImage convolution2d(ColorFloatImage image, GrayscaleFloatImage filter)
        {
            int radius_h=0;
            if (filter.Height%2==0)
            {
                radius_h = filter.Height / 2 - 1;
            }
            else
            {
                radius_h = filter.Height / 2;
            }
            int radius_w = 0;
            if (filter.Width % 2 == 0)
            {
                radius_w = filter.Width / 2 - 1;
            }
            else
            {
                radius_w = filter.Width / 2;
            }
            int radius = Math.Min(radius_h, radius_w);
            float normSum = 0.0f;
            for (int i=0;i<filter.Width;i++)
            {
                for (int j=0;j<filter.Height;j++)
                {
                    normSum += filter[i, j];
                }
            }

            ColorFloatImage imageCopyExtended = new ColorFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            ColorFloatImage imageCopy = new ColorFloatImage(image.Width, image.Height);
            for (int i = radius; i < imageCopyExtended.Width - radius; i++)
            {
                for (int j = radius; j < imageCopyExtended.Height - radius; j++)
                {
                    imageCopyExtended[i, j] = image[i-radius, j-radius];
                }
            }
            ImageIO.ImageToFile(imageCopyExtended, @"Input\filterTest\imageCopyExtended.bmp");
            ImageIO.ImageToFile(filter, @"Input\filterTest\convFilter.bmp");
            for (int i = 0; i < imageCopy.Width; i++)
            {
                for (int j = 0; j < imageCopy.Height; j++)
                {
                    ColorFloatPixel p = new ColorFloatPixel();
                    float sumR = 0;
                    float sumG = 0;
                    float sumB = 0;
                    for (int x = -radius_w; x < radius_w; x++)
                    {
                        for (int y = -radius_h; y < radius_h; y++)
                        {
                            sumR += imageCopyExtended[(i + radius) + x, (j + radius) + y].r / normSum * filter[x + radius_w, y + radius_h];
                            sumG += imageCopyExtended[(i + radius) + x, (j + radius) + y].g / normSum * filter[x + radius_w, y + radius_h];
                            sumB += imageCopyExtended[(i + radius) + x, (j + radius) + y].b / normSum * filter[x + radius_w, y + radius_h];
                         
                        }
                    }
                    p.r = sumR;
                    p.g = sumG;
                    p.b = sumB;
                    imageCopy[i, j] = p;
                }
            }
            return imageCopy;
        }

        static GrayscaleFloatImage convolution2dGray(GrayscaleFloatImage image, GrayscaleFloatImage filter)
        {
            int radius_h = 0;
            if (filter.Height % 2 == 0)
            {
                radius_h = filter.Height / 2 - 1;
            }
            else
            {
                radius_h = filter.Height / 2;
            }
            int radius_w = 0;
            if (filter.Width % 2 == 0)
            {
                radius_w = filter.Width / 2 - 1;
            }
            else
            {
                radius_w = filter.Width / 2;
            }
            int radius = Math.Min(radius_h, radius_w);
            float normSum = 0.0f;
            for (int i = 0; i < filter.Width; i++)
            {
                for (int j = 0; j < filter.Height; j++)
                {
                    normSum += filter[i, j];
                }
            }

            GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopy = new GrayscaleFloatImage(image.Width, image.Height);
            for (int i = radius; i < imageCopyExtended.Width - radius; i++)
            {
                for (int j = radius; j < imageCopyExtended.Height - radius; j++)
                {
                    imageCopyExtended[i, j] = image[i - radius, j - radius];
                }
            }
            ImageIO.ImageToFile(imageCopyExtended, @"Input\filterTest\imageCopyExtended.bmp");
            ImageIO.ImageToFile(filter, @"Input\filterTest\convFilter.bmp");
            for (int i = 0; i < imageCopy.Width; i++)
            {
                for (int j = 0; j < imageCopy.Height; j++)
                {
                    float sumGray = 0;
                    for (int x = -radius_w; x < radius_w; x++)
                    {
                        for (int y = -radius_h; y < radius_h; y++)
                        {
                            sumGray += imageCopyExtended[(i + radius) + x, (j + radius) + y] / normSum * filter[x + radius_w, y + radius_h];
                        }
                    }
                    imageCopy[i, j] = sumGray;
                }
            }
            return imageCopy;
        }


        static void mothionBlur(ColorFloatImage image, float L, float angle)
        {
            int radius = (int)Math.Round(L);
            GrayscaleFloatImage imageCopyExtended_xx = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopy_xx = new GrayscaleFloatImage(image.Width, image.Height);

            GrayscaleFloatImage filterKernel = new GrayscaleFloatImage(2 * radius + 1, 2 * radius + 1);
            float sum = 2 * radius + 1;
            float real_sum = 0;
            for (int i=-radius;i<=radius;i++)
            {
                filterKernel[i + radius, radius] = 255;
                //real_sum += filterKernel[i + radius, radius];
            }
            ImageIO.ImageToFile(filterKernel, @"Input\filterTest\filterMotion.bmp");
            //for (int i = -radius; i <= radius; i++)
            //{
            //    filterKernel[i + radius, radius] /= real_sum; //sum;
            //}

            GrayscaleFloatImage filterMotionLine = rotateImage(filterKernel, true, angle);
            ImageIO.ImageToFile(convolution2d(image, filterMotionLine), @"Input\filterTest\MotionBlured.bmp");

        }
        static void mothionBlurGray(GrayscaleFloatImage image, float L, float angle)
        {
            int radius = (int)Math.Round(L);
            GrayscaleFloatImage imageCopyExtended_xx = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopy_xx = new GrayscaleFloatImage(image.Width, image.Height);

            GrayscaleFloatImage filterKernel = new GrayscaleFloatImage(2 * radius + 1, 2 * radius + 1);
            float sum = 2 * radius + 1;
            float real_sum = 0;
            for (int i = -radius; i <= radius; i++)
            {
                filterKernel[i + radius, radius] = 255;
                //real_sum += filterKernel[i + radius, radius];
            }
            ImageIO.ImageToFile(filterKernel, @"Input\filterTest\filterMotion.bmp");
            //for (int i = -radius; i <= radius; i++)
            //{
            //    filterKernel[i + radius, radius] /= real_sum; //sum;
            //}

            GrayscaleFloatImage filterMotionLine = rotateImage(filterKernel, true, angle);
            ImageIO.ImageToFile(convolution2dGray(image, filterMotionLine), @"Input\filterTest\MotionBluredGray.bmp");

        }
        static double t(int u, int v, float a, float b)
        {
            return Math.PI * (u * a + v * b);
        }

  

        //static void Main(string[] args)
        //{
        //    if (args.Length < 2)
        //        return;
        //    string InputFileName = args[0], OutputFileName = args[1];
        //    if (!File.Exists(InputFileName))
        //        return;

        //    ColorFloatImage image = ImageIO.FileToColorFloatImage(InputFileName);
        //    GrayscaleFloatImage imageGray = ImageIO.FileToGrayscaleFloatImage(InputFileName);

        //    //mothionBlur(image, 10, 45);
        //    //mothionBlurGray(imageGray, 10, 45);

            
        //    //MedianImage(image, 10);
        //    float DEG2RAD = (float)Math.PI / 180;
        //    ImageIO.ImageToFile(image, OutputFileName);
        //    Console.WriteLine("Done:)");
        //    Console.WriteLine(Math.Tan(45 * DEG2RAD));
        //    Console.ReadLine();

        //}

    }
}
