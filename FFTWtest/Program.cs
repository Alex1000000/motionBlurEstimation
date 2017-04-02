using FFTWtest;
using ImageReadCS;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace RidgeDetection
{
    class Program
    {


        static double GaussFunction(float x, float x0, float sigma)//(float sigma, float i)
        {
            //double result = (1 / (Math.Sqrt(2 * Math.PI) * sigma)) * Math.Exp((-i * i) / (2 * sigma * sigma));
            //////return (float)result;
            //return result;
            ////            return Math.Exp (
            ////-(Math.Pow (i, 2)) / (2 * Math.Pow (sigma, 2)));
            return Math.Exp(-Math.Pow(x - x0, 2) / (2 * Math.Pow(sigma, 2)));
        }
        static double GaussFunctionDerivative(float x, float x0, float sigma)//float sigma, float i)
        {
            //////double result = -(1 / (Math.Sqrt(2 * Math.PI) * Math.Pow( sigma,3)))* i * Math.Exp( (-i * i) / (2 * sigma * sigma));
            //////return (float)result;
            //double result = -(1 / (i * Math.Pow(sigma, 2))) *GaussFunction(sigma, i);
            //return result;
            ///////////////return GaussFunction(i,i ,sigma) * (i) / Math.Pow(sigma, 2);//минус был
           // return Math.Exp(-Math.Pow(x - x0, 2) / (2 * Math.Pow(sigma, 2)));
            return -GaussFunction(x, x0, sigma) * (x - x0) / Math.Pow(sigma, 2);
        }
        static double GaussFunctionSecondDerivative(float x, float x0, float sigma)//(float sigma, float i)
        {
            //////double result = -(1 / (Math.Sqrt(2 * Math.PI) * Math.Pow( sigma,3)))* i * Math.Exp( (-i * i) / (2 * sigma * sigma));
            //////return (float)result;
            //double result = -(1 / (i * Math.Pow(sigma, 2))) *GaussFunction(sigma, i);
            //return result;
            //////////////return GaussFunction(sigma, i, i) * (2 * i * i / (sigma * sigma) - 2) / Math.Pow(sigma, 2);//минус был
            return GaussFunction(x, x0, sigma) *
                (Math.Pow(x0, 2) - Math.Pow(sigma, 2) + Math.Pow(x, 2) - 2 * x0 * x) / Math.Pow(sigma, 4);
        }

        static public double  Drop2PI(double phi) {
	        double times = Math.Floor(phi / (2 * Math.PI));
	        return phi - times * 2 * Math.PI;
        }

        static void RidgeDetection( GrayscaleFloatImage image, float sigma, int radius)
        {
            GrayscaleFloatImage imageCopyExtended_xx = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopyExtended_xy = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopyExtended_yy = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
            GrayscaleFloatImage imageCopy_xx = new GrayscaleFloatImage(image.Width , image.Height );
            GrayscaleFloatImage imageCopy_xy = new GrayscaleFloatImage(image.Width , image.Height );
            GrayscaleFloatImage imageCopy_yy = new GrayscaleFloatImage(image.Width , image.Height );
            //mirrow edges
            //copy center
            //копируем в середину
            for (int k = 0; k < image.Width; k++)
            {
                for (int l = 0; l < image.Height; l++)
                {
                    imageCopyExtended_xx[k + radius, l + radius] = image[k, l];
                    imageCopyExtended_xy[k + radius, l + radius] = image[k, l];
                    imageCopyExtended_yy[k + radius, l + radius] = image[k, l];
                }
            }
            //копируем края
            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[k, radius - l - 1] = imageCopyExtended_xx[k, l + radius];
                    imageCopyExtended_xy[k, radius - l - 1] = imageCopyExtended_xy[k, l + radius];
                    imageCopyExtended_yy[k, radius - l - 1] = imageCopyExtended_yy[k, l + radius];

                }
            }

            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k, l - radius];
                    imageCopyExtended_xy[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k, l - radius];
                    imageCopyExtended_yy[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k, l - radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, l] = imageCopyExtended_xx[k + radius, l];
                    imageCopyExtended_xy[radius - k - 1, l] = imageCopyExtended_xy[k + radius, l];
                    imageCopyExtended_yy[radius - k - 1, l] = imageCopyExtended_yy[k + radius, l];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_xx[k - radius, l];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_xy[k - radius, l];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_yy[k - radius, l];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int k = 0; k < radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, radius - l - 1] = imageCopyExtended_xx[k + radius, l + radius];
                    imageCopyExtended_xy[radius - k - 1, radius - l - 1] = imageCopyExtended_xy[k + radius, l + radius];
                    imageCopyExtended_yy[radius - k - 1, radius - l - 1] = imageCopyExtended_yy[k + radius, l + radius];

                }
            }

            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k - radius - 1, l - radius - 1];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k - radius - 1, l - radius - 1];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k - radius - 1, l - radius - 1];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_xx[k - radius, l + radius];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_xy[k - radius, l + radius];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_yy[k - radius, l + radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k + radius, l - radius];
                    imageCopyExtended_xy[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k + radius, l - radius];
                    imageCopyExtended_yy[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k + radius, l - radius];

                }
            }
            ///////////////////////////////////////////////////////////////////////
            float[] Kernel = new float[radius * 2 + 1];
            float SumKernel = 0;
            for (int k = -radius; k < radius + 1; k++)
            {

                Kernel[k + radius] = (float)GaussFunction(k+radius,radius, sigma);


                SumKernel += Kernel[k + radius];
                //Console.WriteLine(Kernel[k + radius] + "|");
            }
            for (int k = -radius; k < radius + 1; k++)
            {
                //Kernel[k + radius] /= SumKernel;

                //Console.WriteLine(Kernel[k + radius] + "|");
            }
            ////////////////////////////////////////////////////////////////////
            float[] KernelDerivative = new float[radius * 2 + 1];
            //Console.WriteLine("GausKernel_without norm");
            for (int k = -radius; k < radius + 1; k++)
            {

                KernelDerivative[k + radius] = (float)GaussFunctionDerivative(k + radius, radius, sigma);
                //SumKernelDerivative += KernelDerivative[k + radius];
                //Console.WriteLine(Kernel[k + radius] + "|");
            }
            
            ///////////////////////////////////////////////////////////////////
            float[] KernelSecondDerivative = new float[radius * 2 + 1];
            float SumKernelSecondDerivative = 0;
            for (int k = -radius; k < radius + 1; k++)
            {

                KernelSecondDerivative[k + radius] = (float)GaussFunctionSecondDerivative(k + radius, radius, sigma);

                //SumKernelSecondDerivative += KernelSecondDerivative[k + radius];
                //Console.WriteLine(Kernel[k + radius] + "|");
            }
            //for (int k = -radius; k < radius + 1; k++)
            //{
            //    KernelSecondDerivative[k + radius] /= SumKernelSecondDerivative;

            //    //Console.WriteLine(Kernel[k + radius] + "|");
            //}
            ///////////////////////////////////////////////////////////////////
            for (int y = radius; y < imageCopyExtended_xx.Height-radius; y++)
            {
                for (int x = radius; x < imageCopyExtended_xx.Width-radius; x++)
                {
                    float[,] hessian = new float[2, 2];
                    //hessian[0,0]=   
                    float sum_xx = 0;
                    float sum_xy = 0;
                    float sum_yy = 0;
                    for (int k = -radius; k < radius + 1; k++)
                    {
                        sum_xx += imageCopyExtended_xx[x + k, y] * KernelSecondDerivative[k + radius];
                        sum_xy += imageCopyExtended_xy[x + k, y] * KernelDerivative[k + radius];
                        sum_yy += imageCopyExtended_yy[x + k, y] * Kernel[k + radius];
                    }
                       //gyy
                    //hessian[1,0]=  //gyx
                    //hessian[0,1]=  //gxy
                    //hessian[1,1]=  //gxx
                    //derivativeY[x - radius, y - radius] = sum; //*100
                    imageCopy_xx[x-radius, y-radius] = sum_xx;
                    imageCopy_xy[x-radius, y-radius] = sum_xy;
                    imageCopy_yy[x-radius, y-radius] = sum_yy;

                }
            }
            ////////////////////////////////////////////////////////////////////////
            //ImageIO.ImageToFile(imageCopy_xx, "Input/xx.bmp");
            //ImageIO.ImageToFile(imageCopy_xy, "Input/xy.bmp");
            //ImageIO.ImageToFile(imageCopy_yy, "Input/yy.bmp");
            ////////////////////////////////////////////////////////////////////////
            //копируем в середину
            for (int k = 0; k < image.Width; k++)
            {
                for (int l = 0; l < image.Height; l++)
                {
                    imageCopyExtended_xx[k + radius, l + radius] = imageCopy_xx[k, l];
                    imageCopyExtended_xy[k + radius, l + radius] = imageCopy_xy[k, l];
                    imageCopyExtended_yy[k + radius, l + radius] = imageCopy_yy[k, l];
                }
            }
            //копируем края
            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[k, radius - l - 1] = imageCopyExtended_xx[k, l + radius];
                    imageCopyExtended_xy[k, radius - l - 1] = imageCopyExtended_xy[k, l + radius];
                    imageCopyExtended_yy[k, radius - l - 1] = imageCopyExtended_yy[k, l + radius];

                }
            }

            for (int k = radius; k < image.Width + radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k, l - radius];
                    imageCopyExtended_xy[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k, l - radius];
                    imageCopyExtended_yy[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k, l - radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, l] = imageCopyExtended_xx[k + radius, l];
                    imageCopyExtended_xy[radius - k - 1, l] = imageCopyExtended_xy[k + radius, l];
                    imageCopyExtended_yy[radius - k - 1, l] = imageCopyExtended_yy[k + radius, l];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = radius; l < image.Height + radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_xx[k - radius, l];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_xy[k - radius, l];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended_yy[k - radius, l];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int k = 0; k < radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, radius - l - 1] = imageCopyExtended_xx[k + radius, l + radius];
                    imageCopyExtended_xy[radius - k - 1, radius - l - 1] = imageCopyExtended_xy[k + radius, l + radius];
                    imageCopyExtended_yy[radius - k - 1, radius - l - 1] = imageCopyExtended_yy[k + radius, l + radius];

                }
            }

            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k - radius - 1, l - radius - 1];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k - radius - 1, l - radius - 1];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k - radius - 1, l - radius - 1];

                }
            }
            for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
            {
                for (int l = 0; l < radius; l++)
                {

                    imageCopyExtended_xx[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_xx[k - radius, l + radius];
                    imageCopyExtended_xy[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_xy[k - radius, l + radius];
                    imageCopyExtended_yy[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended_yy[k - radius, l + radius];

                }
            }
            for (int k = 0; k < radius; k++)
            {
                for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                {

                    imageCopyExtended_xx[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xx[k + radius, l - radius];
                    imageCopyExtended_xy[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_xy[k + radius, l - radius];
                    imageCopyExtended_yy[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended_yy[k + radius, l - radius];

                }
            }



            ////////////////////////////////////////////////////////////////////////
            for (int y = radius; y < imageCopyExtended_xx.Height-radius; y++)
            {
                for (int x = radius; x < imageCopyExtended_xx.Width-radius; x++)
                {
                    float[,] hessian = new float[2, 2];
                    //hessian[0,0]=   
                    float sum_xx = 0;
                    float sum_xy = 0;
                    float sum_yy = 0;
                    for (int k = -radius; k < radius + 1; k++)
                    {
                        sum_xx += imageCopyExtended_xx[x , y + k] * Kernel[k + radius];
                        sum_xy += imageCopyExtended_xy[x , y + k] * KernelDerivative[k + radius];
                        sum_yy += imageCopyExtended_yy[x , y + k ] * KernelSecondDerivative[k + radius];
                    }
                    //gyy
                    //hessian[1,0]=  //gyx
                    //hessian[0,1]=  //gxy
                    //hessian[1,1]=  //gxx
                    //derivativeY[x - radius, y - radius] = sum; //*100
                    imageCopy_xx[x-radius, y-radius] = sum_xx;
                    imageCopy_xy[x-radius, y-radius] = sum_xy;
                    imageCopy_yy[x-radius, y-radius] = sum_yy;

                }
            }
            //ImageIO.ImageToFile(imageCopy_xx, "Input/Output/xx.bmp");
            //ImageIO.ImageToFile(imageCopy_xy, "Input/Output/xy.bmp");
            //ImageIO.ImageToFile(imageCopy_yy, "Input/Output/yy.bmp");
            System.IO.StreamWriter file = new System.IO.StreamWriter("Input/angleAnalisys.txt");
            GrayscaleFloatImage directionMap = new GrayscaleFloatImage(image.Width, image.Height);
            double RAD2DEG = (180 / Math.PI);
            int[] angleMaxHist = new int[181];
            int[] angleMaxHistNonMax = new int[181];
            ColorFloatImage colorRidge = new ColorFloatImage(image.Width, image.Height);
            for (int y = 0; y < imageCopy_xx.Height; y++)
            {
                for (int x = 0; x < imageCopy_xx.Width; x++)
                {

                    float D = ((float)Math.Pow(imageCopy_xx[x, y] + imageCopy_yy[x, y],2) - 4 * (imageCopy_xx[x, y] * imageCopy_yy[x, y] - imageCopy_xy[x, y] * imageCopy_xy[x, y]));
                    if ((D<0)&&(Math.Abs(D-0.0001)<0.1))
                    {
                        D = 0;
                    }
                    float Lambda1 = (imageCopy_xx[x, y] + imageCopy_yy[x, y] + (float)Math.Sqrt(D)) / 2;
                    float Lambda2 = (imageCopy_xx[x, y] + imageCopy_yy[x, y] - (float)Math.Sqrt(D)) / 2;

                    float maxLambda = Math.Max(Lambda1, Lambda2);
                    float minLambda = Math.Min(Lambda1, Lambda2);
                    image[x, y] = maxLambda;


                    directionMap[x, y] = (float)(RAD2DEG* Math.Atan2((double)(imageCopy_xy[x, y]), (double)(maxLambda - imageCopy_xx[x, y])));
                    //if (directionMap[x, y] != 0)
                    {
                        if (Double.IsNaN(directionMap[x, y]))
                        {
                            Console.WriteLine("NaN");
                        }
                        double theta = directionMap[x, y];
                        if (theta < 0)
                        {
                            theta = theta + 180;//Math.PI;
                        }
                        //else
                        //{
                        int angle = (int)Math.Round(theta);

                        //graphGray[angle, 1]++;
                        if ((angle != 0) && (angle <= 180) && (angle > 0))
                        {
                            angleMaxHist[angle]++;
                            file.Write(angle + ", ");
                        }
                        //}

                    }
                }
            }

            //////
            file.Close();
            //угол всх точек
            int maxInd = 0;
            int max_angle = angleMaxHist[0];
            for (int i = 0; i < angleMaxHist.Length; i++)
            {
                if (max_angle < angleMaxHist[i])
                {
                    max_angle = angleMaxHist[i];
                    maxInd = i;
                }
            }
            Console.WriteLine("MAX_Angle=" + (maxInd - 90));
 
            //пока не будем портить т к нон макс
            float max = image[0, 0];
            float min = image[0, 0];
            for (int j = 0; j < image.Height; j++)
            {
                for (int i = 0; i < image.Width; i++)
                {
                    if (image[i, j] > max) { max = image[i, j]; }
                    if (image[i, j] < min) { min = image[i, j]; }
                }
            }

            for (int j = 0; j < image.Height; j++)
            {
                for (int i = 0; i < image.Width; i++)
                {
                    //где проверка на не деление на ноль!?
                    image[i, j] = ((image[i, j]) / (max));//для макс суп уберем пока *256;
                    ColorFloatPixel pixel = new ColorFloatPixel();
                    pixel.g = image[i, j];
                    pixel.r = 255;
                    pixel.b = 255;
                    colorRidge[i, j] = pixel;
                }
            }


            /////////////////////////////////////////////

            max = imageCopy_xx[0, 0];
            min = imageCopy_xx[0, 0];
            for (int j = 0; j < imageCopy_xx.Height; j++)
            {
                for (int i = 0; i < imageCopy_xx.Width; i++)
                {
                    if (imageCopy_xx[i, j] > max) { max = imageCopy_xx[i, j]; }
                    if (imageCopy_xx[i, j] < min) { min = imageCopy_xx[i, j]; }
                }
            }

            for (int j = 0; j < imageCopy_xx.Height; j++)
            {
                for (int i = 0; i < imageCopy_xx.Width; i++)
                {
                    imageCopy_xx[i, j] = ((imageCopy_xx[i, j] ) / (max)) * 256;

                }
            }

            max = imageCopy_xy[0, 0];
            min = imageCopy_xy[0, 0];
            for (int j = 0; j < imageCopy_xy.Height; j++)
            {
                for (int i = 0; i < imageCopy_xy.Width; i++)
                {
                    if (imageCopy_xy[i, j] > max) { max = imageCopy_xy[i, j]; }
                    if (imageCopy_xy[i, j] < min) { min = imageCopy_xy[i, j]; }
                }
            }

            for (int j = 0; j < imageCopy_xy.Height; j++)
            {
                for (int i = 0; i < imageCopy_xy.Width; i++)
                {
                    imageCopy_xy[i, j] = ((imageCopy_xy[i, j] ) / (max)) * 256;

                }
            }

            max = imageCopy_yy[0, 0];
            min = imageCopy_yy[0, 0];
            for (int j = 0; j < imageCopy_yy.Height; j++)
            {
                for (int i = 0; i < imageCopy_yy.Width; i++)
                {
                    if (imageCopy_yy[i, j] > max) { max = imageCopy_yy[i, j]; }
                    if (imageCopy_yy[i, j] < min) { min = imageCopy_yy[i, j]; }
                }
            }

            for (int j = 0; j < imageCopy_yy.Height; j++)
            {
                for (int i = 0; i < imageCopy_yy.Width; i++)
                {
                    imageCopy_yy[i, j] = ((imageCopy_yy[i, j] ) / (max)) * 256;

                }
            }

            ////////////////////////////////////////////
            ImageIO.ImageToFile(image, "Input/Output/image_ridge.bmp");
            //ImageIO.ImageToFile(colorRidge, "Input/Output/image_ClorRidge.bmp");

            ImageIO.ImageToFile(imageCopy_xx, "Input/Output/xx.bmp");
            ImageIO.ImageToFile(imageCopy_xy, "Input/Output/xy.bmp");
            ImageIO.ImageToFile(imageCopy_yy, "Input/Output/yy.bmp");
            /////////////////////////////////////////////////
            //ну а теперь нон максимум сапрешшн ...вообще хорошо бы его в отдельную фнкцию
            
            ColorFloatImage imageNonMax = new ColorFloatImage(image.Width , image.Height );
            GrayscaleFloatImage imageNonMax_gray = new GrayscaleFloatImage(image.Width , image.Height );
            System.IO.StreamWriter file_NonMax = new System.IO.StreamWriter("Input/angleAnalisysFromomMax.txt");

            for (int i = 1; i < image.Width - 1; i++)
            {
                for (int j = 1; j < image.Height - 1; j++)
                {
                    float imageNonMaxValue = image[i, j];
                    int angle_non_max = (int)(Math.Round(RAD2DEG * (Drop2PI(directionMap[i, j] - Math.PI / 2)) / 45)) * 45 % 180;
                    double[] compareTo = new double[2];


                    if (angle_non_max == 0)
                    {
                        compareTo[0] = image[i - 1, j];
                        compareTo[1] = image[i + 1, j];
                    }
                    else if (angle_non_max == 45)
                    {
                        compareTo[0] = image[i + 1, j - 1];
                        compareTo[1] = image[i - 1, j + 1];
                    }
                    else if (angle_non_max == 90)
                    {
                        compareTo[0] = image[i, j - 1];
                        compareTo[1] = image[i, j + 1];
                    }
                    else if (angle_non_max == 135)
                    {
                        compareTo[0] = image[i + 1, j + 1];
                        compareTo[1] = image[i - 1, j - 1];
                    }

                    ///это можно убрать, тогда будут жиирные линии
                    //if ((Math.Abs(imageNonMaxValue) <= Math.Abs(compareTo[0])) || (Math.Abs(imageNonMaxValue) <= Math.Abs(compareTo[1])))
                    //{
                    //    imageNonMaxValue = 0;
                    //}
                    //if (Math.Abs(imageNonMaxValue) < 0.3) {
                    //    imageNonMaxValue = 0;
                    //}
                    ///

                    ColorFloatPixel p = new ColorFloatPixel();
                    if (imageNonMaxValue < 0)
                    {
                        p.r = 0;
                        p.g = Math.Abs(imageNonMaxValue) * 255;
                        p.b = 0;
                        p.a = Math.Abs(imageNonMaxValue) * 255;
                        imageNonMax[i - 1, j - 1] = p;
                    }
                    else if (imageNonMaxValue > 0)
                    {
                        p.r = 0;
                        p.g = 0;
                        p.b = Math.Abs(imageNonMaxValue) * 255;
                        p.a = Math.Abs(imageNonMaxValue) * 255;
                        imageNonMax[i - 1, j - 1] = p;
                    }
                    else
                    {
                        p.r = 0;
                        p.g = 0;
                        p.b = 0;
                        p.a = 0;
                        imageNonMax[i - 1, j - 1] = p;
                    }


                }
            }
            int [] DirectionMap = new int[] {0, 45, 90, 135, 0, 45, 90, 135, 0};
            for (int i = 1; i < image.Width - 1; i++ )
            {
                for (int j=1; j < image.Height - 1; j++)
                {
                    float magnitudeValue = image[i, j];
                    float direction = directionMap[i, j]/(float)RAD2DEG;

                    double[] compareTo= new double[2];
                    float angle_d=0;
                    if (System.Double.IsNaN( direction))
                    { Console.WriteLine("OOOps!"); }
                    else
                    {
                        angle_d = DirectionMap[(int)(Math.Round(direction / Math.PI * 4) + 4)];
                    }


                    if (angle_d == 0) {
                        compareTo[0] = image[i-1, j];
                        compareTo[1] = image[i+1, j];
                    } else if (angle_d == 45) {
                        compareTo[0] = image[i+1, j+1];
                        compareTo[1] = image[i-1, j-1];
                    } else if (angle_d == 90) {
                        compareTo[0] = image[i, j-1];
                        compareTo[1] = image[i, j+1];
                    } else if (angle_d == 135) {
                        compareTo[0] = image[i+1, j-1];
                        compareTo[1] = image[i-1, j+1];
                    }

                    if ((magnitudeValue <= compareTo[0]) || (magnitudeValue <= compareTo[1]))
                    {
                        imageNonMax_gray[i, j] = 0;

                    }
                    else
                    {
                        imageNonMax_gray[i, j] =  magnitudeValue;

                    }

                }
            }
            ///////////////////////////////////////////////////////////////////////////
            //конец нон максимум сапрешна

            max = imageNonMax_gray[0, 0];
            min = imageNonMax_gray[0, 0];
            for (int j = 0; j < imageNonMax_gray.Height; j++)
            {
                for (int i = 0; i < imageNonMax_gray.Width; i++)
                {
                    if (imageNonMax_gray[i, j] > max) { max = imageNonMax_gray[i, j]; }
                    if (imageNonMax_gray[i, j] < min) { min = imageNonMax_gray[i, j]; }
                }
            }

            for (int j = 0; j < imageNonMax_gray.Height; j++)
            {
                for (int i = 0; i < imageNonMax_gray.Width; i++)
                {
                    imageNonMax_gray[i, j] = ((imageNonMax_gray[i, j]) / (max)) * 256;
                    if (imageNonMax_gray[i, j] > 64)
                    {
                        double theta = directionMap[i, j];
                        if (theta < 0)
                        {
                            theta = theta + 180;//Math.PI;
                        }
                        //else
                        //{
                        int angle = (int)Math.Round(theta);

                        //graphGray[angle, 1]++; 
                        if ((angle != 0) && (angle <= 180) && (angle > 0))
                        {
                            angleMaxHistNonMax[angle]++;
                            file_NonMax.Write(angle + ", ");
                        }
                    }

                }
            }

            file_NonMax.Close();
            ImageIO.ImageToFile(imageNonMax, "Input/Output/Ridge.bmp");
            ImageIO.ImageToFile(imageNonMax_gray, "Input/Output/NMS_gray.bmp");

             maxInd = 0;
             max_angle = angleMaxHistNonMax[0];
             for (int i = 0; i < angleMaxHistNonMax.Length; i++)
            {
                if (max_angle < angleMaxHistNonMax[i])
                {
                    max_angle = angleMaxHistNonMax[i];
                    maxInd = i;
                }
            }
             int sumLeft = 0;
             int sumRight = 0;
             Console.WriteLine("MAX_Angle_NonMax=" + (maxInd ));
             for (int i = 0; i < 91; i++)
             {
                 sumLeft += angleMaxHistNonMax[i];
                 //sumRight += angleMaxHistNonMax[i + 90];
             }
             for (int i = 91; i < angleMaxHistNonMax.Length; i++)
             {
                 //sumLeft += angleMaxHistNonMax[i];
                 sumRight += angleMaxHistNonMax[i];
             }
             Console.WriteLine("MAX_Angle_NonMax=" + (maxInd));
             Console.WriteLine("Left=" + sumLeft+" | "+sumRight+"=Right");


        }
        static Complex[,] mothionBlurGray(float L, float angle)
        {
            int radius = (int)Math.Round(L);

            GrayscaleFloatImage filterKernel = new GrayscaleFloatImage(2 * radius + 1, 2 * radius + 1);
            float sum = 2 * radius + 1;
            float real_sum = 0;
            for (int i = -radius; i <= radius; i++)
            {
                filterKernel[i + radius, radius] = 255;
                //real_sum += filterKernel[i + radius, radius];
            }
            //ImageIO.ImageToFile(filterKernel, @"Input\filterTest\filterMotion.bmp");
            //for (int i = -radius; i <= radius; i++)
            //{
            //    filterKernel[i + radius, radius] /= real_sum; //sum;
            //}

            GrayscaleFloatImage filterMotionLine = testOfMotionBlurFilter.rotateImage(filterKernel, true, angle);
            float normSum = 0.0f;
            for (int i = 0; i < filterMotionLine.Width; i++)
            {
                for (int j = 0; j < filterMotionLine.Height; j++)
                {
                    normSum += filterMotionLine[i, j];
                }
            }
            for (int i = 0; i < filterMotionLine.Width; i++)
            {
                for (int j = 0; j < filterMotionLine.Height; j++)
                {
                    filterMotionLine[i, j] /= normSum;
                }
            }

            Complex[,] complexKernell = new Complex[filterMotionLine.Width, filterMotionLine.Height];
            for (int i = 0; i < complexKernell.GetLength(0); i++)
            {
                for (int j = 0; j < complexKernell.GetLength(1); j++)
                {
                    complexKernell[i, j] = (Complex)filterMotionLine[i, j];
                }
            }
            return complexKernell;
            //ImageIO.ImageToFile(convolution2dGray(image, filterMotionLine), @"Input\filterTest\MotionBluredGray.bmp");

        }

        static void Main(string[] args)
        {
            if (args.Length < 2)
            {
                Console.WriteLine("No required arguments");
                Console.ReadLine();
                return;
            }
            string InputFileName = args[0], OutputFileName = args[1];
            if (!File.Exists(InputFileName))
            {
                Console.WriteLine("File does not exists.");
                Console.ReadLine();
                return;
            }
            GrayscaleFloatImage image = ImageIO.FileToGrayscaleFloatImage(InputFileName);
            MotionBlurParametersEstimation motionBlur = new MotionBlurParametersEstimation(image, mothionBlurGray(15, 45));
            image = motionBlur.addMotionBlurEffect(image);
            ImageIO.ImageToFile(image, OutputFileName);
            ////float sig = 3;//sigma=1 circle
            ////Console.WriteLine("File=" + InputFileName);
            ////Console.WriteLine("Sigma=" + sig);
            ////RidgeDetection(image, sig, (int)(3 * sig));
            ////ImageIO.ImageToFile(image, OutputFileName);
            Console.WriteLine("Done");
            Console.ReadLine();

        }

    }
}
