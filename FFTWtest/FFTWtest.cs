using System;
using System.Runtime.InteropServices;
using FFTWSharp;
using System.Linq;
using System.Numerics;
using FFTTools;
using Emgu.CV;
using Emgu.CV.Structure;
using System.Security.Cryptography;
using System.Text;
using FFTWtest;
using ImageReadCS;

namespace FFTWSharp_test
{
    public class FFTWtest
    {
        const int sampleSize = 16384;
        const int repeatPlan = 10000;
        const bool forgetWisdom = false;

        //pointers to unmanaged arrays
        IntPtr pin, pout;

        //managed arrays
        float[] fin, fout;
        double[] din, dout;

        //handles to managed arrays, keeps them pinned in memory
        GCHandle hin, hout, hdin, hdout;

        //pointers to the FFTW plan objects
        IntPtr fplan1, fplan2, fplan3, fplan4, fplan5;

        // and an example of the managed interface
        fftwf_plan mplan1, mplan2, mplan3, mplan4;
        fftwf_complexarray mfin, mfout;
        fftw_plan mplan5;
        fftw_complexarray mdin, mdout;

        // length of arrays
        int fftLength = 0;

        // Initializes FFTW and all arrays
        // n: Logical size of the transform
        public FFTWtest(int n)
        {
            Console.WriteLine("Start testing with n = {n.ToString(#,0)} complex numbers. All plans will be executed {repeatPlan.ToString(#,0)} times on a single thread.");
            Console.WriteLine("Please wait, creating plans...");
            fftLength = n;

            // create two unmanaged arrays, properly aligned
            pin = fftwf.malloc(n * 8);
            pout = fftwf.malloc(n * 8);

            // create two managed arrays, possibly misalinged
            // n*2 because we are dealing with complex numbers
            fin = new float[n * 2];
            fout = new float[n * 2];
            // and two more for double FFTW
            din = new double[n * 2];
            dout = new double[n * 2];

            // get handles and pin arrays so the GC doesn't move them
            hin = GCHandle.Alloc(fin, GCHandleType.Pinned);
            hout = GCHandle.Alloc(fout, GCHandleType.Pinned);
            hdin = GCHandle.Alloc(din, GCHandleType.Pinned);
            hdout = GCHandle.Alloc(dout, GCHandleType.Pinned);

            // create a few test transforms
            fplan1 = fftwf.dft_1d(n, pin, pout, fftw_direction.Forward, fftw_flags.Estimate);
            fplan2 = fftwf.dft_1d(n, hin.AddrOfPinnedObject(), hout.AddrOfPinnedObject(), fftw_direction.Forward, fftw_flags.Estimate);
            fplan3 = fftwf.dft_1d(n, hout.AddrOfPinnedObject(), pin, fftw_direction.Backward, fftw_flags.Measure);
            // end with transforming back to original array
            fplan4 = fftwf.dft_1d(n, hout.AddrOfPinnedObject(), hin.AddrOfPinnedObject(), fftw_direction.Backward, fftw_flags.Estimate);
            // and check a quick one with doubles, just to be sure
            fplan5 = fftw.dft_1d(n, hdin.AddrOfPinnedObject(), hdout.AddrOfPinnedObject(), fftw_direction.Backward, fftw_flags.Measure);

            // create a managed plan as well
            mfin = new fftwf_complexarray(fin);
            mfout = new fftwf_complexarray(fout);
            mdin = new fftw_complexarray(din);
            mdout = new fftw_complexarray(dout);

            mplan1 = fftwf_plan.dft_1d(n, mfin, mfout, fftw_direction.Forward, fftw_flags.Estimate);
            mplan2 = fftwf_plan.dft_1d(n, mfin, mfout, fftw_direction.Forward, fftw_flags.Measure);
            mplan3 = fftwf_plan.dft_1d(n, mfin, mfout, fftw_direction.Forward, fftw_flags.Patient);

            mplan4 = fftwf_plan.dft_1d(n, mfout, mfin, fftw_direction.Backward, fftw_flags.Measure);
            
            mplan5 = fftw_plan.dft_1d(n, mdin, mdout, fftw_direction.Forward, fftw_flags.Measure);


            // fill our arrays with an arbitrary complex sawtooth-like signal
            for (int i = 0; i < n * 2; i++) fin[i] = i % 50;
            for (int i = 0; i < n * 2; i++) fout[i] = i % 50;
            for (int i = 0; i < n * 2; i++) din[i] = i % 50;

            // copy managed arrays to unmanaged arrays
            Marshal.Copy(fin, 0, pin, n * 2);
            Marshal.Copy(fout, 0, pout, n * 2);




            Console.WriteLine();
        }

        public void TestAll()
        {
            TestPlan(fplan1, "single (malloc) | pin,  pout,  Forward,  Estimate");
            TestPlan(fplan2, "single (pinned) | hin,  hout,  Forward,  Estimate");
            TestPlan(fplan3, "single (pinned) | hout, pin,   Backward, Measure ");

            // set fin to 0, and try to refill it from a backwards fft from fout (aka hin/hout)
            for (int i = 0; i < fftLength * 2; i++) fin[i] = 0;

            TestPlan(fplan4, "single (pinned) | hout, hin,   Backward, Estimate");

            // check and see how we did, don't say anyt
            for (int i = 0; i < fftLength * 2; i++)
            {
                // check against original values
                // note that we need to scale down by length, due to FFTW scaling by N
                if (System.Math.Abs(fin[i]/fftLength - (i % 50)) > 1e-3)
                {
                    System.Console.WriteLine("FFTW consistency error!");
                    return;
                }
            }

            TestPlan(fplan5, "double (pinned) | hdin, hdout, Backward, Measure ");


            Console.WriteLine();
            TestPlan(mplan1, "#1 single (managed) | mfin, mfout, Forward,  Estimate");
            TestPlan(mplan2, "#2 single (managed) | mfin, mfout, Forward,  Measure ");
            TestPlan(mplan3, "#3 single (managed) | mfin, mfout, Forward,  Patient ");
            Console.WriteLine();

            // fill our input array with an arbitrary complex sawtooth-like signal
            for (int i = 0; i < fftLength * 2; i++) fin[i] = i % 50;
            for (int i = 0; i < fftLength * 2; i++) fout[i] = 0;

            mfin.SetData(fin);
            mfout.SetData(fout);
            TestPlan(mplan2, "#2 single (managed) | mfin, mfout, Forward,  Measure ");

            fout = mfout.GetData_Float();  // let's see what's in mfout
                                           // at this point mfout contains the FFT'd mfin

            TestPlan(mplan4, "#4 single (managed) | mfout, mfin, Backward, Measure ");
            // at this point we have transfarred backwards the mfout into mfin, so mfin should be very close to the original complex sawtooth-like signal

            fin = mfin.GetData_Float();
            for (int i = 0; i < fftLength * 2; i++) fin[i] /= fftLength;
            // at this point fin should be very close to the original (sawtooth-like) signal

            // check and see how we did, don't say anyt
            for (int i = 0; i < fftLength * 2; i++)
            {
                // check against original values
                if (System.Math.Abs(fin[i] - (i % 50)) > 1e-3)
                {
                    System.Console.WriteLine("FFTW consistency error!");
                    return;
                }
            }

            Console.WriteLine();

            TestPlan(mplan2, "#2 single (managed) | mfin, mfout, Forward,  Measure ");
            TestPlan(mplan5, "#5 double (managed) | mdin, mdout, Forward,  Measure ");
            Console.WriteLine();
        }

        // Tests a single plan, displaying results
        // plan: Pointer to plan to test
        public void TestPlan(object plan, string planName)
        {
            // a: adds, b: muls, c: fmas
            double a = 0, b = 0, c = 0;

            int start = System.Environment.TickCount;

            if (plan is IntPtr)
            {
                IntPtr umplan = (IntPtr)plan;

                for (int i = 0; i < repeatPlan; i++)
                {
                        fftwf.execute(umplan);
                }

                fftwf.flops(umplan, ref a, ref b, ref c);
            }

            if (plan is fftw_plan)
            {
                fftw_plan mplan = (fftw_plan)plan;

                for (int i = 0; i < repeatPlan; i++)
                {
                    mplan.Execute();
                }

                fftw.flops(mplan.Handle, ref a, ref b, ref c);
            }

            if (plan is fftwf_plan)
            {
                fftwf_plan mplan = (fftwf_plan)plan;

                for (int i = 0; i < repeatPlan; i++)
                {
                    mplan.Execute();
                }

                fftwf.flops(mplan.Handle, ref a, ref b, ref c);
            }

            double mflops = (((a + b + 2 * c)) * repeatPlan) / (1024 * 1024);
            long ticks = (System.Environment.TickCount - start);

            //Console.WriteLine($"Plan '{planName}': {ticks.ToString("#,0")} us | mflops: {FormatNumber(mflops)} | mflops/s: {(1000*mflops/ticks).ToString("#,0.0")}");
            Console.WriteLine("Plan '{0}': {1,8:N0} us | mflops: {2,8:N0} | mflops/s: {3,8:N0}", planName, ticks, mflops, (1000 * mflops / ticks));
        }

        // Releases all memory used by FFTW/C#
        ~FFTWtest()
        {
            fftwf.export_wisdom_to_filename("wisdom.wsd");

            // it is essential that you call these after finishing
            // that's why they're in the destructor. See Also: RAII
            fftwf.free(pin);
            fftwf.free(pout);
            fftwf.destroy_plan(fplan1);
            fftwf.destroy_plan(fplan2);
            fftwf.destroy_plan(fplan3);
            fftwf.destroy_plan(fplan4);
            fftwf.destroy_plan(fplan5);
            hin.Free();
            hout.Free();
        }
        static GrayscaleFloatImage rotateImage(GrayscaleFloatImage image, bool isCW, float angleDeg)
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

        static Complex[,] mothionBlurGray( float L, float angle)
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

            GrayscaleFloatImage filterMotionLine = rotateImage(filterKernel, true, angle);
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

        private static readonly Random Rnd = new Random((int)DateTime.Now.Ticks);
        //static void Main(string[] args)
        //{
        //    if (forgetWisdom)
        //    {
        //        fftwf.fftwf_forget_wisdom();
        //    }
        //    else
        //    {
        //        Console.WriteLine("Importing wisdom (wisdom speeds up the plan creation process, if that plan was previously created at least once)");
        //        fftwf.import_wisdom_from_filename("wisdom.wsd");
        //    }

        //    // initialize our FFTW test class
        //    FFTWtest test = new FFTWtest(sampleSize);

        //    // run the tests, print debug output
        //    //test.TestAll();


        //    //const int length = 1000;
        //    //var f = Enumerable.Range(0, length).Select(i => Rnd.NextDouble()).ToArray();
        //    //var complex = f.Select(x => new Complex(x, 0)).ToArray();
        //    //BuilderBase.Fourier(complex, FourierDirection.Forward);
        //    //complex = complex.Select(x => Complex.Conjugate(x) / length).ToArray();
        //    //BuilderBase.Fourier(complex, FourierDirection.Backward);
        //    //var f1 = complex.Select(x => x.Magnitude).ToArray();

        //    //float[,] blurKernel = new float[,] { 
        //    //{1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f}, 
        //    //{1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f},
        //    //{1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f}, 
        //    //{1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f},
        //    //{1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f, 1.0f / 25.0f}};
        //    //{1.0f / 16.0f, 1.0f / 16.0f, 1.0f / 16.0f}, 
        //    //{1.0f / 16.0f, 8.0f / 16.0f, 1.0f / 16.0f}, 
        //    //{1.0f / 16.0f, 1.0f / 16.0f, 1.0f / 16.0f}};
        //    int n0 = 11;
        //    int n1 = 11;
        //    float[,] blurKernel = new float[n0, n1];
        //    for (int i = -n0 / 2; i < n0 / 2; i++)
        //    {
        //        for (int j = -n1 / 2; j < n1 / 2; j++)
        //        {
        //            if ((i != 0) && (j != 0))
        //                blurKernel[i + n0 / 2, j + n1 / 2] = (float)((1 / (Math.PI * 0.1 * (i + j))) * Math.Sin(Math.PI * 0.1 * (i + j)));
        //            else blurKernel[i + n0 / 2, j + n1 / 2] = 0;
        //        }
        //    }
        //    Complex[,] complexKernel = new Complex[blurKernel.GetLength(0), blurKernel.GetLength(1)];
        //    for (int i = 0; i < blurKernel.GetLength(0); i++)
        //    {
        //        for (int j = 0; j < blurKernel.GetLength(1); j++)
        //        {
        //            complexKernel[i, j] = (Complex)blurKernel[i, j];
        //        }
        //    }

        //    const int width = 100;
        //    const int height = 100;
        //    const int filterStep = 1;
        //    RandomNumberGenerator Rng = RandomNumberGenerator.Create();

        //    //var image = new Image<Gray, byte>(height, width);
        //    var image = new Image<Gray, byte>(@"InputImage\lena.bmp");

        //    //var bytes = new byte[image.Data.Length];

        //    //Rng.GetBytes(bytes);

        //    var handle = GCHandle.Alloc(image.Data, GCHandleType.Pinned);
        //    //Marshal.Copy(bytes, 0, handle.AddrOfPinnedObject(), bytes.Length);
        //    handle.Free();

        //    //var average = bytes.Average(x => (double)x);
        //    //var delta = Math.Sqrt(bytes.Average(x => (double)x * x) - average * average);
        //    //var minValue = bytes.Min(x => (double)x);
        //    //var maxValue = bytes.Max(x => (double)x);
        //    //var sb = new StringBuilder();
        //    //sb.AppendLine(string.Format("Length {0}", bytes.Length));
        //    //sb.AppendLine(string.Format("Average {0}", average));
        //    //sb.AppendLine(string.Format("Delta {0}", delta));
        //    //sb.AppendLine(string.Format("minValue {0}", minValue));
        //    //sb.AppendLine(string.Format("maxValue {0}", maxValue));
        //    //Console.WriteLine(sb.ToString());

        //    //using (var filterBuilder = new FilterBuilder(filterStep))
        //    //    image = filterBuilder.Blur(image);
        //    using (var filterBuilder = new FilterBuilder(mothionBlurGray(15, 45)))
        //        image = filterBuilder.BlurMotionBlur(image);
        //    //using (var filterBuilder = new FilterBuilder(complexKernel))
        //    //    image = filterBuilder.BlurMotionBlur(image);
        //    //using (var filterBuilder = new FilterBuilder(complexKernel))
        //    //    image = filterBuilder.BlurFromSource(image.Select(x => new Complex(x, 0)).ToArray());

        //    handle = GCHandle.Alloc(image.Data, GCHandleType.Pinned);
        //    //Marshal.Copy(handle.AddrOfPinnedObject(), bytes, 0, bytes.Length);
        //    handle.Free();

        //    //var average1 = bytes.Average(x => (double)x);
        //    //var delta1 = Math.Sqrt(bytes.Average(x => (double)x * x) - average1 * average1);
        //    //var minValue1 = bytes.Min(x => (double)x);
        //    //var maxValue1 = bytes.Max(x => (double)x);
        //    //var sb1 = new StringBuilder();
        //    //sb1.AppendLine(string.Format("Length {0}", bytes.Length));
        //    //sb1.AppendLine(string.Format("Average {0}", average1));
        //    //sb1.AppendLine(string.Format("Delta {0}", delta1));
        //    //sb1.AppendLine(string.Format("minValue {0}", minValue1));
        //    //sb1.AppendLine(string.Format("maxValue {0}", maxValue1));
        //    //Console.WriteLine(sb1.ToString());

        //    image.Save(@"InputImage\OutImg.bmp");

        //    //Assert.IsTrue(delta1 < delta);
        //    // pause for user input, then quit
        //    System.Console.WriteLine("\nDone. Press any key to exit.");
        //    String str = System.Console.ReadLine();
        //}
    }
}