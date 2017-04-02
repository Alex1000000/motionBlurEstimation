using Emgu.CV;
using FFTTools;
using ImageReadCS;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using FFTWSharp;
using System.Threading;
using System.Drawing;

namespace FFTWtest
{
    class MotionBlurParametersEstimation : fftw
    {
        private GrayscaleFloatImage imageToAnalize;
        private static readonly Mutex FftwLock = new Mutex();
        private readonly FilterMode _filterMode; // builder mode
        private readonly Complex[,] _filterKernel; // filter kernel
        private readonly double _filterPower; // filter power
        private readonly Size _filterSize; // filter size
        private readonly int _filterStep; // filter step
        private readonly KeepOption _keepOption; // energy save options

        public MotionBlurParametersEstimation(GrayscaleFloatImage imageToAnalize, Complex[,] filterKernel)
        {
            this.imageToAnalize = imageToAnalize;
            this._filterMode = FilterMode.FilterKernel;
            _filterMode = FilterMode.FilterKernel;
            _filterKernel = filterKernel;//Complex[,] 
            _filterPower = 1;
            _keepOption = KeepOption.AverageAndDelta;
        }

        public GrayscaleFloatImage ImageToAnalize
        {
            get { return imageToAnalize; }
            set { imageToAnalize = value; }
        }

        public GrayscaleFloatImage addMotionBlurEffect(GrayscaleFloatImage imageToAddEffect)
        {
            GrayscaleFloatImage imageWithMotionBlur = new GrayscaleFloatImage(imageToAddEffect.Width, imageToAddEffect.Height);
            //var image = new Image<Gray, byte>(imageToAddEffect.Height, width);
             
            //var handle = GCHandle.Alloc(imageToAddEffect, GCHandleType.Pinned);
            var f = Complex.One;
            var length = imageToAddEffect.Height*imageToAddEffect.Width*1; // Image length = height*width*colors
            var n0 = imageToAddEffect.Height; // Image height
            var n1 = imageToAddEffect.Width; // Image width
            var n2 = 1; // Image colors

            var kernelFourier = new Complex[n0, n1]; //Filter values
            GetKernelFourier(kernelFourier);//теперь здесь пф (на самом деле здесь модуль пф фильтра)

            // Filter main loop
            //var doubles = new double[length];
            //Marshal.Copy(handle.AddrOfPinnedObject(), doubles, 0, doubles.Length);//кажется копируем исходное изображение в дабл
            Complex[] complex = new Complex[n0 * n1];
            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    complex[i * n1 + j] = imageToAddEffect[i ,j];
                }
            }
            //var complex = doubles.Select(x => new Complex(x, 0)).ToArray();
            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    complex[i * n1 + j] *= Math.Pow(-1, i + j);
                }
            }

            var complexInMatrix = new Complex[n0, n1];
            //Fourier(n0, n1, n2, complex, FourierDirection.Forward);//пф исходного изображения
            Fourier(n0, n1, complex, FourierDirection.Forward);
            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    complexInMatrix[i, j] = kernelFourier[i, j];//complex[i * j + j];
                }
            }
            //Fourier(complexInMatrix, FourierDirection.Forward);

            //var imageF = new Image<Gray, byte>(n0, n1);
            //var imageFconv = new Image<Gray, byte>(n0, n1);

            double[] magnitude = new double[complex.Length];
            for (var i = 0; i < n0; i++)
            {
                magnitude[i] = complex[i].Magnitude;
            }

            double max = magnitude.Max();
            //double max = magnitude.Cast<double>().Max();
            double c = 255 / Math.Log(1 + max, 10);
            Complex[,] magicFunctionkernel = new Complex[n0, n1];
            double[,] magicFunctionkernelMagnitude = new double[n0, n1];
            for (int i = -n0 / 2; i < n0 / 2; i++)
            {
                for (int j = -n1 / 2; j < n1 / 2; j++)
                {
                    magicFunctionkernel[i + n0 / 2, j + n1 / 2] = motionBlurFunctionFT(i, j);
                }
            }
            //for (int i = 0; i < n0 ; i++)
            //{
            //    for (int j = 0; j < n1 ; j++)
            //    {
            //        magicFunctionkernel[i , j ] = motionBlurFunctionFT(i, j);
            //    }
            //}
            

            for (int i = 0; i < n0; i++)
            {
                for (int j = 0; j < n1; j++)
                {
                    magicFunctionkernelMagnitude[i, j] = magicFunctionkernel[i, j].Magnitude;
                }
            }
            double maxMagic = magicFunctionkernelMagnitude.Cast<double>().Max();
            double cMagic = 255 / Math.Log(1 + maxMagic);
            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    Complex valueConv = complex[i * n1 + j] * magicFunctionkernel[i, j];
                    Complex value = complex[i * n1 + j];////magicFunctionkernel[i, j]; //complex[i * n1 + j] * magicFunctionkernel[i, j];//* kernelFourier[i, j];//motionBlurFunctionFT(i, j);
                    complex[i * n1 + j] *= magicFunctionkernel[i, j];//kernelFourier[i, j];//motionBlurFunctionFT(i , j );
                    //imageF[i, j] = new Gray(value.Magnitude / Math.Sqrt(n0 * n1) );
                    //imageF[i, j] =new Gray(value.Magnitude/Math.Sqrt(n0*n1));//работает
                    //imageFconv[i, j] = new Gray(valueConv.Magnitude / Math.Sqrt(n0 * n1));
                    //imageF[i, j] = new Gray(c * Math.Log(value.Magnitude, 10));//работает
                    //imageF[i, j] = new Gray(cMagic * Math.Log(value.Magnitude));
                }
            }
            //imageF.Save(@"InputImage\FT.bmp");
            //imageFconv.Save(@"InputImage\FT_aconv.bmp");

            Fourier(n0, n1, n2, complex, FourierDirection.Backward);
            //doubles = complex.Select(x => x.Magnitude / length).ToArray();
            //var imageFBackward = new Image<Gray, byte>(n0, n1);

            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    Complex value = complex[i * n1 + j] *Math.Pow(-1, i + j);
                    //imageF[i, j] =new Gray( value.Magnitude);
                    //imageFBackward[i, j] = new Gray(value.Real / (n0 * n1));
                    //imageF[i, j] = new Gray(Complex.Pow(value * Complex.Conjugate(value), _filterPower / 2).real );
                    imageWithMotionBlur[i,j]= (float)(value.Real / (n0 * n1));
                }
            }
            //imageFBackward.Save(@"InputImage\FT_back.bmp");
            return imageWithMotionBlur;

        }
        private Complex motionBlurFunctionFT(int u, int v)
        {
            float T = 1.0f;
            float a = 0.01f;
            float b = 0.1f;
            double m = t(u, v, a, b);
            Complex h;
            if (m != 0)
                h = new Complex((T / m) * Math.Sin(m) * Math.Cos(m), -T / m * Math.Sin(m) * Math.Sin(m));
            else h = 1;
            if (Double.IsNaN(h.Real) || (Double.IsNaN(h.Imaginary)))
            {
                return 1;
            }
            //Complex h = new Complex((T / m) * Math.Sin(m) * Math.Cos(m), (-T / m) * Math.Sin(m) * Math.Sin(m));
            return h;
        }
        private double t(int u, int v, float a, float b)
        {
            return Math.PI * (u * a + v * b);
        }
        private void GetKernelFourier(Complex[,] kernelFourier)
        {
            var length = kernelFourier.Length; // Kernel length = height*width
            var n0 = kernelFourier.GetLength(0); // Kernel height
            var n1 = kernelFourier.GetLength(1); // Kernel width
            switch (_filterMode)
            {
                case FilterMode.FilterKernel:
                    Copy(_filterKernel, kernelFourier);
                    break;
                case FilterMode.FilterSize:
                    Fill(kernelFourier, _filterSize, Complex.One);
                    break;
                case FilterMode.FilterStep:
                    var filterStep = _filterStep;
                    var filterSize = new Size(MulDiv(n1, 1, filterStep + 1),
                        MulDiv(n0, 1, filterStep + 1));
                    Fill(kernelFourier, filterSize, Complex.One);//заполнили все 1ми
                    break;
                default:
                    throw new NotImplementedException();
            }
            Fourier(kernelFourier, FourierDirection.Forward);//теперь это пф фильтра (т.е. из 1)
            ////kernelFourier-fft
            //for (var i = 0; i < n0; i++)
            //    for (var j = 0; j < n1; j++)
            //    {
            //        var value = kernelFourier[i, j];
            //        value = Complex.Pow(value * Complex.Conjugate(value), _filterPower / 2);
            //        kernelFourier[i, j] = value;//а тут похоже норм модуль 
            //    }
        }

        /// <summary>
        ///     Fourier transform
        /// </summary>
        /// <param name="array">Array</param>
        /// <param name="direction">Fourier direction</param>
        public static void Fourier(Array array, FourierDirection direction)
        {
            var handle = GCHandle.Alloc(array, GCHandleType.Pinned);
            FftwLock.WaitOne();
            var plan = dft(array.Rank, Enumerable.Range(0, array.Rank).Select(array.GetLength).ToArray(),
                handle.AddrOfPinnedObject(), handle.AddrOfPinnedObject(),
                (fftw_direction)direction,
                fftw_flags.Estimate);
            execute(plan);
            destroy_plan(plan);
            FftwLock.ReleaseMutex();
            handle.Free();
        }
        /// <summary>
        ///     Fourier transform
        /// </summary>
        /// <param name="n0">Array size</param>
        /// <param name="n1">Array size</param>
        /// <param name="array">Array</param>
        /// <param name="direction">Fourier direction</param>
        public static void Fourier(int n0, int n1, Array array, FourierDirection direction)
        {
            var handle = GCHandle.Alloc(array, GCHandleType.Pinned);
            FftwLock.WaitOne();
            var plan = dft_2d(n0, n1,
                handle.AddrOfPinnedObject(), handle.AddrOfPinnedObject(),
                (fftw_direction)direction,
                fftw_flags.Estimate);
            execute(plan);
            destroy_plan(plan);
            FftwLock.ReleaseMutex();
            handle.Free();
        }
        /// <summary>
        ///     Copy 2D array to 2D array (sizes can be different)
        /// </summary>
        /// <param name="input">Input array</param>
        /// <param name="output">Output array</param>
        private static void Copy(Complex[,] input, Complex[,] output)
        {
            var n0 = output.GetLength(0);
            var n1 = output.GetLength(1);
            var m0 = Math.Min(n0, input.GetLength(0));
            var m1 = Math.Min(n1, input.GetLength(1));

            for (var i = 0; i < m0; i++)
                for (var j = 0; j < m1; j++)
                    output[i, j] = input[i, j];
        }


        /// <summary>
        /// </summary>
        /// <param name="number"></param>
        /// <param name="numerator"></param>
        /// <param name="denominator"></param>
        /// <returns></returns>
        public static int MulDiv(int number, int numerator, int denominator)
        {
            return (int)(((long)number * numerator) / denominator);
        }

        /// <summary>
        ///     Fill region of array by value
        /// </summary>
        /// <param name="filter">Output array</param>
        /// <param name="size"></param>
        /// <param name="value">Value to replace copied data</param>
        private static void Fill(Complex[,] filter, Size size, Complex value)
        {
            var n0 = filter.GetLength(0);
            var n1 = filter.GetLength(1);
            var m0 = Math.Min(n0 - 1, size.Height);
            var m1 = Math.Min(n1 - 1, size.Width);

            Array.Clear(filter, 0, filter.Length);

            for (var i = 0; i <= m0; i++)
                for (var j = 0; j <= m1; j++)
                    filter[i, j] = value;
        }
        /// <summary>
        ///     Fourier transform
        /// </summary>
        /// <param name="n0">Array size</param>
        /// <param name="n1">Array size</param>
        /// <param name="n2">Array size</param>
        /// <param name="array">Array</param>
        /// <param name="direction">Fourier direction</param>
        public static void Fourier(int n0, int n1, int n2, Array array, FourierDirection direction)
        {
            var handle = GCHandle.Alloc(array, GCHandleType.Pinned);
            FftwLock.WaitOne();
            var plan = dft_3d(n0, n1, n2,
                handle.AddrOfPinnedObject(), handle.AddrOfPinnedObject(),
                (fftw_direction)direction,
                fftw_flags.Estimate);
            execute(plan);
            destroy_plan(plan);
            FftwLock.ReleaseMutex();
            handle.Free();
        }

    }
}
