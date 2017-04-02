using Emgu.CV;
using Emgu.CV.Structure;
using FFTTools;
using ImageReadCS;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;

namespace FFTWtest
{
    public enum FilterMode
    {
        FilterSize,
        FilterStep,
        FilterKernel
    };
    /// <summary>
    ///     Filter bitmap with the Fastest Fourier Transform
    /// </summary>
    public class FilterBuilder : BuilderBase, IBuilder
    {
        private readonly Complex[,] _filterKernel; // filter kernel
        private readonly FilterMode _filterMode; // builder mode
        private readonly double _filterPower; // filter power
        private readonly Size _filterSize; // filter size
        private readonly int _filterStep; // filter step
        private readonly KeepOption _keepOption; // energy save options

        /// <summary>
        ///     Builder constructor
        /// </summary>
        /// <param name="filterStep">Bitmap filter step</param>
        /// <param name="filterPower"></param>
        /// <param name="keepOption"></param>
        public FilterBuilder(int filterStep, double filterPower = 1,
            KeepOption keepOption = KeepOption.AverageAndDelta)
        {
            _filterMode = FilterMode.FilterStep;
            _filterStep = filterStep;
            _filterPower = filterPower;
            _keepOption = keepOption;
        }

        /// <summary>
        ///     Builder constructor
        /// </summary>
        /// <param name="filterKernel"></param>
        /// <param name="filterPower"></param>
        /// <param name="keepOption"></param>//here
        public FilterBuilder(Complex[,] filterKernel, double filterPower = 1,
            KeepOption keepOption = KeepOption.AverageAndDelta)
        {
            _filterMode = FilterMode.FilterKernel;
            _filterKernel = filterKernel;
            _filterPower = filterPower;
            _keepOption = keepOption;
        }

        /// <summary>
        ///     Builder constructor
        /// </summary>
        /// <param name="filterSize">Bitmap filter size</param>
        /// <param name="filterPower"></param>
        /// <param name="keepOption"></param>
        public FilterBuilder(Size filterSize, double filterPower = 1,
            KeepOption keepOption = KeepOption.AverageAndDelta)
        {
            _filterMode = FilterMode.FilterSize;
            _filterSize = filterSize;
            _filterPower = filterPower;
            _keepOption = keepOption;
        }

        /// <summary>
        ///     Performs application-defined tasks associated with freeing, releasing, or resetting unmanaged resources.
        /// </summary>
        public void Dispose()
        {
        }

        /// <summary>
        ///     Vizualize builder
        /// </summary>
        /// <param name="source"></param>
        /// <returns></returns>
        public Bitmap ToBitmap(Bitmap source)
        {
            var width = source.Width;
            var height = source.Height;
            var length = width * height;

            var kernel = new Complex[height, width]; // Kernel values
            GetKernelFourier(kernel);
            Fourier(kernel, FourierDirection.Backward);

            var doubles = new double[kernel.Length];
            var index = 0;
            for (var i = 0; i < height; i++)
                for (var j = 0; j < width; j++)
                    doubles[index++] = 1 + kernel[i, j].Magnitude;

            var max = doubles.Max();
            doubles = doubles.Select(x => Math.Round(255.0 * x / max)).ToArray();
            using (var image = new Image<Gray, double>(width, height))
            {
                Buffer.BlockCopy(doubles, 0, image.Data, 0, length * sizeof(double));
                return image.Convert<Bgr, byte>().ToBitmap();
            }
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
        private double t(int u, int v, float a, float b)
        {
            return Math.PI * (u * a + v * b);
        }

        private Complex motionBlurFunctionFT(int u, int v)
        {
            float T = 1.0f;
            float a = 0.1f;
            float b = 0.1f;
            double m = t(u, v, a, b);
            Complex h;
            if (m!= 0)
                h = new Complex((T / m) * Math.Sin(m) * Math.Cos(m), -T / m * Math.Sin(m) * Math.Sin(m));
            else h = 1;
            if (Double.IsNaN(h.Real)||(Double.IsNaN(h.Imaginary)))
            {
                return 1;
            }
            //Complex h = new Complex((T / m) * Math.Sin(m) * Math.Cos(m), (-T / m) * Math.Sin(m) * Math.Sin(m));
            return h;
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
        ///     Filter bitmap with the Fastest Fourier Transform
        /// </summary>
        /// <returns>Blurred bitmap</returns>
        private Array Apply(Array data, FilterAction filterAction)
        {
            var handle = GCHandle.Alloc(data, GCHandleType.Pinned);

            var f = Complex.One;
            var length = data.Length; // Image length = height*width*colors
            var n0 = data.GetLength(0); // Image height
            var n1 = data.GetLength(1); // Image width
            var n2 = data.GetLength(2); // Image colors

            
            var kernelFourier = new Complex[n0, n1]; // Filter values
            GetKernelFourier(kernelFourier);//теперь здесь пф (на самом деле здесь модуль пф фильтра)
            
            // Filter main loop

            var doubles = new double[length];

            Marshal.Copy(handle.AddrOfPinnedObject(), doubles, 0, doubles.Length);//кажется копируем исходное изображение в дабл

            double average;
            double delta;
            AverageAndDelta(out average, out delta, doubles, _keepOption);

            var complex = doubles.Select(x => new Complex(x, 0)).ToArray();

            var complexInMatrix = new Complex[n0, n1];
            //Fourier(n0, n1, n2, complex, FourierDirection.Forward);//пф исходного изображения?
            Fourier(n0, n1, complex, FourierDirection.Forward);
            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    complexInMatrix[i, j] = kernelFourier[i, j];//complex[i * j + j];
                    //complexInMatrix[i, j] = complex[i * j + j];
                }
            }
            //Fourier(complexInMatrix, FourierDirection.Forward);

            var imageF = new Image<Gray, byte>(n0, n1);

            double[] magnitude = new double[complex.Length];
            for (var i = 0; i < n0; i++)
            {
                magnitude[i] = complex[i].Magnitude;
            }

            //for (var i = 0; i < n0; i++)
            //{
            //    magnitude[i] = complex[i].Magnitude;
            //}
            double max = magnitude.Max();
            //double max = magnitude.Cast<double>().Max();
            double c = 255/ Math.Log (1+max,10);
            Complex[,] magicFunctionkernel = new Complex[n0, n1];
            double[,] magicFunctionkernelMagnitude = new double[n0, n1];
            for (int i = -n0 / 2; i < n0 / 2; i++)
            {
                for (int j = -n1 / 2; j < n1 / 2; j++)
                {
                    magicFunctionkernel[i + n0 / 2, j + n1 / 2] = motionBlurFunctionFT(i, j);
                }
            }
            //for (int i = n0 - 1; i >= 0; i--)
            //{
            //    for (int j = n1 / 2; j >= 0; j--)
            //    {
            //        magicFunctionkernel[n0-i-1, n1 - j - 1] = magicFunctionkernel[i, j];
            //    }
            //}
            //for (int i = n0 / 2; i >= 0; i--)
            //{
            //    for (int j = n1 / 2; j >= 0; j--)
            //    {
            //        magicFunctionkernel[n0 - 1 - i, j] = magicFunctionkernel[i, j];
            //    }
            //}
            //for (int i = n0 / 2; i >= 0; i--)
            //{
            //    for (int j = n1 / 2; j >= 0; j--)
            //    {
            //        magicFunctionkernel[n0 - 1 - i, n1 - j - 1] = magicFunctionkernel[i, j];
            //    }
            //}

            for (int i = 0; i < n0; i++)
            {
                for (int j = 0; j < n1; j++)
                {
                    magicFunctionkernelMagnitude[i,j]= magicFunctionkernel[i, j].Magnitude;
                }
            }
            double maxMagic = magicFunctionkernelMagnitude.Cast<double>().Max();
            double cMagic = 255 / Math.Log(1 + maxMagic);
                for (var i = 0; i < n0; i++)
                {
                    for (var j = 0; j < n1; j++)
                    {
                        Complex value = magicFunctionkernel[i, j];//complex[i * n1 + j] * motionBlurFunctionFT(i - n0 / 2, j - n1 / 2);//* magicFunctionkernel[i, j];////magicFunctionkernel[i, j]; //complex[i * n1 + j] * magicFunctionkernel[i, j];//* kernelFourier[i, j];//motionBlurFunctionFT(i, j);
                        complex[i * n1 + j] *= magicFunctionkernel[i, j];//magicFunctionkernel[i, j];//kernelFourier[i, j];//motionBlurFunctionFT(i , j );
                        //imageF[i, j] = new Gray(value.Magnitude / Math.Sqrt(n0 * n1) );
                        //imageF[i, j] =new Gray(value.Magnitude/Math.Sqrt(n0*n1));//работает
                        //imageF[i, j] = new Gray(c * Math.Log(value.Magnitude, 10));//работает
                        imageF[i, j] = new Gray(cMagic * Math.Log(value.Magnitude));//new Gray(value.Magnitude / Math.Sqrt(n0 * n1));//new Gray(cMagic * Math.Log(value.Magnitude, 10));
                        //imageF[i, j] = new Gray(Complex.Pow(value * Complex.Conjugate(value), _filterPower / 2).real );
                    }
                }
            imageF.Save(@"InputImage\FT.bmp");

            //// Apply filter values
            //var index = 0;
            //switch (filterAction)
            //{
            //    case FilterAction.Multiply:
            //        for (var i = 0; i < n0; i++)
            //            for (var j = 0; j < n1; j++)
            //            {
            //                var value = kernelFourier[i, j];
            //                for (var k = 0; k < n2; k++)
            //                    complex[index++] *= f + value;
            //            }
            //        break;
            //    case FilterAction.Divide:
            //        for (var i = 0; i < n0; i++)
            //            for (var j = 0; j < n1; j++)
            //            {
            //                var value = kernelFourier[i, j];
            //                for (var k = 0; k < n2; k++)
            //                    complex[index++] /= f + value;
            //            }
            //        break;
            //    default:
            //        throw new NotImplementedException();
            //}

            Fourier(n0, n1, n2, complex, FourierDirection.Backward);
            doubles = complex.Select(x => x.Magnitude / length).ToArray();
            var imageFBackward = new Image<Gray, byte>(n0, n1);

            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    Complex value = complex[i * n1 + j];
                    //imageF[i, j] =new Gray( value.Magnitude);
                    imageFBackward[i, j] = new Gray(value.Magnitude/(n0*n1));
                    //imageF[i, j] = new Gray(Complex.Pow(value * Complex.Conjugate(value), _filterPower / 2).real );
                }
            }
            imageFBackward.Save(@"InputImage\FT_back.bmp");
            double average2;
            double delta2;
            AverageAndDelta(out average2, out delta2, doubles, _keepOption);

            // a*average2 + b == average
            // a*delta2 == delta
            ////var a = (_keepOption == KeepOption.AverageAndDelta) ? (delta / delta2) : (average / average2);
            ////var b = (_keepOption == KeepOption.AverageAndDelta) ? (average - a * average2) : 0;
            ////Debug.Assert(Math.Abs(a * average2 + b - average) < 0.1);
            ////doubles = doubles.Select(x => Math.Round(a * x + b)).ToArray();

            Marshal.Copy(doubles, 0, handle.AddrOfPinnedObject(), doubles.Length);

            handle.Free();

            return data;
        }
        public Array BlurFromSource(Array data)
        {
            var handle = GCHandle.Alloc(data, GCHandleType.Pinned);

            var length = data.Length;
            var n0 = data.GetLength(0);
            var n1 = data.GetLength(1);
            var n2 = data.GetLength(2);

            var doubles = new double[length];

            Marshal.Copy(handle.AddrOfPinnedObject(), doubles, 0, doubles.Length);

            double average;
            double delta;
            AverageAndDelta(out average, out delta, doubles, _keepOption);

            var complex = doubles.Select(x => new Complex(x, 0)).ToArray();

            Fourier(n0, n1, n2, complex, FourierDirection.Forward);
            var level = complex[0];
            var filterSize = _filterSize;
            switch (_filterMode)
            {
                case FilterMode.FilterSize:
                    break;
                case FilterMode.FilterStep:
                    var filterStep = _filterStep;
                    filterSize = new Size(MulDiv(n1, filterStep, filterStep + 1),
                        MulDiv(n0, filterStep, filterStep + 1));
                    break;
                default:
                    throw new NotImplementedException();
            }

            BlindInner(n0, n1, complex, filterSize, n2);
            complex[0] = level;
            Fourier(n0, n1, n2, complex, FourierDirection.Backward);
            doubles = complex.Select(x => x.Magnitude).ToArray();

            double average2;
            double delta2;
            AverageAndDelta(out average2, out delta2, doubles, _keepOption);

            // a*average2 + b == average
            // a*delta2 == delta
            var a = (_keepOption == KeepOption.AverageAndDelta) ? (delta / delta2) : (average / average2);
            var b = (_keepOption == KeepOption.AverageAndDelta) ? (average - a * average2) : 0;
            Debug.Assert(Math.Abs(a * average2 + b - average) < 0.1);
            doubles = doubles.Select(x => Math.Round(a * x + b)).ToArray();

            Marshal.Copy(doubles, 0, handle.AddrOfPinnedObject(), doubles.Length);

            handle.Free();

            return data;
        }
        private Array ApplyMotionBlur(Array data, FilterAction filterAction)
        {
            var handle = GCHandle.Alloc(data, GCHandleType.Pinned);

            var f = Complex.One;
            var length = data.Length; // Image length = height*width*colors
            var n0 = data.GetLength(0); // Image height
            var n1 = data.GetLength(1); // Image width
            var n2 = data.GetLength(2); // Image colors

            var kernelFourier = new Complex[n0, n1]; //Filter values
            GetKernelFourier(kernelFourier);//теперь здесь пф (на самом деле здесь модуль пф фильтра)

            // Filter main loop
            var doubles = new double[length];
            Marshal.Copy(handle.AddrOfPinnedObject(), doubles, 0, doubles.Length);//кажется копируем исходное изображение в дабл

            var complex = doubles.Select(x => new Complex(x, 0)).ToArray();
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

            var imageF = new Image<Gray, byte>(n0, n1);
            var imageFconv = new Image<Gray, byte>(n0, n1);

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
                    imageF[i, j] =new Gray(value.Magnitude/Math.Sqrt(n0*n1));//работает
                    imageFconv[i, j] = new Gray(valueConv.Magnitude / Math.Sqrt(n0 * n1));
                    //imageF[i, j] = new Gray(c * Math.Log(value.Magnitude, 10));//работает
                    //imageF[i, j] = new Gray(cMagic * Math.Log(value.Magnitude));
                }
            }
            imageF.Save(@"InputImage\FT.bmp");
            imageFconv.Save(@"InputImage\FT_aconv.bmp");

            Fourier(n0, n1, n2, complex, FourierDirection.Backward);
            doubles = complex.Select(x => x.Magnitude / length).ToArray();
            var imageFBackward = new Image<Gray, byte>(n0, n1);

            for (var i = 0; i < n0; i++)
            {
                for (var j = 0; j < n1; j++)
                {
                    Complex value = complex[i * n1 + j] *Math.Pow(-1, i + j);
                    //imageF[i, j] =new Gray( value.Magnitude);
                    imageFBackward[i, j] = new Gray(value.Real / (n0 * n1));
                    //imageF[i, j] = new Gray(Complex.Pow(value * Complex.Conjugate(value), _filterPower / 2).real );
                }
            }
            imageFBackward.Save(@"InputImage\FT_back.bmp");

            Marshal.Copy(doubles, 0, handle.AddrOfPinnedObject(), doubles.Length);

            handle.Free();

            return data;
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
        ///     Blur bitmap with the Fastest Fourier Transform
        /// </summary>
        /// <returns>Blurred bitmap</returns>
        public Image<TColor, TDepth> Blur<TColor, TDepth>(Image<TColor, TDepth> bitmap)
            where TColor : struct, IColor
            where TDepth : new()
        {
            using (var image = bitmap.Convert<TColor, double>())
            {
                image.Data = Apply(image.Data, FilterAction.Multiply) as double[, ,];
                return image.Convert<TColor, TDepth>();
            }
        }


        public Image<TColor, TDepth> BlurMotionBlur<TColor, TDepth>(Image<TColor, TDepth> bitmap)
            where TColor : struct, IColor
            where TDepth : new()
        {
            using (var image = bitmap.Convert<TColor, double>())
            {
                image.Data = ApplyMotionBlur(image.Data, FilterAction.Multiply) as double[, ,];
                return image.Convert<TColor, TDepth>();
            }
        }

        /// <summary>
        ///     Blur bitmap with the Fastest Fourier Transform
        /// </summary>
        /// <returns>Blurred bitmap</returns>
        public Bitmap Blur(Bitmap bitmap)
        {
            using (var image = new Image<Bgr, double>(bitmap))
            {
                image.Data = Apply(image.Data, FilterAction.Multiply) as double[, ,];
                return image.ToBitmap();
            }
        }

        /// <summary>
        ///     Sharp bitmap with the Fastest Fourier Transform
        /// </summary>
        /// <returns>Sharped bitmap</returns>
        public Image<TColor, TDepth> Sharp<TColor, TDepth>(Image<TColor, TDepth> bitmap)
            where TColor : struct, IColor
            where TDepth : new()
        {
            using (var image = bitmap.Convert<TColor, double>())
            {
                image.Data = Apply(image.Data, FilterAction.Divide) as double[, ,];
                return image.Convert<TColor, TDepth>();
            }
        }

        /// <summary>
        ///     Sharp bitmap with the Fastest Fourier Transform
        /// </summary>
        /// <returns>Sharped bitmap</returns>
        public Bitmap Sharp(Bitmap bitmap)
        {
            using (var image = new Image<Bgr, double>(bitmap))
            {
                image.Data = Apply(image.Data, FilterAction.Divide) as double[, ,];
                return image.ToBitmap();
            }
        }
    }

}
