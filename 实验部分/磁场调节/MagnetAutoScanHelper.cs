﻿using System;
using System.IO;
using MathNet.Numerics;
using System.Security.Cryptography;
using MathNet.Numerics.RootFinding;
using System.Windows.Interop;
using System.Windows.Documents;
using System.Collections.Generic;
using ODMR_Lab.Python.LbviewHandler;
using static System.Net.WebRequestMethods;
using System.Linq;
using PythonHandler;
using ODMR_Lab.磁场调节;

namespace ODMR_Lab.实验部分.磁场调节
{
    /// <summary>
    /// 圆柱形磁铁类
    /// </summary>
    public class MagnetAutoScanHelper
    {

        // 旋磁比(Mhz/Gauss)
        private static double gammaE = 2.803;

        /// <summary>
        /// 根据谱峰位置计算横向和纵向的磁场大小
        /// </summary>
        /// <param name="freq1"></param>
        /// <param name="freq2"></param>
        /// <param name="D"></param>
        /// <param name="Bp"></param>
        /// <param name="Bv"></param>
        public static void CalculateB(double freq1, double freq2, out double Bp, out double Bv, double D = 2870)
        {
            double detplus = (freq2 + freq1) / 2;
            double detminus = Math.Abs(freq2 - freq1);
            double det1 = 2 * D * (detplus - D) / 3;
            if (det1 < 0)
                det1 = 0;
            Bv = Math.Sqrt(det1) / gammaE;
            double det2 = detminus / Math.Pow(2 * gammaE, 2) - Math.Pow(Bv * Bv * gammaE / (2 * D), 2);
            if (det2 < 0)
                det2 = 0;
            Bp = Math.Sqrt(det2);
        }

        /// <summary>
        /// 寻找Z方向的根
        /// </summary>
        /// <param name="M"></param>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <param name="ratiovalue"></param>
        /// <returns></returns>
        public static double FindRoot(PillarMagnet M, double z1, double z2, double ratiovalue)
        {
            dynamic result = Python_NetInterpretor.ExcuteFunction(Path.Combine(Environment.CurrentDirectory, "Python", "LabviewHandler", "Magnet.py"), "FindRoot", TimeSpan.FromSeconds(5), M.Radius, M.Length, z1, z2, ratiovalue);
            return result;
        }

        /// <summary>
        /// 获取初始偏移量
        /// </summary>
        /// <param name="angle1">初始角度</param>
        /// <param name="x1">初始位置X</param>
        /// <param name="y1">初始位置Y</param>
        /// <param name="angle2">转180度的角度</param>
        /// <param name="x2">之后的位置X</param>
        /// <param name="y2">之后的位置Y</param>
        /// <param name="reverseX">是否反转X</param>
        /// <param name="reverseY">是否反转Y</param>
        /// <param name="reverseAngle">是否反转角度</param>
        /// <param name="initAngle">角度起始点</param>
        /// <param name="offsetx">X偏移量</param>
        /// <param name="offsety">Y偏移量</param>
        public static void GetZeroOffset(double angle1, double x1, double y1, double angle2, double x2, double y2, bool reverseX, bool reverseY, bool reverseAngle, double initAngle, out double offsetx, out double offsety)
        {
            double dx = x2 - x1;
            double dy = y2 - y1;

            double reverseFacX = reverseX ? -1 : 1;
            double reverseFacY = reverseY ? -1 : 1;
            double reverseFacAngle = reverseAngle ? -1 : 1;

            double x = 0.5 * reverseFacX * dx;
            double y = 0.5 * reverseFacY * dy;

            double r = Math.Sqrt(x * x + y * y);

            double theta = Math.Atan2(y, x) * 180 / Math.PI;

            double angle = theta + reverseFacAngle * (initAngle - angle1);


            offsetx = r * Math.Cos(angle * Math.PI / 180);
            offsety = r * Math.Sin(angle * Math.PI / 180);
        }

        public static List<double> GetTargetOffset(MagnetScanParams param, double targetAngle)
        {

            // 计算随体坐标系下的极坐标
            double r = Math.Sqrt(param.OffsetX.Value * param.OffsetX.Value + param.OffsetY.Value * param.OffsetY.Value);
            double the = Math.Atan2(param.OffsetY.Value, param.OffsetX.Value) * 180 / Math.PI;

            double the1 = the + param.ReverseANum.Value * (targetAngle - param.StartAngle.Value);

            double dx = param.ReverseXNum.Value * r * (Math.Cos(the / 180 * Math.PI) - Math.Cos(the1 / 180 * Math.PI));
            double dy = param.ReverseYNum.Value * r * (Math.Sin(the / 180 * Math.PI) - Math.Sin(the1 / 180 * Math.PI));

            return new List<double>() { dx, dy };
        }


        /// <summary>
        /// 指定方向和给定的Z方向距离，找出对应方向上的磁场，相对于角度基准点（磁铁朝向为X轴）要转动的角度，以及XY平面上相对于原点的坐标
        /// </summary>
        /// <param name="r0"></param>
        /// <param name="l0"></param>
        /// <param name="targetthe"></param>
        /// <param name="targetphi"></param>
        /// <param name="distance"></param>
        public static dynamic FindDire(double r0, double l0, double targetthe, double targetphi, double distance)
        {
            dynamic result = Python_NetInterpretor.ExcuteFunction(Path.Combine(Environment.CurrentDirectory, "Python", "LabviewHandler", "Magnet.py"), "FindDire", TimeSpan.FromSeconds(5), r0, l0, targetthe, targetphi, distance);
            return result;
        }

        private static void TotalCW(out List<double> peaks, out List<double> freqs1, out List<double> contracts1, out List<double> freqs2, out List<double> contracts2)
        {
            //粗扫
            LabviewConverter.CW(out List<double> freqs, out List<double> contracts, out List<double> ps, out List<double> fitcontracts, out Exception exc, 2600, 3200, 3, 500, 2);

            if (ps.Count == 1)
            {
                peaks = new List<double>();
                //细扫
                LabviewConverter.CW(out freqs1, out contracts1, out List<double> peaks1, out List<double> fitcontracts1, out exc, ps[0] - 10, ps[0] + 10, 1, 2000, 1, isReverse: false);
                if (peaks1.Count != 0)
                {
                    peaks.Add(peaks1.Min());
                }
                LabviewConverter.CW(out freqs2, out contracts2, out peaks1, out fitcontracts1, out exc, ps[0] - 10, ps[0] + 10, 1, 2000, 1, isReverse: true);
                if (peaks1.Count != 0)
                {
                    peaks.Add(peaks1.Max());
                }
                return;
            }
            if (ps.Count == 2)
            {
                peaks = new List<double>();
                double minfreq = Math.Min(ps[0], ps[1]);
                double maxfreq = Math.Max(ps[0], ps[1]);
                LabviewConverter.CW(out freqs1, out contracts1, out List<double> peaks1, out List<double> fitcontracts1, out exc, minfreq - 5, minfreq + 5, 1, 2000, 1, scanRangeAfterstop: 10, isReverse: false);
                if (peaks1.Count != 0)
                {
                    peaks.Add(peaks1.Min());
                }
                LabviewConverter.CW(out freqs2, out contracts2, out peaks1, out fitcontracts1, out exc, maxfreq - 5, maxfreq + 5, 1, 2000, 1, scanRangeAfterstop: 10, isReverse: false);
                if (peaks1.Count != 0)
                {
                    peaks.Add(peaks1.Max());
                }
                return;
            }
            peaks = new List<double>();
            freqs1 = new List<double>();
            contracts1 = new List<double>();
            freqs2 = new List<double>();
            contracts2 = new List<double>();
            exc = null;
            return;
        }

        /// <summary>
        /// 全谱扫两个共振峰，扫不全则报错
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="freqs1"></param>
        /// <param name="contracts1"></param>
        /// <param name="freqs2"></param>
        /// <param name="contracts2"></param>
        /// <exception cref="Exception"></exception>
        public static void TotalCWPeaks2OrException(out List<double> peaks, out List<double> freqs1, out List<double> contracts1, out List<double> freqs2, out List<double> contracts2)
        {
            TotalCW(out peaks, out freqs1, out contracts1, out freqs2, out contracts2);
            if (peaks.Count == 0)
            {
                //没有扫到峰，报错
                throw new Exception("未扫描到完整的谱峰");
            }
            if (peaks.Count == 1)
            {
                peaks = new List<double>() { peaks[0], peaks[0] };
                return;
            }
        }

        // 扫描出CW谱的准确峰位置
        // 先粗扫后细扫(步长先设置为2MHz,确定峰位置之后按步长0.2MHz扫描)
        // 范围限制,扫描范围在2400Mhz到3400Mhz之间
        private static void ScanCW(out List<double> peaks, out List<double> freqvalues, out List<double> contractvalues, double peakapprox, double scanWidth = 30, bool isReverse = false)
        {
            if (peakapprox < 2100 || peakapprox > 3500)
            {
                peaks = new List<double>();
                freqvalues = new List<double>();
                contractvalues = new List<double>();
                return;
            }
            peaks = new List<double>();
            LabviewConverter.CW(out freqvalues, out contractvalues, out List<double> freqs, out List<double> fitcontracts, out Exception exc, startFreq: peakapprox - scanWidth, endFreq: peakapprox + scanWidth, step: 2, averagetime: 500, peakcount: 1, isReverse: isReverse);

            double f1 = 0;
            if (freqs.Count != 0)
            {
                f1 = freqs[0];
                if (freqs.Count == 2)
                {
                    f1 = isReverse ? Math.Max(freqs[0], freqs[1]) : Math.Min(freqs[0], freqs[1]);
                }
                //如果扫到的峰在范围边缘则扩大范围扫描
                if (Math.Abs(f1 - peakapprox + scanWidth) < 3)
                {
                    LabviewConverter.CW(out freqvalues, out contractvalues, out freqs, out fitcontracts, out exc, startFreq: peakapprox - scanWidth - 20, endFreq: peakapprox - scanWidth + 20, step: 2, averagetime: 500, peakcount: 1, isReverse: isReverse);
                    if (freqs.Count != 0)
                    {
                        //细扫
                        LabviewConverter.CW(out freqvalues, out contractvalues, out List<double> f1s, out fitcontracts, out exc, startFreq: freqs[0] - 5, endFreq: freqs[0] + 5, step: 1, averagetime: 2000, peakcount: 1, scanRangeAfterstop: 10, isReverse: isReverse);
                        if (f1s.Count != 0)
                        {
                            peaks.Add(f1s[0]);
                        }
                    }
                }
                if (Math.Abs(f1 - peakapprox - scanWidth) < 3)
                {
                    LabviewConverter.CW(out freqvalues, out contractvalues, out freqs, out fitcontracts, out exc, startFreq: peakapprox + scanWidth - 20, endFreq: peakapprox + scanWidth + 20, step: 2, averagetime: 500, peakcount: 1, isReverse: isReverse);
                    if (freqs.Count != 0)
                    {
                        LabviewConverter.CW(out freqvalues, out fitcontracts, out List<double> f1s, out fitcontracts, out exc, startFreq: freqs[0] - 5, endFreq: freqs[0] + 5, step: 1, averagetime: 2000, peakcount: 1, scanRangeAfterstop: 10, isReverse: isReverse);
                        if (f1s.Count != 0)
                        {
                            peaks.Add(f1s[0]);
                        }
                    }
                }
                if (Math.Abs(f1 - peakapprox + scanWidth) >= 3 && Math.Abs(f1 - peakapprox - scanWidth) >= 3)
                {
                    LabviewConverter.CW(out freqvalues, out contractvalues, out List<double> f1s, out fitcontracts, out exc, startFreq: f1 - 5, endFreq: f1 + 5, step: 1, averagetime: 2000, peakcount: 1, scanRangeAfterstop: 10, isReverse: isReverse);
                    if (f1s.Count != 0)
                    {
                        peaks.Add(f1s[0]);
                    }
                }
            }
        }


        /// <summary>
        /// 扫两频点的CW谱
        /// </summary>
        public static void ScanCW2(out double peakout1, out double peakout2, out List<double> freqsout1, out List<double> contrastout1, out List<double> freqsout2, out List<double> contrastout2, double peakapprox1, double peakapprox2, double scanWidth = 30)
        {
            //寻找峰值较小值，正向扫描
            double p1 = Math.Min(peakapprox1, peakapprox2);
            double p2 = Math.Max(peakapprox1, peakapprox2);
            ScanCW(out List<double> peak1, out freqsout1, out contrastout1, p1, scanWidth, isReverse: false);
            ScanCW(out List<double> peak2, out freqsout2, out contrastout2, p2, scanWidth, isReverse: true);

            if (peak1.Count == 0 || peak2.Count == 0)
            {
                peakout1 = 0;
                peakout2 = 0;
                freqsout1 = new List<double>();
                freqsout2 = new List<double>();
                contrastout1 = new List<double>();
                contrastout2 = new List<double>();
                return;
            }
            peakout1 = peak1.Min();
            peakout2 = peak2.Max();
        }


        /// <summary>
        /// 用二次函数拟合
        /// </summary>
        public static List<double> FitDataWithPow2(List<double> xdata, List<double> ydata)
        {
            if (xdata.Count == 0 || ydata.Count == 0)
            {
                return new List<double>();
            }
            List<double> result = Fit.Polynomial(xdata.ToArray(), ydata.ToArray(), 2).ToList();
            return result;
        }

        public static int GetNearestDataIndex(List<double> dataarray, double targetvalue)
        {
            List<double> cal = new List<double>();
            for (int i = 0; i < dataarray.Count; i++)
            {
                cal.Add(Math.Abs(dataarray[i] - targetvalue));
            }
            return cal.IndexOf(cal.Min());
        }

        /// <summary>
        /// 三角函数拟合
        /// </summary>
        /// <returns></returns>
        public static List<double> FitSinCurve(List<double> x, List<double> y)
        {
            return Python_NetInterpretor.ExcuteFunction(Path.Combine(Environment.CurrentDirectory, "Python", "LabviewHandler", "Math.py"), "FitSinData", TimeSpan.FromMilliseconds(100000), x, y);
        }
    }
}
